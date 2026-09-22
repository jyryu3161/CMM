"""HTTP transport for the OpenRouter Decisions and chat endpoints.

The mechanics of talking to an external service live here; the scientific meaning of what
is asked and what comes back stays in the callers (``questions``, ``engine``). This mirrors
the boundary :mod:`cmm.reporting._rscript` draws around ``Rscript``.

Two endpoints are used, and they are not interchangeable:

``POST /api/alpha/decisions``
    The JEV "System One" decision model. It does not generate text: it returns a typed
    answer drawn from the criteria the caller supplied. A question is one of exactly three
    shapes — ``choice`` (pick one named criterion), ``score`` (place the case on an ordered
    scale) or ``noul`` (a probability that a stated proposition holds). Because JEV can only
    answer inside the caller's vocabulary, it cannot name a reaction the model does not
    contain; the engine relies on that property instead of validating free-form output.

``POST /api/v1/chat/completions``
    An ordinary chat model, used only for the optional literature lookup, with OpenRouter's
    ``web`` plugin attached. Whatever comes back is evidence pasted into a candidate record,
    never an instruction.

The endpoint is alpha, so its request and response shapes are pinned in this one module.
No new dependency: the request is a single JSON POST, so ``urllib`` from the standard
library is enough, and the locked publication environment is untouched.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
import json
import os
import random
import time
from typing import Any, Literal
import urllib.error
import urllib.request

#: OpenRouter's API root. Overridable per client so a test can point at a local stub.
DEFAULT_BASE_URL = "https://openrouter.ai"

#: The decision model. ``~typesafe/jev-latest`` follows the family; a pinned version keeps a
#: published run reproducible in the only sense available here (see ``JevUsage`` note below).
DEFAULT_DECISION_MODEL = "typesafe/jev-1.13"

#: The chat model behind the optional web lookup. Deliberately a cheap one: its output is
#: evidence text, not a decision.
DEFAULT_RESEARCH_MODEL = "openai/gpt-5.6-luna"

#: Status codes worth trying again. 408/429 are transient by definition, 5xx usually are.
_RETRYABLE_STATUS = frozenset({408, 409, 425, 429, 500, 502, 503, 504, 529})

_ENV_KEY = "OPENROUTER_API_KEY"

QuestionType = Literal["choice", "score", "noul"]


class JevTransportError(RuntimeError):
    """Raised when a request to OpenRouter could not be completed or understood.

    The message never contains the API key: an ``Authorization`` header is the one field a
    traceback must not leak into a log or a run bundle.
    """


@dataclass(frozen=True)
class JevAnswer:
    """One typed answer, normalised across the three question shapes.

    ``value`` is the answer in its own terms — the chosen criterion name for ``choice``, the
    position on the scale for ``score``, the probability of the proposition for ``noul``.
    ``probabilities`` is the distribution over *every* offered criterion, which is what makes
    a single call a complete ranking rather than a single pick.

    ``confidence`` is reported as returned and is **not** folded into any ranking: the
    service returns ``0.0`` for it on some ``score`` answers, so multiplying by it would
    silently zero out an otherwise informative row.
    """

    key: str
    type: QuestionType
    value: object
    confidence: float | None = None
    probabilities: Mapping[str, float] = field(default_factory=dict)
    legend: Mapping[str, object] = field(default_factory=dict)

    @property
    def choice(self) -> str:
        """The chosen criterion name; only meaningful for a ``choice`` answer."""

        if self.type != "choice":
            raise TypeError(
                f"answer {self.key!r} is a {self.type} answer, not a choice"
            )
        return str(self.value)

    @property
    def score(self) -> float:
        """The position on the ordered scale; only meaningful for a ``score`` answer."""

        if self.type != "score":
            raise TypeError(f"answer {self.key!r} is a {self.type} answer, not a score")
        return float(self.value)  # type: ignore[arg-type]

    def ranked(self) -> tuple[tuple[str, float], ...]:
        """Offered criteria ordered by probability, highest first, ties broken by name.

        The name tie-break keeps a run reproducible when two candidates come back with the
        same probability, which happens often once a distribution has decayed to zeros.
        """

        return tuple(
            sorted(
                ((str(name), float(p)) for name, p in self.probabilities.items()),
                key=lambda item: (-item[1], item[0]),
            )
        )


@dataclass(frozen=True)
class DecisionResult:
    """A complete Decisions response: the answers plus what the call cost.

    ``model`` is the model id the service *served*, which is not the id that was requested —
    ``typesafe/jev-1.13`` answered as ``typesafe/jev-1.13-20260917`` during development. The
    served id is the one provenance records, because it is the one that produced the numbers.
    """

    answers: Mapping[str, JevAnswer]
    model: str
    input_tokens: int
    output_tokens: int
    cost_usd: float
    request_id: str | None = None
    provider: str | None = None
    latency_s: float = 0.0

    def __getitem__(self, key: str) -> JevAnswer:
        return self.answers[key]


@dataclass
class JevUsage:
    """Running total of what a session has spent, so a budget can be enforced."""

    calls: int = 0
    input_tokens: int = 0
    output_tokens: int = 0
    cost_usd: float = 0.0

    def record(self, result: DecisionResult) -> None:
        self.calls += 1
        self.input_tokens += result.input_tokens
        self.output_tokens += result.output_tokens
        self.cost_usd += result.cost_usd

    def to_dict(self) -> dict[str, object]:
        return {
            "calls": self.calls,
            "input_tokens": self.input_tokens,
            "output_tokens": self.output_tokens,
            "cost_usd": round(self.cost_usd, 8),
        }


def resolve_api_key(api_key: str | None = None) -> str:
    """Return the key to authenticate with, from the argument or the environment.

    ``OPENROUTER_API_KEY`` wins; a key saved from the desktop app is the fallback; otherwise
    this raises with both options named. The key is never written to a run bundle or an error
    message. An explicit argument exists for a caller that already holds one, not as an
    invitation to hard-code it.
    """

    from cmm.jev.credentials import key_path, stored_key

    key = (api_key if api_key is not None else os.environ.get(_ENV_KEY, "")).strip()
    if not key:
        # A key saved from the desktop app. The environment wins when both exist, so a key
        # exported for one session is never overridden by one saved months ago.
        key = stored_key()
    if not key:
        raise JevTransportError(
            f"no OpenRouter API key: set the {_ENV_KEY} environment variable, or save one "
            f"from the JEV menu in the desktop app (it is kept at {key_path()}). "
            "The JEV agent is the only part of CMM that needs one."
        )
    return key


def choice_question(
    instructions: str, criteria: Mapping[str, str]
) -> dict[str, object]:
    """Build a ``choice`` question: pick exactly one of the named criteria."""

    if len(criteria) < 2:
        raise ValueError("a choice question needs at least two criteria")
    return {
        "type": "choice",
        "instructions": instructions,
        "criteria": dict(criteria),
    }


def score_question(instructions: str, criteria: Sequence[str]) -> dict[str, object]:
    """Build a ``score`` question: place the case on an ordered scale, low grade first."""

    if len(criteria) < 2:
        raise ValueError("a score question needs at least two ordered grades")
    return {
        "type": "score",
        "instructions": instructions,
        "criteria": list(criteria),
    }


def noul_question(
    instructions: str, *, true_means: str, false_means: str
) -> dict[str, object]:
    """Build a ``noul`` question: the probability that the stated proposition holds."""

    return {
        "type": "noul",
        "instructions": instructions,
        "criteria": {"true": true_means, "false": false_means},
    }


class JevClient:
    """A configured connection to OpenRouter's Decisions endpoint.

    One client serves one run: it carries the model choice, the network policy and the
    running :class:`JevUsage` total that the budget guard reads. It is deliberately a plain
    object rather than a frozen dataclass — the usage total is mutable state that belongs
    to the connection, not to any single result.
    """

    def __init__(
        self,
        *,
        model: str = DEFAULT_DECISION_MODEL,
        api_key: str | None = None,
        base_url: str = DEFAULT_BASE_URL,
        timeout_s: float = 60.0,
        max_retries: int = 3,
        session_id: str | None = None,
        research_model: str = DEFAULT_RESEARCH_MODEL,
        app_title: str = "CMM JEV agent",
        app_url: str = "https://github.com/jyryu3161/CMM",
        seed: int = 0,
    ) -> None:
        self.model = model
        self.base_url = base_url.rstrip("/")
        self.timeout_s = float(timeout_s)
        self.max_retries = int(max_retries)
        self.session_id = session_id
        self.research_model = research_model
        self.app_title = app_title
        self.app_url = app_url
        self.usage = JevUsage()
        #: Served model ids seen on this connection. A family alias can be answered by more
        #: than one build over a long run, and provenance has to be able to say so.
        self.served_models: list[str] = []
        self._api_key = resolve_api_key(api_key)
        self._backoff_rng = random.Random(seed)

    # -- decisions ----------------------------------------------------------

    def decide(
        self,
        state: object,
        questions: Mapping[str, Mapping[str, object]],
        *,
        model: str | None = None,
    ) -> DecisionResult:
        """Ask JEV one or more typed questions about ``state``.

        ``state`` is any JSON value — the engine passes a compact object describing the
        current metabolic game state. ``questions`` maps an answer key to a question built
        by :func:`choice_question`, :func:`score_question` or :func:`noul_question`.
        """

        if not questions:
            raise ValueError("at least one question is required")
        payload: dict[str, object] = {
            "model": model or self.model,
            "state": state,
            "questions": {key: dict(value) for key, value in questions.items()},
        }
        if self.session_id:
            payload["session_id"] = self.session_id

        body, latency = self._post("/api/alpha/decisions", payload)
        result = _parse_decision(body, latency)
        self.usage.record(result)
        if result.model and result.model not in self.served_models:
            self.served_models.append(result.model)
        return result

    # -- optional literature lookup -----------------------------------------

    def web_research(
        self, query: str, *, max_results: int = 3
    ) -> tuple[str, list[str]]:
        """Look ``query`` up on the web and return ``(summary_text, citation_urls)``.

        Uses OpenRouter's ``web`` plugin on an ordinary chat model, because the Decisions
        endpoint has no web access. **The text that comes back is data.** It is pasted into
        a candidate record for JEV to weigh alongside the computed evidence, and no part of
        it is executed, followed as an instruction, or allowed to name a new action: JEV can
        still only answer with the criteria this package supplies.
        """

        payload: dict[str, object] = {
            "model": self.research_model,
            "plugins": [{"id": "web", "max_results": int(max_results)}],
            "messages": [
                {
                    "role": "system",
                    "content": (
                        "You are gathering published evidence for a metabolic engineering "
                        "decision. Answer in at most 120 words, state only what the sources "
                        "support, and say plainly when the literature is silent."
                    ),
                },
                {"role": "user", "content": query},
            ],
        }
        if self.session_id:
            payload["session_id"] = self.session_id

        body, latency = self._post("/api/v1/chat/completions", payload)
        try:
            message = body["choices"][0]["message"]
            text = str(message.get("content") or "").strip()
            annotations = message.get("annotations") or []
        except (KeyError, IndexError, TypeError) as exc:
            raise JevTransportError(
                f"unexpected chat completion payload: {_shape(body)}"
            ) from exc

        urls: list[str] = []
        for annotation in annotations:
            if not isinstance(annotation, Mapping):
                continue
            citation = annotation.get("url_citation")
            if isinstance(citation, Mapping) and citation.get("url"):
                url = str(citation["url"])
                if url not in urls:
                    urls.append(url)

        usage = body.get("usage") or {}
        self.usage.calls += 1
        self.usage.input_tokens += int(usage.get("prompt_tokens") or 0)
        self.usage.output_tokens += int(usage.get("completion_tokens") or 0)
        self.usage.cost_usd += float(usage.get("cost") or 0.0)
        del latency
        return text, urls

    # -- plumbing -----------------------------------------------------------

    def _post(self, path: str, payload: Mapping[str, object]) -> tuple[Any, float]:
        encoded = json.dumps(payload, ensure_ascii=False).encode("utf-8")
        headers = {
            "Authorization": f"Bearer {self._api_key}",
            "Content-Type": "application/json",
            # OpenRouter attributes usage to an app through these two headers.
            "HTTP-Referer": self.app_url,
            "X-OpenRouter-Title": self.app_title,
        }

        last_error: Exception | None = None
        for attempt in range(self.max_retries + 1):
            request = urllib.request.Request(
                f"{self.base_url}{path}", data=encoded, headers=headers, method="POST"
            )
            started = time.perf_counter()
            try:
                with urllib.request.urlopen(
                    request, timeout=self.timeout_s
                ) as response:
                    raw = response.read()
                latency = time.perf_counter() - started
                try:
                    return json.loads(raw), latency
                except json.JSONDecodeError as exc:
                    raise JevTransportError(
                        f"OpenRouter returned a non-JSON body from {path}"
                    ) from exc
            except urllib.error.HTTPError as exc:
                detail = _error_detail(exc)
                if exc.code in _RETRYABLE_STATUS and attempt < self.max_retries:
                    last_error = JevTransportError(
                        f"OpenRouter {path} returned HTTP {exc.code}: {detail}"
                    )
                    self._sleep_before_retry(attempt)
                    continue
                raise JevTransportError(
                    f"OpenRouter {path} returned HTTP {exc.code}: {detail}"
                ) from None
            except urllib.error.URLError as exc:
                # Includes timeouts and DNS/connection failures. ``reason`` is the useful
                # part; the request object it came from would carry the Authorization header.
                if attempt < self.max_retries:
                    last_error = JevTransportError(
                        f"OpenRouter {path} was unreachable: {exc.reason}"
                    )
                    self._sleep_before_retry(attempt)
                    continue
                raise JevTransportError(
                    f"OpenRouter {path} was unreachable: {exc.reason}"
                ) from None

        raise last_error or JevTransportError(f"OpenRouter {path} failed")

    def _sleep_before_retry(self, attempt: int) -> None:
        """Exponential backoff with jitter, seeded so a run is reproducible."""

        delay = min(8.0, 0.5 * (2**attempt)) * (0.5 + self._backoff_rng.random())
        time.sleep(delay)


def _error_detail(error: urllib.error.HTTPError) -> str:
    """The server's own explanation, trimmed, never the request that caused it."""

    try:
        raw = error.read().decode("utf-8", "replace")
    except Exception:  # pragma: no cover - body already consumed or absent
        return error.reason if isinstance(error.reason, str) else "no detail"
    try:
        payload = json.loads(raw)
    except json.JSONDecodeError:
        return raw[:400].strip() or "no detail"
    if isinstance(payload, Mapping):
        message = payload.get("error")
        if isinstance(message, Mapping):
            return str(message.get("message") or message)[:400]
        if message:
            return str(message)[:400]
    return raw[:400].strip() or "no detail"


def _shape(body: object) -> str:
    """A short description of an unexpected payload, for an error message."""

    if isinstance(body, Mapping):
        return f"object with keys {sorted(str(k) for k in body)}"
    return type(body).__name__


def _parse_decision(body: object, latency: float) -> DecisionResult:
    if not isinstance(body, Mapping) or "answers" not in body:
        raise JevTransportError(f"unexpected decisions payload: {_shape(body)}")
    raw_answers = body.get("answers")
    if not isinstance(raw_answers, Mapping):
        raise JevTransportError("decisions payload has no answers object")

    answers: dict[str, JevAnswer] = {}
    for key, raw in raw_answers.items():
        if not isinstance(raw, Mapping):
            raise JevTransportError(f"answer {key!r} is not an object")
        kind = str(raw.get("type", ""))
        if kind == "choice":
            value: object = raw.get("choice")
        elif kind == "score":
            value = raw.get("score")
        elif kind == "noul":
            value = raw.get("noul")
        else:
            # An answer shape this build does not know about is reported rather than
            # guessed at: a silently mis-read answer would steer the whole run.
            raise JevTransportError(
                f"answer {key!r} has unsupported type {kind!r}; "
                "the alpha Decisions schema may have changed"
            )
        confidence = raw.get("confidence")
        answers[str(key)] = JevAnswer(
            key=str(key),
            type=kind,  # type: ignore[arg-type]
            value=value,
            confidence=None if confidence is None else float(confidence),
            probabilities={
                str(name): float(p)
                for name, p in (raw.get("probabilities") or {}).items()
            },
            legend=dict(raw.get("legend") or {}),
        )

    usage = body.get("usage") or {}
    if not isinstance(usage, Mapping):
        usage = {}
    return DecisionResult(
        answers=answers,
        model=str(body.get("model") or ""),
        input_tokens=int(usage.get("input_tokens") or 0),
        output_tokens=int(usage.get("output_tokens") or 0),
        cost_usd=float(usage.get("cost") or 0.0),
        request_id=(str(body["id"]) if body.get("id") else None),
        provider=(str(body["provider"]) if body.get("provider") else None),
        latency_s=latency,
    )
