"""Where the OpenRouter key comes from, and where it is kept if the user asks.

CMM needs a credential for exactly one feature, so this is the only place in the project that
handles one. Three rules, in order:

1. ``OPENROUTER_API_KEY`` in the environment wins. A key passed in for one session should not
   be quietly overridden by one saved months ago.
2. Otherwise a key the user saved from the desktop app, in a file only they can read.
3. Otherwise there is no key, and the JEV tab says so and disables itself.

**The saved key is stored in plain text.** There is no way around that without taking on an
OS keyring dependency and its per-platform behaviour, so the honest thing is to say so where
the user saves it, keep the file at ``0600``, and make deleting it one menu item. The file
holds nothing but the key, so a reader knows what they have found.
"""

from __future__ import annotations

import os
from pathlib import Path

ENV_VAR = "OPENROUTER_API_KEY"


def key_path() -> Path:
    """Where a saved key lives. Follows ``XDG_CONFIG_HOME`` when it is set."""

    configured = os.environ.get("XDG_CONFIG_HOME")
    root = Path(configured).expanduser() if configured else Path.home() / ".config"
    return root / "cmm" / "openrouter.key"


def stored_key() -> str:
    """The saved key, or an empty string. Never raises: a missing file is the normal case."""

    path = key_path()
    try:
        return path.read_text(encoding="utf-8").strip()
    except (OSError, UnicodeDecodeError):
        return ""


def save_key(key: str) -> Path:
    """Write the key where :func:`stored_key` will find it, readable only by this user.

    The permissions are set before the key is written, so it is never briefly world-readable
    on a filesystem that honours them.
    """

    key = key.strip()
    if not key:
        raise ValueError("refusing to save an empty API key")
    path = key_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.touch(mode=0o600, exist_ok=True)
    try:
        path.chmod(0o600)
    except OSError:  # pragma: no cover - filesystems without POSIX permissions
        pass
    path.write_text(key + "\n", encoding="utf-8")
    return path


def clear_key() -> bool:
    """Delete the saved key. Returns whether there was one to delete."""

    path = key_path()
    try:
        path.unlink()
    except FileNotFoundError:
        return False
    except OSError:  # pragma: no cover - permission or locking failure
        return False
    return True


def key_source() -> str:
    """Which of the three rules is in force, for the user interface to report."""

    if os.environ.get(ENV_VAR, "").strip():
        return "environment"
    if stored_key():
        return "saved"
    return "none"


def masked(key: str) -> str:
    """A key rendered so it can be shown on screen: the last four characters only."""

    key = key.strip()
    if len(key) <= 4:
        return "•" * len(key)
    return "•" * 8 + key[-4:]


__all__ = [
    "ENV_VAR",
    "clear_key",
    "key_path",
    "key_source",
    "masked",
    "save_key",
    "stored_key",
]
