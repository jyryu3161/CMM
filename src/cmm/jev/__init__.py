"""JEV agent: a decision model plays CMM to design a production strain.

CMM's existing target-discovery methods are deterministic optimisations — FSEOF, OptKnock,
MOMA/ROOM. They see the mathematics of the model and nothing else. This package adds a
judgement layer on top of them, driven by TypeSafe's JEV decision model through OpenRouter.

JEV does not generate text. It returns a typed answer chosen from criteria the caller
supplied, which is why the loop needs no output parsing and cannot be steered outside its own
vocabulary: a reaction the model does not contain is not something to reject, it is something
that cannot be said. Each answer also carries a probability for every option, so one call
ranks the whole board.

The loop is a game. CMM renders the metabolic state, JEV picks a move, CMM executes it and
re-solves, and the new flux distribution is the next frame. CMM owns the rules — a move that
makes the model infeasible or drops growth below the floor is undone whatever the agent
predicted — and JEV owns the strategy.

Entry points::

    from cmm.jev import JevConfig, run_jev_design

    config = JevConfig(model_path="e_coli_core.xml", product="EX_succ_e", rounds=5)
    result = run_jev_design(config)
    print(result.summary())

The agent needs ``OPENROUTER_API_KEY`` in the environment. Nothing else in CMM does, and
nothing else in CMM changes when it is absent.
"""

from cmm.jev._transport import (
    DecisionResult,
    JevAnswer,
    JevClient,
    JevTransportError,
    JevUsage,
    choice_question,
    noul_question,
    score_question,
)
from cmm.jev.actions import (
    ACT_ACTIONS,
    ACTION_CATALOGUE,
    LOOK_ACTIONS,
    Action,
    ActionNotApplicable,
    Intervention,
    applicable_actions,
    build_intervention,
)
from cmm.jev.engine import (
    JevConfig,
    JevResult,
    JevWorkflowError,
    RoundRecord,
    TickRecord,
    run_jev_design,
)
from cmm.jev.questions import (
    DEFAULT_QUESTION_SET,
    QUESTION_SETS,
    QuestionSet,
    get_question_set,
)
from cmm.jev.state import (
    CandidateEvidence,
    CofactorBalance,
    GameState,
    ScanCache,
    build_candidates,
    cofactor_balance,
    product_distances,
)

__all__ = [
    "ACTION_CATALOGUE",
    "ACT_ACTIONS",
    "Action",
    "ActionNotApplicable",
    "CandidateEvidence",
    "CofactorBalance",
    "DEFAULT_QUESTION_SET",
    "DecisionResult",
    "GameState",
    "Intervention",
    "JevAnswer",
    "JevClient",
    "JevConfig",
    "JevResult",
    "JevTransportError",
    "JevUsage",
    "JevWorkflowError",
    "LOOK_ACTIONS",
    "QUESTION_SETS",
    "QuestionSet",
    "RoundRecord",
    "ScanCache",
    "TickRecord",
    "applicable_actions",
    "build_candidates",
    "build_intervention",
    "choice_question",
    "cofactor_balance",
    "get_question_set",
    "noul_question",
    "product_distances",
    "run_jev_design",
    "score_question",
]
