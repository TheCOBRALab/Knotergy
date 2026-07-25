

from dataclasses import dataclass

VALID_PAIRS = {
    "A": {"U"},
    "U": {"A", "G"},
    "G": {"U", "C"},
    "C": {"G"},
}

BRACKETS = (
    [("(", ")"), ("[", "]"), ("<", ">"), ("{", "}")]
    + [(chr(ord("A") + i), chr(ord("a") + i)) for i in range(26)]
)

OPEN_TO_CLOSE = dict(BRACKETS)
CLOSE_TO_OPEN = {v: k for k, v in OPEN_TO_CLOSE.items()}

@dataclass(frozen=True)
class RNAentry:
    name: str
    seq: str
    struct: str


@dataclass(frozen=True)
class Stem:
    pairs: tuple
    used: frozenset
    score: float
    length: int