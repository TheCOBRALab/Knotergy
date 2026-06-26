import argparse

from src.data_structures import BRACKETS

SEED = None

PAIR_SCORE = {
    ("G", "C"): 4.0,
    ("C", "G"): 4.0,
    ("A", "U"): 2.5,
    ("U", "A"): 2.5,
    ("G", "U"): 1.0,
    ("U", "G"): 1.0,
}

def make_parser():
    p = argparse.ArgumentParser(description="Condensed stochastic pseudoknot predictor using Knotergy batch scoring.")

    src = p.add_mutually_exclusive_group(required=False)
    src.add_argument("-i", "--input", help="Input file: >name / sequence / structure")
    src.add_argument("-s", "--sequence", help="Raw sequence. T is converted to U.")

    p.add_argument("-r", "--structure", help="Optional starting structure for --sequence")
    p.add_argument("--name", default="sequence_0")
    p.add_argument("-o", "--output")
    p.add_argument("-e", "--round", action="store_true")
    p.add_argument("--knotergy", default="Knotergy")
    p.add_argument("--param-file", default="none", help="Use 'none' to omit -p")

    p.add_argument("--objective", choices=["energy", "density"], default="energy")
    p.add_argument("--seed", type=int, default=42)

    p.add_argument("--min-loop-size", type=int, default=3)
    p.add_argument("--min-stem-length", type=int, default=2)
    p.add_argument("--max-stem-length", type=int, default=40)
    p.add_argument("--top-stems", type=int, default=2500)
    p.add_argument("--max-layers", type=int, default=len(BRACKETS))

    p.add_argument("--restarts", type=int, default=60)
    p.add_argument("--steps-per-restart", type=int, default=500)
    p.add_argument("--sample-every", type=int, default=20)
    p.add_argument("--max-candidates", type=int, default=300)
    p.add_argument("--batch-size", type=int, default=500)

    p.add_argument("--refine-passes", type=int, default=2)
    p.add_argument("--refine-add-limit", type=int, default=300)
    p.add_argument("--refine-swap-limit", type=int, default=300)

    p.add_argument("--random-add-attempts", type=int, default=80)
    p.add_argument("--stem-init-penalty", type=float, default=7.0)
    p.add_argument("--stack-bonus", type=float, default=3.0)
    p.add_argument("--long-stem-bonus", type=float, default=0.4)
    p.add_argument("--crossing-stem-penalty", type=float, default=0.35)

    p.add_argument("--start-temp", type=float, default=12.0)
    p.add_argument("--cooling", type=float, default=0.997)
    p.add_argument("--add-probability", type=float, default=0.45)
    p.add_argument("--remove-probability", type=float, default=0.25)

    p.add_argument("--top", type=int, default=10)

    return p