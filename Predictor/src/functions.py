import argparse
import math
import os
import random
import re
import subprocess
import tempfile

from dataclasses import dataclass

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(x=None, **kwargs):
        return x if x is not None else range(kwargs.get("total", 0))

# Local imports
import src.data_structures as ds
import src.settings as s

def norm(seq):
    return seq.upper().replace("T", "U")

def can_pair(seq, i, j, min_loop):
    return i < j and j - i - 1 >= min_loop and seq[j] in ds.VALID_PAIRS.get(seq[i], set())

def crosses(a, b):
    i, j = a
    k, l = b
    return i < k < j < l or k < i < l < j

# Stores all positions used in pairs into a frozenset
def used_positions(pairs):
    return frozenset(x for p in pairs for x in p)

# Parses a structure string into a list of pairs
def parse_pairs(struct):
    stacks = {op: [] for op in ds.OPEN_TO_CLOSE}
    pairs = []

    for i, ch in enumerate(struct):
        if ch in ds.OPEN_TO_CLOSE:
            stacks[ch].append(i)
        elif ch in ds.CLOSE_TO_OPEN:
            op = ds.CLOSE_TO_OPEN[ch]
            if not stacks[op]:
                raise ValueError(f"Unbalanced structure near position {i}")
            pairs.append((stacks[op].pop(), i))
        elif ch != ".":
            raise ValueError(f"Invalid structure character {ch!r} at position {i}")

    leftovers = [op for op, st in stacks.items() if st]
    if leftovers:
        raise ValueError(f"Unbalanced structure: leftover opener(s) {leftovers}")

    return tuple(sorted(pairs))

# Returns a list of layers, where each layer is a list of pairs that do not cross each other.
def layers_for(pairs, max_layers):
    layers = []

    for p in sorted(pairs):
        for layer in layers:
            if not any(crosses(p, q) for q in layer):
                layer.append(p)
                break
        else:
            layers.append([p])
            if len(layers) > max_layers:
                return None

    return layers

# Converts a list of pairs into a structure string, using different brackets for different layers.
def pairs_to_struct(n, pairs, max_layers):
    layers = layers_for(pairs, max_layers)
    if layers is None:
        return None

    s = ["."] * n
    for layer_id, layer in enumerate(layers):
        left, right = ds.BRACKETS[layer_id]
        for i, j in layer:
            s[i] = left
            s[j] = right

    return "".join(s)

# Computes statistics about a structure, including the number of pairs, density, and number of layers.
def stats(struct, max_layers):
    pairs = parse_pairs(struct)
    layers = layers_for(pairs, max_layers)
    return {
        "pairs": len(pairs),
        "density": 2 * len(pairs) / len(struct),
        "layers": None if layers is None else len(layers),
    }


def read_fasta(path):
    with open(path) as f:
        lines = [x.strip() for x in f if x.strip()]

    entries = []
    i = 0

    while i < len(lines):
        if not lines[i].startswith(">"):
            raise ValueError(f"Expected header starting with '>', got {lines[i]!r}")

        name = lines[i][1:].strip()
        seq = norm(lines[i + 1])
        struct = lines[i + 2].strip()
        i += 3

        if len(seq) != len(struct):
            raise ValueError(f"{name}: sequence length {len(seq)} != structure length {len(struct)}")

        parse_pairs(struct)
        entries.append(ds.RNAentry(name, seq, struct))

    return entries

def entries_from_args(args):
    if args.input:
        return read_fasta(args.input)

    seq = norm(args.sequence)
    struct = args.structure.strip() if args.structure else "." * len(seq)

    if len(seq) != len(struct):
        raise ValueError(f"sequence length {len(seq)} != structure length {len(struct)}")

    parse_pairs(struct)
    return [ds.RNAentry(args.name, seq, struct)]

# Computes the score of a stem based on the sequence, pairs, and scoring parameters.
def stem_score(seq, pairs, args):
    if args.objective == "density":
        return 100.0 * len(pairs) + args.stack_bonus * (len(pairs) - 1)

    score = -args.stem_init_penalty
    score += sum(s.PAIR_SCORE.get((seq[i], seq[j]), 0.0) for i, j in pairs)
    score += args.stack_bonus * (len(pairs) - 1)
    score += args.long_stem_bonus * max(0, len(pairs) - 5)
    return score

# Computes the list of pairs in a stem starting from positions i and j, given the sequence and minimum loop size.
def stem_from(seq, i, j, args):
    pairs = []
    while i < j and can_pair(seq, i, j, args.min_loop_size):
        pairs.append((i, j))
        i += 1
        j -= 1
    return pairs

# Enumerates all possible stems in the sequence, scoring them and returning the top stems based on the scoring criteria.
def enumerate_stems(seq, args):
    stems = {}
    n = len(seq)

    for i in tqdm(range(n), desc="Enumerating stems", unit="nt"):
        for j in range(i + args.min_loop_size + 1, n):
            full_stem = stem_from(seq, i, j, args)
            if len(full_stem) < args.min_stem_length:
                continue

            for L in range(args.min_stem_length, min(len(full_stem), args.max_stem_length) + 1):
                pairs = tuple(full_stem[:L])
                score = stem_score(seq, pairs, args)
                if score <= 0:
                    continue

                stem = ds.Stem(pairs, used_positions(pairs), score, L)
                if pairs not in stems or score > stems[pairs].score:
                    stems[pairs] = stem

    out = sorted(stems.values(), key=lambda s: (s.score, s.length), reverse=True)
    return out[:args.top_stems]

# Computes the penalty for crossing stems based on the selected stems and the scoring parameters.
def crossing_penalty(stems, indices, args):
    chosen = [stems[i] for i in indices]
    penalty = 0.0

    for a in range(len(chosen)):
        for b in range(a + 1, len(chosen)):
            if any(crosses(p, q) for p in chosen[a].pairs for q in chosen[b].pairs):
                penalty += args.crossing_stem_penalty

    return penalty

# Creates a state representation of the selected stems, including the used positions, pairs, score, and structure.
def make_state(seq, stems, indices, args):
    indices = frozenset(indices)
    used = set()
    pairs = []
    score = 0.0

    for idx in indices:
        stem = stems[idx]
        if used & stem.used:
            return None
        used |= stem.used
        pairs.extend(stem.pairs)
        score += stem.score

    pairs = tuple(sorted(pairs))
    struct = pairs_to_struct(len(seq), pairs, args.max_layers)
    if struct is None:
        return None

    score -= crossing_penalty(stems, indices, args)

    return {
        "indices": indices,
        "used": frozenset(used),
        "pairs": pairs,
        "score": score,
        "structure": struct,
    }

def empty_state():
    return {
        "indices": frozenset(),
        "used": frozenset(),
        "pairs": tuple(),
        "score": 0.0,
        "structure": None,
    }

# Attempts to add a random stem to the current state, returning a new state if successful, or the original state if no valid addition is found.
def add_move(seq, stems, state, args):
    current = set(state["indices"])

    for _ in range(args.random_add_attempts):
        idx = random.randrange(len(stems))
        if idx in current or state["used"] & stems[idx].used:
            continue

        trial = make_state(seq, stems, current | {idx}, args)
        if trial:
            return trial

    return state

# Attempts to remove a random stem from the current state, returning a new state with one less stem.
def remove_move(seq, stems, state, args):
    if not state["indices"]:
        return state

    current = set(state["indices"])
    current.remove(random.choice(tuple(current)))
    return make_state(seq, stems, current, args) or empty_state()

# Attempts to swap a random stem in the current state with another stem, returning a new state if successful, or the original state if no valid swap is found.
def swap_move(seq, stems, state, args):
    if not state["indices"]:
        return add_move(seq, stems, state, args)

    current = set(state["indices"])
    current.remove(random.choice(tuple(current)))
    base = make_state(seq, stems, current, args) or empty_state()
    return add_move(seq, stems, base, args)

# Proposes a move (add, remove, or swap) based on the specified probabilities, returning a new state after the move.
def propose(seq, stems, state, args):
    r = random.random()

    if r < args.add_probability:
        return add_move(seq, stems, state, args)

    if r < args.add_probability + args.remove_probability:
        return remove_move(seq, stems, state, args)

    return swap_move(seq, stems, state, args)


# Performs a stochastic search for the best RNA structure by repeatedly proposing moves and accepting them based on their score and a temperature parameter, returning the top candidate structures found.
def stochastic_search(seq, stems, args):
    survivors = {}

    for _ in tqdm(range(args.restarts), desc="Stochastic restarts", unit="restart"):
        state = empty_state()
        temp = args.start_temp

        for step in range(args.steps_per_restart):
            trial = propose(seq, stems, state, args)
            delta = trial["score"] - state["score"]

            if delta >= 0 or random.random() < math.exp(delta / max(temp, 1e-9)):
                state = trial

            temp *= args.cooling

            if step % args.sample_every == 0 and state["structure"]:
                old = survivors.get(state["structure"])
                if old is None or state["score"] > old["score"]:
                    survivors[state["structure"]] = state

        if state["structure"]:
            old = survivors.get(state["structure"])
            if old is None or state["score"] > old["score"]:
                survivors[state["structure"]] = state

    candidates = sorted(
        survivors.values(),
        key=lambda s: (s["score"], len(s["pairs"])),
        reverse=True,
    )

    return candidates[:args.max_candidates]

# Writes a batch of RNA entries to a file in FASTA format, for knotergy to process
def write_batch(entries, path):
    with open(path, "w") as f:
        for name, seq, struct in entries:
            f.write(f">{name}\n{seq}\n{struct}\n")


# Runs knotergy on a batch of RNA entries, capturing the output energies and returning them as a list. If knotergy fails, it returns None for the corresponding entries.
def run_knotergy_batch(entries, args):
    energies = []

    for start in tqdm(range(0, len(entries), args.batch_size), desc="Knotergy batches", unit="batch"):
        chunk = entries[start:start + args.batch_size]

        with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as tmp:
            path = tmp.name

        try:
            write_batch(chunk, path)

            cmd = [args.knotergy, "-i", path]
            if args.param_file and args.param_file.lower() != "none":
                cmd += ["-P", args.param_file]
            if args.round:
                cmd += ["-e"]

            result = subprocess.run(cmd, capture_output=True, text=True)
            output = result.stdout + "\n" + result.stderr

            if result.returncode != 0:
                print("\nKnotergy batch failed:")
                print(output)
                energies.extend([None] * len(chunk))
                continue

            found = [float(x) for x in re.findall(r"ENERGY:\s*([-+]?\d*\.?\d+)\s*kcal/mol", output)]

            if len(found) != len(chunk):
                print(f"\nWarning: expected {len(chunk)} energies, got {len(found)}")
                found += [None] * (len(chunk) - len(found))
                found = found[:len(chunk)]

            energies.extend(found)

        finally:
            try:
                os.remove(path)
            except OSError:
                print(f"Error occurred while trying to remove temporary file: {path}")

    return energies

# Compares two candidate states based on the specified objective (density or energy), returning True if the first state is better than the second.
def better(a, b, args):
    if b is None:
        return True

    if args.objective == "density":
        return (
            a["pair_count"] > b["pair_count"]
            or (
                a["pair_count"] == b["pair_count"]
                and a["energy"] < b["energy"]
            )
        )

    return a["energy"] < b["energy"]

# Creates a result object containing the source, state, structure, energy, score, pair count, density, and number of layers for a given RNA entry.
def result_obj(source, state, struct, energy, args):
    st = stats(struct, args.max_layers)

    return {
        "energy": energy,
        "source": source,
        "state": state,
        "structure": struct,
        "score": None if state is None else state["score"],
        "pair_count": st["pairs"],
        "density": st["density"],
        "layers": st["layers"],
    }

# Evaluates a record against a list of candidate states, running knotergy to compute energies and returning a sorted list of results based on the specified objective (density or energy).
def evaluate(record, candidates, args):
    entries = [(f"{record.name}_input", record.seq, record.struct)]
    items = [("input", None, record.struct)]
    seen = {record.struct}

    for i, state in enumerate(candidates):
        struct = state["structure"]
        if struct in seen:
            continue
        seen.add(struct)
        entries.append((f"{record.name}_cand_{i}", record.seq, struct))
        items.append(("generated", state, struct))

    energies = run_knotergy_batch(entries, args)
    results = []

    for (source, state, struct), energy in zip(items, energies):
        if energy is not None:
            results.append(result_obj(source, state, struct, energy, args))

    results.sort(
        key=(
            (lambda x: (-x["pair_count"], x["energy"]))
            if args.objective == "density"
            else (lambda x: x["energy"])
        )
    )

    return results

# Refines the best candidate state by iteratively adding, removing, or swapping stems, running knotergy to evaluate energies, and returning the best refined result.
def refine(record, stems, best, args):
    if best["state"] is None:
        return best

    current = best["state"]
    current_best = best

    for pass_id in tqdm(range(args.refine_passes), desc="Refinement passes", unit="pass"):
        selected = set(current["indices"])
        trials = []

        # remove one
        for idx in list(selected):
            s = make_state(record.seq, stems, selected - {idx}, args)
            if s and s["pairs"]:
                trials.append(s)

        # add one
        add_count = 0
        for idx, stem in enumerate(stems):
            if idx in selected or current["used"] & stem.used:
                continue

            s = make_state(record.seq, stems, selected | {idx}, args)
            if s:
                trials.append(s)
                add_count += 1

            if add_count >= args.refine_add_limit:
                break

        # swap one
        swap_count = 0
        for remove_idx in list(selected):
            base = selected - {remove_idx}

            for add_idx in range(len(stems)):
                if add_idx in base:
                    continue

                s = make_state(record.seq, stems, base | {add_idx}, args)
                if s:
                    trials.append(s)
                    swap_count += 1

                if swap_count >= args.refine_swap_limit:
                    break

            if swap_count >= args.refine_swap_limit:
                break

        if not trials:
            break

        entries = []
        by_struct = {}
        seen = set()

        for i, s in enumerate(trials):
            struct = s["structure"]
            if struct in seen:
                continue
            seen.add(struct)
            by_struct[struct] = s
            entries.append((f"{record.name}_refine_{pass_id}_{i}", record.seq, struct))

        energies = run_knotergy_batch(entries, args)

        new_best = current_best

        for (_, _, struct), energy in zip(entries, energies):
            if energy is None:
                continue

            candidate = result_obj("refined", by_struct[struct], struct, energy, args)

            if better(candidate, new_best, args):
                new_best = candidate

        if new_best is current_best:
            break

        current_best = new_best
        current = new_best["state"]
        print(f"Refined to {new_best['energy']:.4f}")

    return current_best

# Predicts the best RNA structure for a given record, enumerating stems, performing stochastic search, 
# evaluating candidates, and refining the best result. 
# Returns a dictionary containing the record, the best result, and all evaluated results.
def predict(record, args):
    print("\n" + "=" * 80)
    print(f"Record: {record.name}")
    print(f"Length: {len(record.seq)}")
    print("=" * 80)

    random.seed(args.seed)

    stems = enumerate_stems(record.seq, args)
    print(f"Stems kept: {len(stems)}")

    if not stems:
        return {"record": record, "best": None, "results": []}

    candidates = stochastic_search(record.seq, stems, args)
    print(f"Generated candidate structures: {len(candidates)}")

    results = evaluate(record, candidates, args)

    if not results:
        return {"record": record, "best": None, "results": []}

    best = results[0]

    print("\nBest before refinement:")
    print(f"Energy:    {best['energy']}")
    print(f"Pairs:     {best['pair_count']}")
    print(f"Density:   {best['density']:.3f}")
    print(f"Layers:    {best['layers']}")
    print(f"Source:    {best['source']}")
    print(f"Structure: {best['structure']}")

    if args.refine_passes > 0:
        refined = refine(record, stems, best, args)
        if better(refined, best, args):
            best = refined

    return {"record": record, "best": best, "results": results}

# Writes the final predictions to an output file in FASTA format, including the best structure and its energy for each RNA entry.
def write_output(path, predictions):
    with open(path, "w") as f:
        for pred in predictions:
            rec = pred["record"]
            best = pred["best"]
            energy = None if best is None else best["energy"]
            struct = rec.struct if best is None else best["structure"]

            f.write(f">{rec.name} ENERGY={energy}\n{rec.seq}\n{struct}\n")


