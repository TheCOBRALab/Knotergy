import src.settings as s
import src.functions as f

def main():
    parser = s.make_parser()
    args = parser.parse_args()

    # If user did not provide --input or --sequence, ask for sequence
    if not args.input and not args.sequence:
        args.sequence = input("Enter RNA sequence: ").strip().upper().replace("T", "U")

        if not args.sequence:
            parser.error("No input file or sequence provided.")

    if args.max_layers > len(s.BRACKETS):
        raise ValueError(f"--max-layers cannot exceed {len(s.BRACKETS)}")

    entries = f.entries_from_args(args)
    predictions = []

    for rec in f.tqdm(entries, desc="Input records", unit="record"):
        pred = f.predict(rec, args)
        predictions.append(pred)

        best = pred["best"]

        print("\nBest result")
        print("Name:     ", rec.name)

        if best is None:
            print("No result")
            continue

        print("Energy:   ", best["energy"])
        print("Pairs:    ", best["pair_count"])
        print("Density:  ", f"{best['density']:.3f}")
        print("Layers:   ", best["layers"])
        print("Source:   ", best["source"])
        print("Sequence: ", rec.seq)
        print("Structure:", best["structure"])

        print(f"\nTop {args.top} Knotergy candidates:")
        for r in pred["results"][:args.top]:
            print(
                f"{r['energy']:10.4f}  "
                f"pairs={r['pair_count']:4d}  "
                f"density={r['density']:.3f}  "
                f"layers={r['layers']}  "
                f"{r['source']:10s}  "
                f"{r['structure']}"
            )

    if args.output:
        f.write_output(args.output, predictions)
        print(f"\nWrote output to {args.output}")


if __name__ == "__main__":
    main()