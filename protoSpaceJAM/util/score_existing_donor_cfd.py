import argparse
import csv
import re
from pathlib import Path

try:
    from protoSpaceJAM.util.cfdscore import cfd_score
except ImportError:  # pragma: no cover
    from cfdscore import cfd_score


DNA_RE = re.compile(r"[^ACGTNacgtn]")


def reverse_complement(seq):
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def normalize_dna(seq):
    return DNA_RE.sub("", seq).upper()


def load_sequence(sequence=None, sequence_file=None):
    if bool(sequence) == bool(sequence_file):
        raise ValueError("Provide exactly one of --sequence or --sequence-file.")
    if sequence:
        return normalize_dna(sequence)
    text = Path(sequence_file).read_text()
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        raise ValueError(f"No sequence found in {sequence_file}")
    if lines[0].startswith(">"):
        lines = [line for line in lines if not line.startswith(">")]
    return normalize_dna("".join(lines))


def split_guide_and_pam(guide_seq, pam):
    guide_seq = normalize_dna(guide_seq)
    pam = normalize_dna(pam)
    if not guide_seq:
        raise ValueError("Guide sequence is empty after normalization.")
    if len(guide_seq) == 23:
        return guide_seq[:20], guide_seq[20:]
    if len(guide_seq) == 20:
        if not pam:
            raise ValueError("Provide --pam when using a 20 nt guide sequence.")
        return guide_seq, pam
    raise ValueError(
        "Guide sequence must be either 20 nt protospacer-only or 23 nt protospacer+PAM."
    )


def iter_scored_windows(search_seq, guide_seq, pam, sequence_name):
    protospacer, pam = split_guide_and_pam(guide_seq, pam)
    window_size = len(protospacer) + len(pam)
    if len(search_seq) < window_size:
        raise ValueError(
            f"Sequence is shorter than the required scoring window ({window_size} bp)."
        )

    for strand, seq in (("+", search_seq), ("-", reverse_complement(search_seq))):
        for idx in range(0, len(seq) - window_size + 1):
            window = seq[idx : idx + window_size]
            if "N" in window:
                continue
            score = cfd_score(protospacer, window, pam=pam)
            if strand == "+":
                start_1based = idx + 1
                end_1based = idx + window_size
            else:
                start_1based = len(search_seq) - (idx + window_size) + 1
                end_1based = len(search_seq) - idx
            yield {
                "sequence_name": sequence_name,
                "guide_seq": protospacer,
                "pam": pam,
                "strand": strand,
                "window_start_1based": start_1based,
                "window_end_1based": end_1based,
                "window_seq": window,
                "cfd_score": score,
            }


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Score an existing donor or homology-arm-containing sequence against a guide "
            "using protoSpaceJAM's CFD scoring function."
        )
    )
    parser.add_argument("--guide-seq", required=True, help="20 nt protospacer or 23 nt protospacer+PAM.")
    parser.add_argument("--pam", default="NGG", help="PAM to use when --guide-seq is 20 nt. Default: NGG")
    parser.add_argument("--sequence", help="Raw donor / HA sequence to scan.")
    parser.add_argument("--sequence-file", help="Path to a raw-text or FASTA sequence file to scan.")
    parser.add_argument("--sequence-name", default="input_sequence", help="Label to use in output.")
    parser.add_argument("--top-n", type=int, default=10, help="Number of highest-scoring windows to report. Default: 10")
    parser.add_argument("--output-csv", help="Optional CSV path for the top-scoring windows.")
    args = parser.parse_args()

    search_seq = load_sequence(sequence=args.sequence, sequence_file=args.sequence_file)
    rows = sorted(
        iter_scored_windows(search_seq, args.guide_seq, args.pam, args.sequence_name),
        key=lambda row: row["cfd_score"],
        reverse=True,
    )

    if not rows:
        raise ValueError("No valid windows were scored.")

    top_rows = rows[: max(1, args.top_n)]
    top_hit = top_rows[0]

    print(f"Sequence: {args.sequence_name}")
    print(f"Guide: {top_hit['guide_seq']}")
    print(f"PAM: {top_hit['pam']}")
    print(f"Top CFD: {top_hit['cfd_score']:.6f}")
    print(
        "Top window: "
        f"{top_hit['window_seq']} "
        f"({top_hit['strand']} strand, "
        f"{top_hit['window_start_1based']}-{top_hit['window_end_1based']})"
    )
    print(f"Windows scored: {len(rows)}")

    if args.output_csv:
        fieldnames = [
            "sequence_name",
            "guide_seq",
            "pam",
            "strand",
            "window_start_1based",
            "window_end_1based",
            "window_seq",
            "cfd_score",
        ]
        with open(args.output_csv, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(top_rows)
        print(f"Wrote top {len(top_rows)} hits to {args.output_csv}")


if __name__ == "__main__":
    main()
