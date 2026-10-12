#!/usr/bin/env python3
"""
BLASTp local con NCBI BLAST+ contra Swiss-Prot (fallback reproducible).

Se usa cuando el BLAST web vs nr queda en cola o el XML de 20000 hits
es inmanejable. Misma matriz (BLOSUM62), mismo Expect=100 y mismo filtro
de cobertura 70%. La tendencia E-value vs %id / Score vs %id es la misma,
solo cambia el N total (Swiss-Prot << nr).

Requiere: brew install blast + base swissprot descargada de
https://ftp.ncbi.nlm.nih.gov/blast/db/ (o generada con makeblastdb
desde UniProt Swiss-Prot FASTA).

Uso:
    python blast_local.py --db /tmp/blastdb/swissprot --out resultados_local
"""

import argparse
import subprocess
import sys
from pathlib import Path

SEQ_FULL = "VVGGLGGYMLGSAMSRPIIHFGSDYEDRYYRENMHRYPNQVYYRPMDEYSNQNNFVHDCVNITIKQHTVTTTTKGENFTETDVKMMERVVEQMCITQYERESQAYYQRGSSMVLFSSPPVILLISFLIFLIVG"
SEQ_HALF = SEQ_FULL[:len(SEQ_FULL)//2]
SEQ_QUARTER = SEQ_FULL[:len(SEQ_FULL)//4]
SEQ2 = "MIDKSAFVHPTAIVEEGASIGANAHIGPFCIVGPHVEIGEGTVLKSHVVVNGHTKIGRDNEIYQFASIGEVNQDLKYAGEPTRVEIGDRNRIRESVTIHRGTVQGGGLTKVGSDNLLMINAHIAHDCTVGNRCILANNATLAGHVSVDDFAIIGGMTAVHQFCIIGAHVMVGGCSGVAQDVPPYVIAQGNHATPFGVNIEGLKRRGFSREAITAIRNAYKLIYRSGKTLDEVKPEIAELAETYPEVKAFTDFFARSTRGLIR"


def run_blastp(db: str, seq: str, evalue: float, max_targets: int, out: Path, num_threads: int = 4):
    q = out.with_suffix(".query.fasta")
    q.write_text(f">query\n{seq}\n")
    cmd = ["blastp", "-query", str(q), "-db", db, "-evalue", str(evalue),
           "-max_target_seqs", str(max_targets), "-outfmt", "5",
           "-out", str(out), "-num_threads", str(num_threads)]
    print(f"Corriendo: {' '.join(cmd)}", flush=True)
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stderr[:2000], file=sys.stderr)
        raise RuntimeError(f"blastp falló (rc={r.returncode})")
    print(f"OK {out} ({out.stat().st_size} bytes)", flush=True)


def main():
    ap = argparse.ArgumentParser(description="BLAST local vs Swiss-Prot")
    ap.add_argument("--db", required=True, help="Prefijo base BLAST (ej. /tmp/blastdb/swissprot)")
    ap.add_argument("--out", default="resultados_local", help="Carpeta salida XML")
    ap.add_argument("--evalue", type=float, default=100)
    ap.add_argument("--hits", type=int, default=5000)
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    run_blastp(args.db, SEQ_FULL, args.evalue, args.hits, out / "VII_full_swissprot.xml")
    run_blastp(args.db, SEQ_HALF, args.evalue, 2000, out / "VIII_half_swissprot.xml")
    run_blastp(args.db, SEQ_QUARTER, args.evalue, 2000, out / "VIII_quarter_swissprot.xml")
    run_blastp(args.db, SEQ2, 10, 1000, out / "IX_seq2_swissprot.xml")
    print(f"Listo en {out}/")


if __name__ == "__main__":
    main()
