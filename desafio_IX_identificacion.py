#!/usr/bin/env python3
"""
Desafío IX - Identificar la proteína problema (262 aa) por BLASTp vs nr.

Consigna: identificar función, grupo taxonómico y cuántas secuencias
similares hay a un nivel de significancia adecuado.

Uso offline:
    python desafio_IX_identificacion.py --xml blast_results/IX_seq2_nr.xml --out-dir resultados_IX

Para correr BLAST nuevo:
    python desafio_IX_identificacion.py --run-blast --hits 1000 --expect 10 --out-dir resultados_IX

Criterio de significancia: E-value <= 1e-10 (estricto) y también se
reporta <= 1e-5 como referencia. Ver informe para justificación.
"""

import argparse
import sys
from pathlib import Path

SEQ2 = "MIDKSAFVHPTAIVEEGASIGANAHIGPFCIVGPHVEIGEGTVLKSHVVVNGHTKIGRDNEIYQFASIGEVNQDLKYAGEPTRVEIGDRNRIRESVTIHRGTVQGGGLTKVGSDNLLMINAHIAHDCTVGNRCILANNATLAGHVSVDDFAIIGGMTAVHQFCIIGAHVMVGGCSGVAQDVPPYVIAQGNHATPFGVNIEGLKRRGFSREAITAIRNAYKLIYRSGKTLDEVKPEIAELAETYPEVKAFTDFFARSTRGLIR"

try:
    from blast_common import parse_blast_xml
except ImportError:
    sys.path.insert(0, str(Path(__file__).parent))
    from blast_common import parse_blast_xml

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description="Desafío IX: identificar proteína")
    ap.add_argument("--xml", help="XML BLASTp vs nr")
    ap.add_argument("--out-dir", default="resultados_IX")
    ap.add_argument("--run-blast", action="store_true")
    ap.add_argument("--hits", type=int, default=1000)
    ap.add_argument("--expect", type=float, default=10)
    ap.add_argument("--evalue-corte", type=float, default=1e-10)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.run_blast:
        from Bio.Blast import NCBIWWW
        print(f"BLASTp nr len={len(SEQ2)} expect={args.expect} hits={args.hits} ...", flush=True)
        h = NCBIWWW.qblast("blastp", "nr", SEQ2, expect=args.expect, hitlist_size=args.hits, format_type="XML")
        (out_dir / "IX_seq2_nr.xml").write_text(h.read())
        print("Guardado.", flush=True)
        return

    if not args.xml:
        ap.print_help()
        sys.exit(1)

    df = parse_blast_xml(args.xml, len(SEQ2))
    print(f"Hits totales: {len(df)}")
    if df.empty:
        print("Sin hits.", file=sys.stderr)
        sys.exit(1)

    df.to_csv(out_dir / "IX_todos.tsv", sep="\t", index=False)
    top20 = df.head(20)
    top20.to_csv(out_dir / "IX_top20.tsv", sep="\t", index=False)
    print("\nTop 10:")
    print(df[["accession", "identidad_%", "cobertura_%", "evalue", "bit_score"]].head(10).to_string(index=False))

    for corte in [1e-5, 1e-10]:
        n = int((df["evalue"] <= corte).sum())
        print(f"E-value <= {corte:g}: {n} secuencias similares")
    # Cobertura alta + significativo
    sig = df[df["evalue"] <= args.evalue_corte]
    print(f"\nCon E<={args.evalue_corte:g}: {len(sig)} hits. Rango identidad: "
          f"{sig['identidad_%'].min():.1f}-{sig['identidad_%'].max():.1f}% " if len(sig) else "sin hits significativos.")
    sig.to_csv(out_dir / f"IX_significativos_E{args.evalue_corte:g}.tsv", sep="\t", index=False)
    print(f"Tablas en {out_dir}/")


if __name__ == "__main__":
    main()
