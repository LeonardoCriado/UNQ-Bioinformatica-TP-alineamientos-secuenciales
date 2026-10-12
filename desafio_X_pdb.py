#!/usr/bin/env python3
"""
Desafío X - BLASTp de la misma secuencia (262 aa) contra PDB.

Consigna: comparar con el resultado vs nr. ¿Mismos hits? ¿Qué tipo
de hits se recuperan? ¿Cuándo conviene usar PDB?

Uso offline:
    python desafio_X_pdb.py --xml-nr blast_results/IX_seq2_nr.xml --xml-pdb blast_results/X_seq2_pdb.xml --out-dir resultados_X

Para correr BLAST nuevo vs PDB:
    python desafio_X_pdb.py --run-blast --hits 500 --expect 10 --out-dir resultados_X
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
    ap = argparse.ArgumentParser(description="Desafío X: nr vs PDB")
    ap.add_argument("--xml-nr", help="XML vs nr (del desafío IX)")
    ap.add_argument("--xml-pdb", help="XML vs pdb")
    ap.add_argument("--out-dir", default="resultados_X")
    ap.add_argument("--run-blast", action="store_true")
    ap.add_argument("--hits", type=int, default=500)
    ap.add_argument("--expect", type=float, default=10)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.run_blast:
        from Bio.Blast import NCBIWWW
        print(f"BLASTp pdb len={len(SEQ2)} ...", flush=True)
        h = NCBIWWW.qblast("blastp", "pdb", SEQ2, expect=args.expect, hitlist_size=args.hits, format_type="XML")
        (out_dir / "X_seq2_pdb.xml").write_text(h.read())
        print("Guardado.", flush=True)
        return

    if not (args.xml_nr and args.xml_pdb):
        ap.print_help()
        sys.exit(1)

    df_nr = parse_blast_xml(args.xml_nr, len(SEQ2))
    df_pdb = parse_blast_xml(args.xml_pdb, len(SEQ2))
    print(f"nr: {len(df_nr)} hits | pdb: {len(df_pdb)} hits")
    df_pdb.to_csv(out_dir / "X_pdb_todos.tsv", sep="\t", index=False)
    df_pdb.head(20).to_csv(out_dir / "X_pdb_top20.tsv", sep="\t", index=False)
    print("\nTop 10 PDB:")
    if len(df_pdb):
        print(df_pdb[["accession", "identidad_%", "cobertura_%", "evalue", "bit_score"]].head(10).to_string(index=False))
        print(f"\nPDB E<=1e-10: {int((df_pdb['evalue'] <= 1e-10).sum())} | E<=1e-5: {int((df_pdb['evalue'] <= 1e-5).sum())}")
    else:
        print("(sin hits en PDB con esos parámetros)")

    # ¿Accesiones en común?
    comunes = set(df_nr["accession"]) & set(df_pdb["accession"]) if len(df_nr) and len(df_pdb) else set()
    print(f"\nAccesiones en común nr∩pdb: {len(comunes)}")
    if comunes:
        print(sorted(list(comunes))[:10])
    print(f"Tablas en {out_dir}/")


if __name__ == "__main__":
    main()
