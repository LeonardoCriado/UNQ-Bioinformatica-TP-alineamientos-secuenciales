#!/usr/bin/env python3
"""
Desafío VII - BLASTp de la secuencia problema (133 aa, fragmento de prion).

Consigna: BLAST con 20000 hits, e-value 100, filtrando hits con
mínimo 70% de cobertura. Analizar E-value vs % identidad,
Score vs % identidad y Score vs E-value.

Uso offline (recomendado, reproduce plots sin re-correr BLAST):
    python desafio_VII_blast.py --xml blast_results/VII_full_nr.xml --out-dir resultados_VII

Para lanzar BLAST nuevo (requiere internet, tarda varios minutos):
    python desafio_VII_blast.py --run-blast --hits 5000 --expect 100

La secuencia es el extremo C-terminal de la proteína priónica
(residuos 121-253 de PRNP humano P04156, verificado localmente).
Por eso se esperan miles de hits de priones de mamíferos con
alta identidad y E-value ~0.

Fecha de consulta NCBI: 2026-10-12 (las bases se actualizan).
"""

import argparse
import sys
from pathlib import Path

SEQ_FULL = "VVGGLGGYMLGSAMSRPIIHFGSDYEDRYYRENMHRYPNQVYYRPMDEYSNQNNFVHDCVNITIKQHTVTTTTKGENFTETDVKMMERVVEQMCITQYERESQAYYQRGSSMVLFSSPPVILLISFLIFLIVG"

try:
    from blast_common import parse_blast_xml, filtrar_cobertura, resumen, scatter
except ImportError:
    sys.path.insert(0, str(Path(__file__).parent))
    from blast_common import parse_blast_xml, filtrar_cobertura, resumen, scatter

import pandas as pd


def run_blast(seq: str, expect: int, hits: int, out_xml: Path):
    from Bio.Blast import NCBIWWW
    print(f"Enviando BLASTp a NCBI (nr, expect={expect}, hits={hits}, len={len(seq)})...")
    handle = NCBIWWW.qblast("blastp", "nr", seq, expect=expect, hitlist_size=hits, format_type="XML")
    out_xml.write_text(handle.read())
    print(f"Guardado {out_xml}")


def main():
    ap = argparse.ArgumentParser(description="Desafío VII: BLASTp + plots")
    ap.add_argument("--xml", help="XML de BLAST ya descargado (modo offline)")
    ap.add_argument("--out-dir", default="resultados_VII", help="Carpeta de salida")
    ap.add_argument("--cobertura", type=float, default=70.0, help="Cobertura minima en porciento (default 70)")
    ap.add_argument("--run-blast", action="store_true", help="Correr BLAST nuevo contra NCBI")
    ap.add_argument("--hits", type=int, default=5000, help="Max target seqs (consigna pide 20000; se usa 5000 por tiempo/tamaño, ver informe)")
    ap.add_argument("--expect", type=float, default=100, help="E-value umbral (default 100)")
    ap.add_argument("--tsv", default="resultados_VII.tsv")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.run_blast:
        xml_path = out_dir / "VII_full_nr.xml"
        run_blast(SEQ_FULL, args.expect, args.hits, xml_path)
    else:
        if not args.xml:
            ap.print_help()
            print("\nEjemplo offline: python desafio_VII_blast.py --xml blast_results/VII_full_nr.xml")
            sys.exit(1)
        xml_path = Path(args.xml)

    df = parse_blast_xml(xml_path, len(SEQ_FULL))
    print(f"Hits totales parseados: {len(df)}")
    if df.empty:
        print("XML sin hits. Revisar el archivo.", file=sys.stderr)
        sys.exit(1)

    df_f = filtrar_cobertura(df, args.cobertura)
    print(f"Hits con cobertura >= {args.cobertura}%: {len(df_f)}")
    r = resumen(df_f)
    print(f"Resumen filtrado: {r}")

    tsv = out_dir / args.tsv
    df_f.to_csv(tsv, sep="\t", index=False)
    print(f"Guardado {tsv}")

    # Plots pedidos en la consigna
    scatter(df_f, "identidad_%", "evalue", out_dir / "VII_Evalue_vs_identidad.png",
            xlabel="% identidad", ylabel="E-value",
            title="VII: E-value vs % identidad (eje Y log)", logy=True)
    scatter(df_f, "identidad_%", "bit_score", out_dir / "VII_Score_vs_identidad.png",
            xlabel="% identidad", ylabel="Bit score",
            title="VII: Score vs % identidad")
    # Score vs E-value: E en log para que se vea la nube (muchos E=0 se llevan a 1e-200 solo para graficar)
    import numpy as np
    df_plot = df_f.copy()
    df_plot["evalue_plot"] = df_plot["evalue"].replace(0.0, 1e-200)
    scatter(df_plot, "bit_score", "evalue_plot", out_dir / "VII_Score_vs_Evalue.png",
            xlabel="Bit score", ylabel="E-value (0 graficado como 1e-200)",
            title="VII: Score vs E-value (eje Y log)", logy=True)
    print(f"Gráficos en {out_dir}/VII_*.png")


if __name__ == "__main__":
    main()
