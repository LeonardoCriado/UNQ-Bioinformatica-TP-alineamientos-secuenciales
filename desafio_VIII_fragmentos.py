#!/usr/bin/env python3
"""
Desafío VIII - Efecto de acortar la query (completa / mitad / cuarto).

Consigna: nuevas búsquedas con la mitad y un cuarto de la secuencia
del desafío VII, comparar gráficos y sacar conclusiones.

Uso:
    python desafio_VIII_fragmentos.py --xml-full X.xml --xml-half Y.xml --xml-quarter Z.xml --out-dir resultados_VIII

Para lanzar BLAST nuevos (lento):
    python desafio_VIII_fragmentos.py --run-blast --hits 2000 --expect 100 --out-dir resultados_VIII

Mitad = primeros 66 aa, cuarto = primeros 33 aa (corte N-terminal del
fragmento, ver informe). Al acortar cae el score máximo posible,
sube el E-value y se pierde especificidad.
"""

import argparse
import sys
from pathlib import Path

SEQ_FULL = "VVGGLGGYMLGSAMSRPIIHFGSDYEDRYYRENMHRYPNQVYYRPMDEYSNQNNFVHDCVNITIKQHTVTTTTKGENFTETDVKMMERVVEQMCITQYERESQAYYQRGSSMVLFSSPPVILLISFLIFLIVG"
SEQ_HALF = SEQ_FULL[:len(SEQ_FULL)//2]
SEQ_QUARTER = SEQ_FULL[:len(SEQ_FULL)//4]

try:
    from blast_common import parse_blast_xml, filtrar_cobertura, scatter
except ImportError:
    sys.path.insert(0, str(Path(__file__).parent))
    from blast_common import parse_blast_xml, filtrar_cobertura, scatter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def run_one(seq, expect, hits, db, out):
    from Bio.Blast import NCBIWWW
    print(f"BLASTp {db} len={len(seq)} ...", flush=True)
    h = NCBIWWW.qblast("blastp", db, seq, expect=expect, hitlist_size=hits, format_type="XML")
    Path(out).write_text(h.read())
    print(f"Guardado {out}", flush=True)


def main():
    ap = argparse.ArgumentParser(description="Desafío VIII: full vs mitad vs cuarto")
    ap.add_argument("--xml-full", help="XML query completa")
    ap.add_argument("--xml-half", help="XML mitad")
    ap.add_argument("--xml-quarter", help="XML cuarto")
    ap.add_argument("--out-dir", default="resultados_VIII")
    ap.add_argument("--cobertura", type=float, default=70.0)
    ap.add_argument("--run-blast", action="store_true")
    ap.add_argument("--hits", type=int, default=2000)
    ap.add_argument("--expect", type=float, default=100)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.run_blast:
        run_one(SEQ_FULL, args.expect, args.hits, "nr", out_dir / "VIII_full.xml")
        run_one(SEQ_HALF, args.expect, args.hits, "nr", out_dir / "VIII_half.xml")
        run_one(SEQ_QUARTER, args.expect, args.hits, "nr", out_dir / "VIII_quarter.xml")
        return

    if not (args.xml_full and args.xml_half and args.xml_quarter):
        ap.print_help()
        print("\nEjemplo: python desafio_VIII_fragmentos.py --xml-full a.xml --xml-half b.xml --xml-quarter c.xml")
        sys.exit(1)

    datas = {}
    for nombre, xml, seq in [("completa", args.xml_full, SEQ_FULL),
                             ("mitad", args.xml_half, SEQ_HALF),
                             ("cuarto", args.xml_quarter, SEQ_QUARTER)]:
        df = parse_blast_xml(xml, len(seq))
        df_f = filtrar_cobertura(df, args.cobertura)
        datas[nombre] = (df, df_f, len(seq))
        print(f"{nombre} (len={len(seq)}): {len(df)} hits totales, {len(df_f)} con cov>={args.cobertura}%")
        df_f.to_csv(out_dir / f"VIII_{nombre}.tsv", sep="\t", index=False)

    # Comparativa: E-value vs %id superpuestos + boxplot de scores
    fig, ax = plt.subplots(figsize=(7, 5))
    for nombre, (df, df_f, _) in datas.items():
        ax.scatter(df_f["identidad_%"], df_f["evalue"].replace(0.0, 1e-200),
                   s=10, alpha=0.5, label=nombre)
    ax.set_yscale("log")
    ax.set_xlabel("% identidad")
    ax.set_ylabel("E-value (0 como 1e-200)")
    ax.set_title("VIII: E-value vs % identidad según largo de query")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "VIII_comparativa_Evalue_vs_identidad.png", dpi=170)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 5))
    for nombre, (df, df_f, _) in datas.items():
        ax.scatter(df_f["identidad_%"], df_f["bit_score"], s=10, alpha=0.5, label=nombre)
    ax.set_xlabel("% identidad")
    ax.set_ylabel("Bit score")
    ax.set_title("VIII: Score vs % identidad según largo de query")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "VIII_comparativa_Score_vs_identidad.png", dpi=170)
    plt.close(fig)

    # Tabla resumen en consola + CSV
    filas = []
    for nombre, (df, df_f, L) in datas.items():
        filas.append({
            "query": nombre, "largo": L, "n_total": len(df), "n_filtrado": len(df_f),
            "evalue_min": float(df_f["evalue"].min()) if len(df_f) else None,
            "evalue_mediana": float(df_f["evalue"].median()) if len(df_f) else None,
            "ident_mediana": float(df_f["identidad_%"].median()) if len(df_f) else None,
            "score_max": float(df_f["bit_score"].max()) if len(df_f) else None,
        })
    pd.DataFrame(filas).to_csv(out_dir / "VIII_resumen.csv", index=False)
    print(pd.DataFrame(filas).to_string(index=False))
    print(f"Salidas en {out_dir}/")


if __name__ == "__main__":
    main()
