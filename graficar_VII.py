#!/usr/bin/env python3
"""
Genera los 3 gráficos del Desafío VII en figures/ a partir del TSV guardado.

Lee `resultados_VII/resultados_VII.tsv` (78 hits, BLASTp vs Swiss-Prot,
Expect=100, cobertura >= 70%) y escribe:
    figures/desafio-vii-evalue-vs-identidad.png
    figures/desafio-vii-score-vs-identidad.png
    figures/desafio-vii-score-vs-evalue.png

Sobre E=0: si el TSV trajera E-values 0.0 (por debajo de la resolución
numérica), se grafican en un piso (mínimo E positivo / 10) y el gráfico lo
anota con el conteo. En los datos actuales no hay ningún E=0
(mínimo 3.8e-97), así que el gráfico muestra los valores tal cual.

Uso:
    python graficar_VII.py
    python graficar_VII.py --tsv resultados_VII/resultados_VII.tsv --fig-dir figures
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

BASE_DEFAULT = "BLASTp vs Swiss-Prot (análisis alternativo, n=78)"


def main():
    ap = argparse.ArgumentParser(description="VII: figuras desde el TSV guardado")
    ap.add_argument("--tsv", default="resultados_VII/resultados_VII.tsv")
    ap.add_argument("--fig-dir", default="figures")
    ap.add_argument("--base", default=None,
                    help="Texto de fuente bajo el título (default: Swiss-Prot n=78)")
    args = ap.parse_args()

    df = pd.read_csv(args.tsv, sep="\t")
    base = args.base or f"BLASTp vs Swiss-Prot (análisis alternativo, n={len(df)})"
    fig_dir = Path(args.fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)

    pos = df.loc[df["evalue"] > 0, "evalue"]
    floor = float(pos.min()) / 10.0 if len(pos) else 1e-300
    n_cero = int((df["evalue"] == 0.0).sum())
    y = df["evalue"].replace(0.0, floor)
    nota = (f"E=0 bajo resolución numérica: n={n_cero}, "
            f"graficados en el piso {floor:.1e}") if n_cero else None

    def titulo(ax, linea1):
        ax.set_title(linea1, fontsize=12, pad=18)
        ax.text(0.5, 1.02, base, transform=ax.transAxes, fontsize=9,
                color="dimgray", ha="center", va="bottom")
        if nota:
            ax.text(0.02, 0.02, nota, transform=ax.transAxes, fontsize=7,
                    va="bottom", ha="left",
                    bbox=dict(boxstyle="round,pad=0.3", fc="white",
                              ec="gray", alpha=0.9))

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(df["identidad_%"], y, s=10, alpha=0.6)
    ax.set_yscale("log")
    titulo(ax, "VII: E-value vs % identidad")
    ax.set_xlabel("% identidad (mejor HSP)")
    ax.set_ylabel("E-value (adimensional, escala log)")
    fig.tight_layout()
    fig.savefig(fig_dir / "desafio-vii-evalue-vs-identidad.png", dpi=170)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(df["identidad_%"], df["bit_score"], s=10, alpha=0.6)
    titulo(ax, "VII: Score vs % identidad")
    ax.set_xlabel("% identidad (mejor HSP)")
    ax.set_ylabel("Bit score")
    fig.tight_layout()
    fig.savefig(fig_dir / "desafio-vii-score-vs-identidad.png", dpi=170)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(df["bit_score"], y, s=10, alpha=0.6)
    ax.set_yscale("log")
    titulo(ax, "VII: Score vs E-value")
    ax.set_xlabel("Bit score (mejor HSP)")
    ax.set_ylabel("E-value (adimensional, escala log)")
    fig.tight_layout()
    fig.savefig(fig_dir / "desafio-vii-score-vs-evalue.png", dpi=170)
    plt.close(fig)
    print(f"n={len(df)}, E=0: {n_cero} | Figuras en {fig_dir}/desafio-vii-*.png")


if __name__ == "__main__":
    main()
