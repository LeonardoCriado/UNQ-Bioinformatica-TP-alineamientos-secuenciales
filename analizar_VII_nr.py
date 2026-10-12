#!/usr/bin/env python3
"""
Análisis del Desafío VII a partir del XML de BLASTp vs NCBI nr.

Parsea el XML por streaming (iterparse, apto para decenas de MB),
calcula por hit:
  - mejor HSP (menor E-value): evalue, bit score, % identidad
  - cobertura del HIT: unión de los intervalos de query de todos sus HSP
    / largo de la query (es la lectura más fiel de "cobertura" de la
    consigna; también se reporta la cobertura del mejor HSP para comparar)
Aplica los filtros de la consigna (E-value <= 100, cobertura >= 70%),
guarda el TSV y genera los 3 gráficos en figures/.

Uso:
    python analizar_VII_nr.py --xml resultados_VII_nr/VII_nr.xml --out-tsv resultados_VII_nr/VII_nr_cob70.tsv --fig-dir figures

Sobre E=0 en los gráficos: el XML trae 0.0 cuando el E-value está por debajo
de la resolución numérica (~1e-180). Esos puntos se grafican en un piso
(mínimo E positivo / 10) y el gráfico lo dice explícitamente con el conteo.
No se inventa ningún valor: el TSV conserva el 0.0 original.
"""

import argparse
import sys
from pathlib import Path
from xml.etree.ElementTree import iterparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

QUERY_LEN = 133


def _union_len(intervals):
    """Largo de la unión de intervalos (inicio, fin) 1-indexed inclusivos."""
    ivs = sorted(intervals)
    total, cur_s, cur_e = 0, None, None
    for s, e in ivs:
        if cur_s is None:
            cur_s, cur_e = s, e
        elif s <= cur_e + 1:
            cur_e = max(cur_e, e)
        else:
            total += cur_e - cur_s + 1
            cur_s, cur_e = s, e
    if cur_s is not None:
        total += cur_e - cur_s + 1
    return total


def parse_nr_xml(xml_path: Path, query_len: int = QUERY_LEN) -> pd.DataFrame:
    rows = []
    n_hits_xml = 0
    for event, hit in iterparse(str(xml_path), events=("end",)):
        if hit.tag != "Hit":
            continue
        n_hits_xml += 1
        acc = hit.findtext("Hit_accession", default="")
        descr = (hit.findtext("Hit_def", default="") or "")[:120]
        hlen = int(hit.findtext("Hit_len", default="0") or 0)
        hsps = []
        for hsp in hit.iter("Hsp"):
            try:
                hsps.append({
                    "e": float(hsp.findtext("Hsp_evalue", default="10") or 10),
                    "bits": float(hsp.findtext("Hsp_bit-score", default="0") or 0),
                    "ident": int(float(hsp.findtext("Hsp_identity", default="0") or 0)),
                    "alen": int(hsp.findtext("Hsp_align-len", default="0") or 0),
                    "qs": int(hsp.findtext("Hsp_query-from", default="0") or 0),
                    "qe": int(hsp.findtext("Hsp_query-to", default="0") or 0),
                })
            except (TypeError, ValueError):
                continue
        hit.clear()
        if not hsps:
            continue
        best = min(hsps, key=lambda h: h["e"])
        cov_hit = 100.0 * _union_len([(h["qs"], h["qe"]) for h in hsps]) / query_len
        cov_best = 100.0 * best["alen"] / query_len
        ident_pct = 100.0 * best["ident"] / best["alen"] if best["alen"] else 0.0
        rows.append({
            "accession": acc, "descripcion": descr, "longitud": hlen,
            "n_hsps": len(hsps), "evalue": best["e"], "bit_score": best["bits"],
            "identidades": best["ident"], "alineamiento_len": best["alen"],
            "identidad_%": round(ident_pct, 2),
            "cobertura_hit_%": round(cov_hit, 2),
            "cobertura_mejor_hsp_%": round(cov_best, 2),
        })
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values("evalue").reset_index(drop=True)
    return df, n_hits_xml


def scatter_log(ax, x, y, floor_note=None):
    ax.scatter(x, y, s=8, alpha=0.5)
    ax.set_yscale("log")
    if floor_note:
        ax.text(0.02, 0.02, floor_note, transform=ax.transAxes, fontsize=7,
                va="bottom", ha="left",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.9))


def main():
    ap = argparse.ArgumentParser(description="VII vs nr: TSV + figures/")
    ap.add_argument("--xml", required=True)
    ap.add_argument("--out-tsv", default="resultados_VII_nr/VII_nr_cob70.tsv")
    ap.add_argument("--fig-dir", default="figures")
    ap.add_argument("--cobertura", type=float, default=70.0)
    ap.add_argument("--expect", type=float, default=100.0)
    ap.add_argument("--query-len", type=int, default=QUERY_LEN)
    args = ap.parse_args()

    df, n_xml = parse_nr_xml(Path(args.xml), args.query_len)
    print(f"Hits en el XML: {n_xml} | filas parseadas: {len(df)}")
    n_e = int((df["evalue"] <= args.expect).sum())
    print(f"Con E-value <= {args.expect:g}: {n_e}")
    df_f = df[(df["evalue"] <= args.expect)
              & (df["cobertura_hit_%"] >= args.cobertura)].reset_index(drop=True)
    print(f"Con cobertura_hit >= {args.cobertura:g}%: {len(df_f)} (final en gráficos)")
    print(f"Hits con >1 HSP: {int((df['n_hsps'] > 1).sum())} "
          f"(en esos, la cobertura por hit puede superar a la del mejor HSP)")
    n_cob_best = int(((df["evalue"] <= args.expect)
                      & (df["cobertura_mejor_hsp_%"] >= args.cobertura)).sum())
    print(f"(Referencia: con cobertura del mejor HSP >= {args.cobertura:g}%: {n_cob_best})")

    out_tsv = Path(args.out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    df_f.to_csv(out_tsv, sep="\t", index=False)
    print(f"Guardado {out_tsv}")

    fig_dir = Path(args.fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)
    n = len(df_f)
    # Piso para E=0 en escala log (transparente y anotado)
    pos = df_f.loc[df_f["evalue"] > 0, "evalue"]
    floor = float(pos.min()) / 10.0 if len(pos) else 1e-300
    n_cero = int((df_f["evalue"] == 0.0).sum())
    y = df_f["evalue"].replace(0.0, floor)
    nota = (f"E=0 bajo resolución numérica: n={n_cero}, "
            f"graficados en el piso {floor:.1e}") if n_cero else None
    base = f"BLASTp vs nr, Expect=100, cob.hit>=70% (n={n})"

    fig, ax = plt.subplots(figsize=(7, 5))
    scatter_log(ax, df_f["identidad_%"], y, nota)
    ax.set_xlabel("% identidad (mejor HSP)")
    ax.set_ylabel("E-value (escala log)")
    ax.set_title(f"VII: E-value vs % identidad\n{base}")
    fig.tight_layout()
    fig.savefig(fig_dir / "desafio-vii-evalue-vs-identidad.png", dpi=170)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(df_f["identidad_%"], df_f["bit_score"], s=8, alpha=0.5)
    ax.set_xlabel("% identidad (mejor HSP)")
    ax.set_ylabel("Bit score")
    ax.set_title(f"VII: Score vs % identidad\n{base}")
    fig.tight_layout()
    fig.savefig(fig_dir / "desafio-vii-score-vs-identidad.png", dpi=170)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 5))
    scatter_log(ax, df_f["bit_score"], y, nota)
    ax.set_xlabel("Bit score (mejor HSP)")
    ax.set_ylabel("E-value (escala log)")
    ax.set_title(f"VII: Score vs E-value\n{base}")
    fig.tight_layout()
    fig.savefig(fig_dir / "desafio-vii-score-vs-evalue.png", dpi=170)
    plt.close(fig)
    print(f"Gráficos en {fig_dir}/desafio-vii-*.png")


if __name__ == "__main__":
    main()
