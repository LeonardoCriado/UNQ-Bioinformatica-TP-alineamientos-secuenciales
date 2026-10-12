#!/usr/bin/env python3
"""
Funciones comunes para los desafíos VII-X (TP Alineamientos).

Parsea XML de BLASTp (Biopython NCBIXML), filtra por cobertura
y genera tablas/plots. Pensado para trabajar offline a partir
de un XML ya descargado, sin re-correr BLAST cada vez.

Formato TSV de salida:
  accession | descripcion | longitud | evalue | bit_score | score | identidad_% | cobertura_% | alineamiento_len
"""

from pathlib import Path
from xml.etree import ElementTree as ET

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def parse_blast_xml(xml_path: str | Path, query_len: int) -> pd.DataFrame:
    """Parsea BLAST XML (formato NCBI) y devuelve DataFrame con un fila por HSP.

    Cobra importancia el mejor HSP de cada hit (el de menor e-value).
    Calcula cobertura como alin_len / query_len * 100.
    """
    xml_path = Path(xml_path)
    if not xml_path.exists():
        raise FileNotFoundError(f"No existe {xml_path}")
    # Intentar con Biopython primero (más robusto con variantes de XML)
    try:
        from Bio.Blast import NCBIXML
        with open(xml_path) as f:
            records = list(NCBIXML.parse(f))
        rows = []
        for rec in records:
            for aln in rec.alignments:
                # mejor HSP = menor e-value
                best = min(aln.hsps, key=lambda h: h.expect)
                ident_pct = 100.0 * best.identities / best.align_length if best.align_length else 0.0
                cov_pct = 100.0 * best.align_length / query_len if query_len else 0.0
                rows.append({
                    "accession": aln.accession,
                    "descripcion": aln.title[:120] if aln.title else aln.hit_id,
                    "longitud": aln.length,
                    "evalue": float(best.expect),
                    "bit_score": float(best.bits),
                    "score": float(best.score),
                    "identidades": int(best.identities),
                    "alineamiento_len": int(best.align_length),
                    "identidad_%": round(ident_pct, 2),
                    "cobertura_%": round(cov_pct, 2),
                })
        df = pd.DataFrame(rows)
        if not df.empty:
            df = df.sort_values("evalue").reset_index(drop=True)
        return df
    except Exception:
        pass
    # Fallback: parseo manual mínimo del XML
    tree = ET.parse(xml_path)
    root = tree.getroot()
    rows = []
    for hit in root.iter("Hit"):
        acc = hit.findtext("Hit_accession", default="")
        descr = hit.findtext("Hit_def", default="")[:120]
        hlen = int(hit.findtext("Hit_len", default="0") or 0)
        best = None
        for hsp in hit.iter("Hsp"):
            e = float(hsp.findtext("Hsp_evalue", default="10") or 10)
            if best is None or e < best[0]:
                best = (e,
                        float(hsp.findtext("Hsp_bit-score", default="0") or 0),
                        float(hsp.findtext("Hsp_score", default="0") or 0),
                        int(float(hsp.findtext("Hsp_identity", default="0") or 0)),
                        int(hsp.findtext("Hsp_align-len", default="0") or 0))
        if best is None:
            continue
        e, bits, score, ident, alen = best
        ident_pct = 100.0 * ident / alen if alen else 0.0
        cov_pct = 100.0 * alen / query_len if query_len else 0.0
        rows.append({
            "accession": acc, "descripcion": descr, "longitud": hlen,
            "evalue": e, "bit_score": bits, "score": score,
            "identidades": ident, "alineamiento_len": alen,
            "identidad_%": round(ident_pct, 2), "cobertura_%": round(cov_pct, 2),
        })
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values("evalue").reset_index(drop=True)
    return df


def filtrar_cobertura(df: pd.DataFrame, minimo: float = 70.0) -> pd.DataFrame:
    """Filtra hits con cobertura >= minimo. Devuelve copia."""
    if df.empty:
        return df.copy()
    return df[df["cobertura_%"] >= minimo].reset_index(drop=True)


def resumen(df: pd.DataFrame) -> dict:
    """Estadísticas rápidas para el informe."""
    if df.empty:
        return {"n": 0}
    return {
        "n": len(df),
        "evalue_min": float(df["evalue"].min()),
        "evalue_max": float(df["evalue"].max()),
        "ident_min": float(df["identidad_%"].min()),
        "ident_max": float(df["identidad_%"].max()),
        "score_min": float(df["bit_score"].min()),
        "score_max": float(df["bit_score"].max()),
        "n_ident100": int((df["identidad_%"] >= 99.5).sum()),
        "n_evalue0": int((df["evalue"] == 0.0).sum()),
    }


def scatter(df: pd.DataFrame, x: str, y: str, out: str | Path,
            xlabel: str | None = None, ylabel: str | None = None,
            title: str | None = None, logy: bool = False, logx: bool = False):
    """Scatter simple y legible para el TP."""
    out = Path(out)
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(df[x], df[y], s=10, alpha=0.6)
    ax.set_xlabel(xlabel or x)
    ax.set_ylabel(ylabel or y)
    if title:
        ax.set_title(title)
    if logy:
        ax.set_yscale("log")
    if logx:
        ax.set_xscale("log")
    fig.tight_layout()
    fig.savefig(out, dpi=170)
    plt.close(fig)
    return str(out)
