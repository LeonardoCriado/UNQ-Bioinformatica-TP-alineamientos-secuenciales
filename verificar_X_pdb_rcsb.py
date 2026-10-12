#!/usr/bin/env python3
"""
Verificación del Desafío X por vía directa (RCSB Data API).

Compara la query de 262 aa (LpxA) contra las secuencias depositadas de
varias entradas PDB asociadas a P0A722, sin pasar por BLAST.
Descarga el `pdbx_seq_one_letter_code_can` de cada polymer entity y
calcula identidad y cobertura exactas.

Uso:
    python verificar_X_pdb_rcsb.py
    python verificar_X_pdb_rcsb.py --pdbs 1LXA 2AQ9 2QIA 4J09 --out resultados_X/X_pdb_verificacion.tsv

Requiere: requests (o urllib estándar si no está instalada).
Fuente: RCSB Data API https://data.rcsb.org/ (consultado 2026-10-12).
"""

import argparse
import json
import sys
from pathlib import Path

SEQ2 = ("MIDKSAFVHPTAIVEEGASIGANAHIGPFCIVGPHVEIGEGTVLKSHVVVNGHTKIGRDNEIYQFASIGEVNQDLKYAGEPTRVEIGDRNRIRESVTIHRGTVQGGGLTKVGSDNLLMINAHIAHDCTVGNRCILANNATLAGHVSVDDFAIIGGMTAVHQFCIIGAHVMVGGCSGVAQDVPPYVIAQGNHATPFGVNIEGLKRRGFSREAITAIRNAYKLIYRSGKTLDEVKPEIAELAETYPEVKAFTDFFARSTRGLIR")


def get_json(url: str, timeout: int = 30) -> dict:
    try:
        import requests
        r = requests.get(url, headers={"User-Agent": "TP-UNQ/1.0"}, timeout=timeout)
        r.raise_for_status()
        return r.json()
    except ImportError:
        import urllib.request
        req = urllib.request.Request(url, headers={"User-Agent": "TP-UNQ/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as fh:
            return json.load(fh)


def main():
    ap = argparse.ArgumentParser(description="Desafío X: verificación directa vs PDB (RCSB API)")
    ap.add_argument("--pdbs", nargs="+", default=["1LXA", "2AQ9", "2QIA", "4J09"],
                    help="IDs PDB a comparar (default: 4 LpxA de E. coli)")
    ap.add_argument("--out", default="resultados_X/X_pdb_verificacion.tsv",
                    help="TSV de salida")
    args = ap.parse_args()

    filas = []
    for pdb in args.pdbs:
        pdb = pdb.strip().upper()
        try:
            entry = get_json(f"https://data.rcsb.org/rest/v1/core/entry/{pdb}")
        except Exception as e:
            print(f"# Aviso: no se pudo obtener {pdb}: {e}", file=sys.stderr)
            continue
        title = (entry.get("struct", {}).get("title", "?") or "?")[:100]
        info = entry.get("rcsb_entry_info", {}) or {}
        res = info.get("diffrn_resolution_high", "?")
        if isinstance(res, dict):
            res = res.get("value", "?")
        exp = ((entry.get("exptl") or [{}])[0] or {}).get("method", "?")
        ids = (entry.get("rcsb_entry_container_identifiers", {}) or {}).get("polymer_entity_ids", [])
        for ent in ids:
            try:
                pe = get_json(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb}/{ent}")
            except Exception as e:
                print(f"# Aviso: {pdb}_{ent}: {e}", file=sys.stderr)
                continue
            seq = (pe.get("entity_poly", {}) or {}).get("pdbx_seq_one_letter_code_can", "") or ""
            if len(seq) != len(SEQ2):
                continue  # péptidos ligandos u otras cadenas, no la LpxA completa
            ident = sum(a == b for a, b in zip(seq, SEQ2))
            orgs = ";".join(o.get("scientific_name", "?")
                            for o in pe.get("rcsb_entity_source_organism", []))
            filas.append({
                "pdb_chain": f"{pdb}_{ent}", "titulo": title,
                "metodo": exp, "resolucion_A": res, "organismo": orgs,
                "largo": len(seq), "identicas": ident,
                "identidad_porc": round(100.0 * ident / len(SEQ2), 2),
                "cobertura_porc": 100.0,
                "rcsb_url": f"https://www.rcsb.org/structure/{pdb}",
            })
            print(f"{pdb}_{ent}: {ident}/{len(SEQ2)} idénticas ({100.0*ident/len(SEQ2):.1f}%) {orgs}")

    if not filas:
        print("Sin filas para guardar.", file=sys.stderr)
        sys.exit(1)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    import csv
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(filas[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(filas)
    print(f"Guardado {out} ({len(filas)} filas)")


if __name__ == "__main__":
    main()
