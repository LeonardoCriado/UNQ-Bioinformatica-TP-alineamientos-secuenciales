#!/usr/bin/env python3
"""
BLASTp contra NCBI nr para el Desafío VII (20000 hits, Expect=100).

Usa la API URL de NCBI BLAST (PUT -> poll -> GET) con urllib de la
librería estándar, sin Biopython, para poder reanudar el polling si
la sesión se interrumpe: el RID queda guardado en un archivo.

Uso:
    python blast_nr_job.py put      # envía y guarda VII_nr.rid
    python blast_nr_job.py poll     # espera hasta READY (reanudable)
    python blast_nr_job.py get      # descarga el XML (streaming)

Archivos (en --work-dir, default resultados_VII_nr/):
    VII_nr.rid  -> RID + metadatos del envío
    VII_nr.xml  -> resultado completo (puede ser decenas de MB)

Parámetros fijos de la consigna VII:
    PROGRAM=blastp, DATABASE=nr, EXPECT=100, HITLIST_SIZE=20000.
La consigna no nombra base de datos; se usa nr por ser el default de BLASTp.
"""

import argparse
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

SEQ_FULL = ("VVGGLGGYMLGSAMSRPIIHFGSDYEDRYYRENMHRYPNQVYYRPMDEYSNQNNFVHDCVNITIKQHTVTTTTKGEN"
            "FTETDVKMMERVVEQMCITQYERESQAYYQRGSSMVLFSSPPVILLISFLIFLIVG")
BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"


def _post(params: dict) -> str:
    data = urllib.parse.urlencode(params).encode()
    req = urllib.request.Request(BLAST_URL, data=data,
                                 headers={"User-Agent": "TP-UNQ-Bioinfo/1.0"})
    with urllib.request.urlopen(req, timeout=60) as fh:
        return fh.read().decode("utf-8", "replace")


def _parse_qblast_info(text: str) -> dict:
    # Tolerante a indentación y a las variantes con/sin "<!--" (el PUT y el
    # GET devuelven el bloque QBlastInfo con formatos distintos).
    import re
    m = re.search(r"QBlastInfoBegin(.*?)QBlastInfoEnd", text, re.S)
    if not m:
        return {}
    info = {}
    for line in m.group(1).splitlines():
        s = line.strip()
        if "=" in s:
            k, v = s.split("=", 1)
            info[k.strip()] = v.strip()
    return info


def cmd_put(args) -> int:
    print("Enviando PUT a NCBI (blastp, nr, EXPECT=100, HITLIST_SIZE=20000)...", flush=True)
    resp = _post({"CMD": "Put", "PROGRAM": "blastp", "DATABASE": "nr",
                  "QUERY": SEQ_FULL, "EXPECT": "100", "HITLIST_SIZE": "20000",
                  "FILTER": "m", "FORMAT_TYPE": "XML"})
    info = _parse_qblast_info(resp)
    rid = info.get("RID", "")
    rtoe = info.get("RTOE", "?")
    if not rid:
        print("ERROR: NCBI no devolvió RID. Respuesta:", file=sys.stderr)
        print(resp[:2000], file=sys.stderr)
        return 1
    rid_file = args.work_dir / "VII_nr.rid"
    rid_file.write_text(f"RID={rid}\nRTOE={rtoe}\nSUBMITTED=2026-10-12\n"
                        f"PROGRAM=blastp DATABASE=nr EXPECT=100 HITLIST_SIZE=20000\n")
    print(f"RID={rid} RTOE={rtoe}s -> {rid_file}", flush=True)
    return 0


def _status(rid: str) -> str:
    resp = _post({"CMD": "Get", "RID": rid, "FORMAT_TYPE": "XML"})
    info = _parse_qblast_info(resp)
    return info.get("Status", "UNKNOWN")


def cmd_poll(args) -> int:
    rid = args.work_dir.joinpath("VII_nr.rid").read_text().split("RID=")[1].split()[0]
    print(f"Polling RID={rid} ...", flush=True)
    waited = 0
    step = 60
    while True:
        try:
            st = _status(rid)
        except Exception as e:
            print(f"[{waited}s] error de red: {e} (reintento)", flush=True)
            time.sleep(step)
            waited += step
            continue
        print(f"[{waited}s] Status={st}", flush=True)
        if st == "READY":
            print("READY.", flush=True)
            return 0
        if st in ("FAILED", "UNKNOWN"):
            print(f"Estado terminal: {st}", file=sys.stderr)
            return 1
        time.sleep(step)
        waited += step


def cmd_get(args) -> int:
    rid = args.work_dir.joinpath("VII_nr.rid").read_text().split("RID=")[1].split()[0]
    out = args.work_dir / "VII_nr.xml"
    print(f"Descargando resultado RID={rid} -> {out} ...", flush=True)
    params = urllib.parse.urlencode({"CMD": "Get", "RID": rid, "FORMAT_TYPE": "XML"}).encode()
    req = urllib.request.Request(BLAST_URL, data=params,
                                 headers={"User-Agent": "TP-UNQ-Bioinfo/1.0"})
    total = 0
    with urllib.request.urlopen(req, timeout=120) as fh, open(out, "wb") as f:
        while True:
            chunk = fh.read(1024 * 256)
            if not chunk:
                break
            f.write(chunk)
            total += len(chunk)
    print(f"Descargado {total} bytes en {out}", flush=True)
    return 0


def main():
    ap = argparse.ArgumentParser(description="BLASTp VII vs NCBI nr (PUT/poll/GET)")
    ap.add_argument("cmd", choices=["put", "poll", "get"])
    ap.add_argument("--work-dir", default="resultados_VII_nr")
    args = ap.parse_args()
    args.work_dir = Path(args.work_dir)
    args.work_dir.mkdir(parents=True, exist_ok=True)
    fn = {"put": cmd_put, "poll": cmd_poll, "get": cmd_get}[args.cmd]
    sys.exit(fn(args))


if __name__ == "__main__":
    main()
