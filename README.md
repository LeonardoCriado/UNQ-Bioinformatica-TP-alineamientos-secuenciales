# TP Alineamientos – Introducción a la Bioinformática

Respuestas Desafíos VII–X (Búsqueda de similitud secuencial) + scripts Python.

## Contenido

- `TP_Alineamientos_Respuestas.md` — respuestas teóricas listas para GDocs
- `blast_common.py` — parseo de XML BLAST + filtros + plots (común a VII–X)
- `desafio_VII_blast.py` — BLASTp secuencia priónica completa + 3 gráficos
- `desafio_VIII_fragmentos.py` — comparativa query completa / mitad / cuarto
- `desafio_IX_identificacion.py` — identificación proteína 262 aa vs nr
- `desafio_X_pdb.py` — misma secuencia vs PDB + comparativa nr vs PDB
- `verificar_X_pdb_rcsb.py` — verificación directa de X por RCSB Data API (sin BLAST)
- `blast_local.py` — BLASTp local con NCBI BLAST+ contra Swiss-Prot (reproducible)
- `blast_nr_job.py` — BLASTp VII vs NCBI nr por API URL (PUT/poll/GET reanudable)
- `analizar_VII_nr.py` — parseo por streaming del XML de nr + figuras con E=0 anotado
- `graficar_VII.py` — regenera `figures/` desde el TSV guardado
- `resultados_VII/`, `resultados_VIII/`, `resultados_IX/`, `resultados_X/` — TSVs + PNGs
- `figures/` — los 3 gráficos del desafío VII embebidos en el informe
- `resultados_VII_nr/` — intento de BLAST vs nr (RID + log; sin XML: no completó)
- `requirements.txt`

## Uso rápido

```bash
pip install -r requirements.txt

# VII (offline, con XML ya descargado):
python desafio_VII_blast.py --xml blast_results/VII_full_nr.xml --out-dir resultados_VII

# VIII:
python desafio_VIII_fragmentos.py --xml-full blast_results/VIII_full.xml --xml-half blast_results/VIII_half.xml --xml-quarter blast_results/VIII_quarter.xml --out-dir resultados_VIII

# IX:
python desafio_IX_identificacion.py --xml blast_results/IX_seq2_nr.xml --out-dir resultados_IX

# X (vía BLAST):
python desafio_X_pdb.py --xml-nr blast_results/IX_seq2_nr.xml --xml-pdb blast_results/X_seq2_pdb.xml --out-dir resultados_X

# X (vía directa RCSB, sin BLAST):
python verificar_X_pdb_rcsb.py --out resultados_X/X_pdb_verificacion.tsv

# VII, figuras del informe desde el TSV guardado:
python graficar_VII.py --out-tsv resultados_VII/resultados_VII.tsv --fig-dir figures
```

Para correr BLAST nuevos contra NCBI (requiere internet, tarda varios minutos):

```bash
python desafio_VII_blast.py --run-blast --hits 5000 --expect 100 --out-dir resultados_VII
python desafio_IX_identificacion.py --run-blast --hits 1000 --expect 10 --out-dir resultados_IX
python desafio_X_pdb.py --run-blast --hits 500 --expect 10 --out-dir resultados_X
```

Para BLAST local reproducible (requiere `blast+`, ej. `brew install blast`):

```bash
# Descargar Swiss-Prot y armar la base:
curl -L -o uniprot_sprot.fasta.gz https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
makeblastdb -in uniprot_sprot.fasta.gz -dbtype prot -out uniprot_sprot -parse_seqids
# Correr las 4 queries:
python blast_local.py --db uniprot_sprot --out resultados_local
```

Nota sobre `hits` y base de datos (VII): la consigna pide 20000 hits pero no
nombra base (se asumió `nr`, el default de BLASTp). La búsqueda contra `nr`
no pudo completarse: los RIDs vencieron en cola y el reenvío quedó en
`WAITING` (RTOE ≈ 6 h) hasta que se canceló el polling; el BLAST local contra
una copia de `nr` es inviable en esta máquina (>100 GB de descarga). Los
resultados incluidos son de BLASTp local vs Swiss-Prot (575.748 seqs,
Expect=100, Max targets=5000, cobertura 70%: 138 hits, 78 en gráficos),
identificados como análisis alternativo, no como cumplimiento de la consigna.
Para correr el VII contra `nr` cuando haya tiempo de cola:
`python blast_nr_job.py put/poll/get` y `python analizar_VII_nr.py --xml ...`.
Ver informe para el detalle del bloqueo.

Fuente TP: https://github.com/AJVelezRueda/Introduccion_a_la_Bioinformatica/blob/master/Teorico_Practicos/Bioinform%C3%A1tica_Primeros_pasos/Alineamientos_secuenciales.md
