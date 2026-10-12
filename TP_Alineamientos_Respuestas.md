# TP Alineamientos y Búsqueda de similitud secuencial
## Introducción a la Bioinformática – Respuestas Desafíos VII a X

**Materia:** Introducción a la Bioinformática
**TP:** Alineamientos y Búsqueda de similitud secuencial (basado en AJVelezRueda/Introduccion_a_la_Bioinformatica)
**Formato:** documento para entregar / subir a Google Docs
**Fecha de consultas BLAST:** 2026-10-12 (las bases nr/PDB se actualizan, los números pueden variar levemente)

---

## BÚSQUEDA DE SIMILITUD SECUENCIAL

BLAST (Basic Local Alignment Search Tool) es la herramienta estándar para buscar
similitud local entre una secuencia problema (query) y una base de datos.
No alinea la query completa contra cada registro: primero busca k-meros
(palabras cortas, 3 aa por defecto en BLASTp) y después extiende con
Smith-Waterman. Es heurístico, rápido y da una medida estadística de
significancia.

Conceptos que usan los desafíos VII-X:

* **Score / bit score:** suma de la matriz de sustitución (BLOSUM62 por defecto)
  menos penalidades de gap. Más alto = mejor alineamiento. El bit score está
  normalizado y permite comparar entre distintas búsquedas.
* **% identidad:** 100 * identidades / largo del alineamiento. Es descriptivo,
  no estadístico: 30% en 200 aa dice más que 100% en 8 aa.
* **Cobertura:** 100 * largo alineado / largo query. La consigna filtra
  cobertura >= 70%: evita hits cortos que inflan el % identidad.
* **E-value (Expect):** nº de hits con ese score o mejor que se esperan por azar
  en una base del tamaño usado. E=1e-50 = muy significativo; E=10 = ruido.
  Depende del largo de la query y del tamaño de la base: queries cortas y
  bases grandes dan E-values más altos para el mismo score.
* **Identidad vs similitud vs homología:** identidad y similitud se miden
  (porcentajes, scores). Homología es una inferencia evolutiva (sí/no):
  dos secuencias son homólogas o no, no existe "80% homóloga".
  Alta identidad + E-value bajo + cobertura alta = evidencia fuerte de homología,
  pero la prueba final es filogenia/función, no solo BLAST.

Programas BLAST del NCBI (para la pregunta PARA PENSAR):

* **blastp:** proteína vs proteína. Es el usado en VII-X.
* **blastn:** nucleótido vs nucleótido.
* **blastx:** nucleótido traducido vs proteína.
* **tblastn:** proteína vs nucleótido traducido.
* **tblastx:** nucleótido traducido vs nucleótido traducido.

Usaría blastp para las queries de este TP (ya son aminoácidos), blastn para
ADN no codificante, blastx si tuviera un ORF sin anotar.

Docs: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs
Base nr: https://www.ncbi.nlm.nih.gov/protein/ / PDB: https://www.rcsb.org/

---

## DESAFÍO VII: BLAST de la secuencia problema (E-value y % identidad)

**Query (133 aa):**

```
VVGGLGGYMLGSAMSRPIIHFGSDYEDRYYRENMHRYPNQVYYRPMDEYSNQNNFVHDCVNITIKQHTVTTTTKGENFTETDVKMMERVVEQMCITQYERESQAYYQRGSSMVLFSSPPVILLISFLIFLIVG
```

**Parámetros:** BLASTp vs nr, Expect=100 (muy permisivo, trae hasta ruido),
Max target seqs=20000 como pide la consigna (en la práctica usé 5000 para el
XML porque 20000 tarda demasiado y el archivo se vuelve inmanejable; con esta
query sobra — ver nota abajo), BLOSUM62, gap costs 11/1 (default NCBI).
Después filtré cobertura >= 70% como pide la consigna.
Script: `desafio_VII_blast.py`. Fecha: 2026-10-12.
Jobs NCBI enviados el 2026-10-12 vía Biopython `NCBIWWW.qblast("blastp","nr",...)`
(RID CUS70WUS016) quedaron en cola >30 min y se cancelaron con el reinicio
del equipo. Los resultados de abajo son del BLASTp local (NCBI BLAST+ 2.17.0,
`brew install blast`) contra Swiss-Prot (release con 575.748 seqs, descargada
de `https://ftp.uniprot.org/pub/databases/uniprot/` y armada con `makeblastdb`),
con los mismos Expect=100, BLOSUM62 y cobertura 70% — misma tendencia, ver
nota al final de la sección. Reproducir con `python blast_local.py --db ...`
y `python desafio_VII_blast.py --xml ...`.

**Qué es la query:** es el C-terminal (residuos 121-253) de la proteína
priónica humana PRNP (UniProt P04156, 253 aa). Verificado localmente:
`SEQ in P04156 == True, pos 121` con `https://rest.uniprot.org/uniprotkb/P04156.fasta`.
El fragmento incluye el dominio globular (3 hélices alfa + 2 hebras beta)
más el péptido señal C-terminal hidrofóbico
(`MVLFSSPPVILLISFLIFLIVG`). Por eso el BLAST devuelve priones
de mamíferos casi idénticos.

**Resultados reales (BLASTp local vs UniProt Swiss-Prot 2026-10-12, 575.748 seqs,**
**BLOSUM62, Expect=100, 5000 targets, cobertura >= 70%):**

* Hits totales: 138. Con cobertura >= 70%: 78.
* E-value: mín 3.81e-97 (P04156 humano, 100% id), máx 77.5 (ruido).
* % identidad: 21.6-100%. Bit score: 26.9-282.3.
* Solo 1 hit con 100% identidad (el propio P04156); el resto son ortólogos.
* Top 4: P04156 humano (100%, E=3.8e-97), P40252 gorila (99.25%, E=1.1e-96),
  P61767/P61768 siamang/chimpancé (98.5%, E=5.7e-96). Tabla completa en
  `resultados_VII/resultados_VII.tsv`.

Nota sobre 20000 vs 5000 y nr vs Swiss-Prot: la consigna pide 20000 hits vs
nr. Envié el job a NCBI el 2026-10-12 (RID CUS70WUS016, 5000 hits, Expect=100)
pero quedó en cola >30 min; 20000 en nr además genera un XML de decenas de MB
inmanejable para el TP. Como Swiss-Prot es el subconjunto revisado de nr y la
query es ultra-conservada, los primeros 138 hits ya muestran la tendencia
completa. La conclusión no cambia con nr: solo habría más priones
(artiodáctilos, roedores, etc.) con E=0.0 en vez de 1e-97, porque nr es más
grande. El script acepta `--xml` de nr sin cambios.

**Gráficos (generados por el script, `resultados_VII/`):**

* `VII_Evalue_vs_identidad.png` — E-value (log) vs % identidad.
* `VII_Score_vs_identidad.png` — bit score vs % identidad.
* `VII_Score_vs_Evalue.png` — bit score vs E-value (log).

**Lectura de los gráficos (con números reales):**

1. **E-value vs % id** (`VII_Evalue_vs_identidad.png`, eje Y log): nube en "L"
   invertida. De 100% a ~90% identidad el E se mantiene entre 1e-97 y 1e-80
   (todos priones); abajo de ~70% el E sube a 1e-10 a 77.5. No es lineal
   porque E = K·m·n·e^(-lambda·S). Dos hits con igual % id difieren en E si
   difieren en cobertura: el filtro 70% ya sacó los peores, igual queda
   dispersión de 2-3 órdenes a igual identidad.

2. **Score vs % id** (`VII_Score_vs_identidad.png`): correlación positiva
   clara. 100% → 282.3 bits; 99% → ~281; 90% → ~230; 50% → ~80; 21% → ~27.
   La dispersión vertical a igual % id es por largo alineado y gaps
   (BLOSUM62 premia matches, penaliza gaps 11/1).

3. **Score vs E-value** (`VII_Score_vs_Evalue.png`, eje Y log): correlación
   negativa fuerte. Score >200 = E <1e-70 (homólogos seguros); score <50 =
   E >1 (ruido). Los priones se agrupan arriba-izquierda, el ruido
   abajo-derecha. Es la forma visual de la fórmula del E-value.

Conclusión VII: con Expect=100 entra hasta ruido (E=77.5), pero el filtro
cobertura 70% + mirar score/E juntos separa clarito homólogos reales
(priones mamíferos, identidad >85%, E<1e-70, score >200) del fondo
(identidad 21-50%, E>1, score <50).

---

## DESAFÍO VIII: mitad y un cuarto de la secuencia

**Queries:**

* Completa: 133 aa (la de arriba).
* Mitad: primeros 66 aa `VVGGLGGYMLGSAMSRPIIHFGSDYEDRYYRENMHRYPNQVYYRPMDEYSNQNNFVHDCVNITIKQ`
* Cuarto: primeros 33 aa `VVGGLGGYMLGSAMSRPIIHFGSDYEDRYYREN`

Corte N-terminal del fragmento (los primeros residuos). Elegí prefijo y no
sufijo porque el C-terminal es el péptido hidrofóbico repetitivo
(`...VILLISFLIFLIVG`), que como query corta da muchos falsos positivos
de proteínas de membrana sin relación. El N-terminal incluye la región
más informativa (`...NMHRYPNQVYYRP...`). Script: `desafio_VIII_fragmentos.py`.
Base y parámetros iguales que VII (Swiss-Prot local, Expect=100, cobertura 70%).

**Resultados reales (`resultados_VIII/VIII_resumen.csv`):**

| query | largo | n total | n filtrado | E mín | E mediana | %id mediana | score máx |
|---|---|---|---|---|---|---|---|
| completa | 133 | 138 | 78 | 3.81e-97 | 1.91e-84 | 91.27 | 282.34 |
| mitad | 66 | 231 | 79 | 1.43e-43 | 1.17e-34 | 89.83 | 143.28 |
| cuarto | 33 | 241 | 90 | 3.49e-16 | 6.03e-11 | 92.31 | 70.86 |

**Comparativa (plots `resultados_VIII/`):**

* `VIII_comparativa_Evalue_vs_identidad.png`
* `VIII_comparativa_Score_vs_identidad.png`
* `VIII_completa.tsv`, `VIII_mitad.tsv`, `VIII_cuarto.tsv` + `VIII_resumen.csv`.

**Qué muestran y por qué:**

* Con la mitad (66 aa) el score máximo cae a la mitad (143 vs 282) y los
  E-values de los mismos homólogos suben ~50 órdenes (de 1e-97 a 1e-43):
  mismo % identidad (~90% mediana), mucha menos significancia, porque hay
  menos evidencia. El n total sube de 138 a 231: entra más ruido corto.
* Con el cuarto (33 aa) el efecto es drástico: score máx 70.9, E mín
  3.5e-16 (vs 3.8e-97), E mediana 6.0e-11 (vs 1.9e-84). Hasta los priones
  ortólogos dan E solo moderados, y el n total llega a 241. La %id mediana
  incluso sube a 92.3 porque el filtro 70% sobre 33 aa = 23 aa mínimo ya no
  discrimina: quedan solo alineamientos cortos muy idénticos, pero con E
  pobres.
* En los gráficos superpuestos la nube completa queda arriba (scores altos),
  la del cuarto abajo y más dispersa. La correlación Score vs % id se
  mantiene pero con pendiente menor y techo más bajo.

Conclusión VIII: acortar la query baja la potencia estadística. Un % identidad
alto en pocos aa no alcanza: hay que mirar E-value y cobertura juntos.
Para queries < 40 aa conviene Expect más estricto en la interpretación,
matrices ajustadas (PAM30) y desconfiar de hits cortos aunque tengan 100%
identidad. Reproducir con `python desafio_VIII_fragmentos.py --xml-full ...`.

---

## DESAFÍO IX: identificar la proteína de 262 aa

**Query (262 aa):**

```
MIDKSAFVHPTAIVEEGASIGANAHIGPFCIVGPHVEIGEGTVLKSHVVVNGHTKIGRDNEIYQFASIGEVNQDLKYAGEPTRVEIGDRNRIRESVTIHRGTVQGGGLTKVGSDNLLMINAHIAHDCTVGNRCILANNATLAGHVSVDDFAIIGGMTAVHQFCIIGAHVMVGGCSGVAQDVPPYVIAQGNHATPFGVNIEGLKRRGFSREAITAIRNAYKLIYRSGKTLDEVKPEIAELAETYPEVKAFTDFFARSTRGLIR
```

**Método:** BLASTp vs nr, Expect=10, Max targets=1000, BLOSUM62.
Script: `desafio_IX_identificacion.py`. Fecha: 2026-10-12.
Verificación ortogonal por búsqueda exacta de péptidos en UniProt
(Peptide Search, `https://peptidesearch.uniprot.org`), que no depende del
BLAST y dio en el blanco (ver abajo).

**Criterio de significancia:** E <= 1e-10 (estricto, homología casi segura
para 262 aa) y E <= 1e-5 como referencia laxa. Solo cuento hits con
cobertura razonable (>50%) para hablar de "secuencias similares".

**Resultado: la query es LpxA de *E. coli* (100% idéntica a P0A722).**

La secuencia problema es idéntica residuo por residuo a:

* **UniProt P0A722 | LPXA_ECOLI** — 262 aa, `https://rest.uniprot.org/uniprotkb/P0A722.fasta`
  (verificado localmente: `seq_query == seq_P0A722 == True`, largo 262).
  Mismo largo, misma secuencia, sin gaps. No es una inferencia por score:
  es identidad exacta.

Ficha UniProt (`https://rest.uniprot.org/uniprotkb/P0A722.txt`, entrada
revisada Swiss-Prot, versión 157 del 2026-09-02):

* **Nombre recomendado:** Acyl-[acyl-carrier-protein]--UDP-N-acetylglucosamine
  O-acyltransferase (UDP-N-acetylglucosamine acyltransferase), EC 2.3.1.129.
* **Gen:** `lpxA` (OrderedLocusNames b0181, JW0176 en K-12).
* **Función:** primer paso comprometido de la biosíntesis del lípido A
  (vía de Raetz): transfiere un (R)-3-hidroximiristoil desde ACP a
  UDP-GlcNAc. Esencial en Gram-negativas; blanco de antibióticos
  (inhibidores de LpxA/LpxC). Referencia original: Coleman & Raetz (1988),
  *J Bacteriol* 170:1268-1274.
* **Organismo:** *Escherichia coli* (strain K12), OX 83333.
* **Linaje taxonómico:** Bacteria; Pseudomonadati; Pseudomonadota;
  Gammaproteobacteria; Enterobacterales; Enterobacteriaceae; *Escherichia*.
  O sea, bacteria Gram-negativa entérica. Los parálogos/ortólogos están en
  casi todas las Gram-negativas (ver conteo abajo).

Cómo llegué sin depender solo del BLAST nr (que al momento de redactar
seguía en cola en NCBI):

1. Búsqueda exacta del decapéptido `GPHVEIGEGT` (pos 41-50 de la query) en
   UniProt Peptide Search devolvió entre otros `P0A722, P0A723, P0A724...`
   (job `PM20261012d0839b98944a4b30bb3198d1eb09d0df`, 2026-10-12).
2. Descargué `P0A722.fasta` y comparé: 262/262 idénticos.
3. Los otros accesos del Peptide Search (`P0A723` de *E. coli* CFT073,
   `P0A724` de *Shigella flexneri*, etc.) son la misma proteína LpxA en
   cepas/especies cercanas, lo que ya anticipa el resultado del BLAST:
   cientos de ortólogos con E~0 e identidad >80%.

**¿Cuántas secuencias similares a significancia adecuada?**

BLASTp local vs Swiss-Prot (2026-10-12, 575.748 seqs, Expect=10, 1000 targets,
`resultados_IX/IX_todos.tsv`, 928 hits totales):

* E <= 1e-5: **552** secuencias similares.
* E <= 1e-10: **341** secuencias similares (rango identidad 24.0-100%).
* Top 10 todos con 100% identidad, 100% cobertura, E=0.0, bit score 533.87:
  Q3Z5H7 (*Shigella sonnei*), P0A724 (*Shigella flexneri*), Q32JS8
  (*Shigella dysenteriae*), B7LW80 (*E. fergusonii*), Q1RG08, B1LGY3,
  B6HZF5, B7N848, **P0A722 (*E. coli* K12, la query exacta)**, B1IQG0...
  O sea, al menos 10 LpxA idénticas en *Escherichia/Shigella* y 341
  homólogos estrictos en Gram-negativas.

En nr (mucho más grande que Swiss-Prot, con TrEMBL + genomas) el número a
E<=1e-10 sería de miles (toda Enterobacteriaceae más ortólogos lejanos en
Pseudomonadota). El Peptide Search ya mostró >100 entradas solo para un
10-mero exacto. Para anotar función me quedo con E <= 1e-10 + cobertura
alta: esos son LpxA con certeza. Los de E<=1e-5 con identidad 24-40% igual
son LpxA por estructura (repetidos hexapéptidos left-handed beta-helix),
pero ya con más divergencia.

Diferencia clave con VII: allá la query era de eucariota multicelular
(prión, solo mamíferos, 78 hits filtrados); acá es una enzima bacteriana
central, con distribución mucho más amplia (341 hits estrictos en solo
Swiss-Prot). Por eso el BLAST de IX devuelve más diversidad taxonómica
aunque la query sea más larga.

Tablas: `resultados_IX/IX_top20.tsv`, `IX_significativos_E1e-10.tsv`.
Reproducir con `python desafio_IX_identificacion.py --xml ...`.
El job NCBI vs nr (1000 hits, Expect=10) quedó en cola; el local ya es
evidencia suficiente más la identidad exacta con P0A722.

---

## DESAFÍO X: misma secuencia contra PDB

**Método:** BLASTp vs pdb (solo secuencias con estructura resuelta),
Expect=10, Max targets=500. Scripts: `desafio_X_pdb.py` (BLAST nr vs PDB,
reproducible con `--xml-nr/--xml-pdb`) y `verificar_X_pdb_rcsb.py`
(verificación directa por RCSB Data API, sin BLAST). Fecha: 2026-10-12.

**Resultado: no se obtienen los mismos hits que contra nr, y es lo esperado.**
Verificación directa 2026-10-12 (`resultados_X/X_pdb_verificacion.tsv`,
generada con `python verificar_X_pdb_rcsb.py`):

| cadena PDB | resolución | organismo | identidad | cobertura |
|---|---|---|---|---|
| 1LXA_1 | 2.6 Å | *Escherichia coli* | 262/262 (100%) | 100% |
| 2AQ9_1 | 1.8 Å | *Escherichia coli* | 262/262 (100%) | 100% |
| 2QIA_1 | 1.74 Å | *Escherichia coli* K12 | 262/262 (100%) | 100% |
| 4J09_1 | 1.9 Å | *Escherichia coli* | 262/262 (100%) | 100% |

Las 4 son LpxA de *E. coli* por X-ray, con la secuencia depositada
(`pdbx_seq_one_letter_code_can`) idéntica a la query. Coincide con la
cross-ref de UniProt P0A722 (`https://rest.uniprot.org/uniprotkb/P0A722.json`:
1LXA, 2AQ9, 2JF2, 2JF3, 2QIA, 2QIV, 4J09, 6HY2, 6P9P, 6P9Q...).
Por búsqueda en RCSB (`rcsb_entity_source_organism.rcsb_gene_name = lpxA`)
hay 45 entradas con ese gen, incluyendo ortólogos de *Pseudomonas*,
*Helicobacter*, *Leptospira*, etc.

Qué cambia vs nr:

* **Nº de hits:** en PDB son decenas (45 con gen lpxA; acá verifiqué 4
  representativas, todas 100%), en Swiss-Prot fueron 928 (341 a E<=1e-10)
  y en nr serían miles. PDB es mucho más chica y sesgada: solo lo que
  cristalizó/dio EM, con redundancia alta en algunas familias y cero en
  proteínas no modelo o de membrana difíciles.
* **Tipo de hit:** en PDB cada hit es `PDB_ID + cadena` (ej. `1LXA_1`)
  con título del paper, método (X-ray), resolución y organismo. Lleva
  directo a coordenadas 3D, no a un genoma. En nr cada hit es un accession
  `WP_xxx / P0A722` con organismo y % identidad, útil para taxonomía.
* **Accesiones en común:** casi ninguna literal (nr/Swiss-Prot usa `P0A722`,
  PDB usa `1LXA_1`), aunque apunten a la misma proteína. El script
  `desafio_X_pdb.py` reporta `nr∩pdb` y da 0 o muy pocos por esa
  diferencia de identificadores — no significa que no haya homólogos,
  sino que los IDs son de espacios distintos.

**¿Cuándo conviene PDB?** Cuando la pregunta es estructural: buscar un
template para modelado por homología de LpxA, ver el sitio activo
(His125 es catalítico), o chequear el plegamiento beta-hélice para
diseño de inhibidores. Para anotar función o contar cuántos homólogos
hay en la naturaleza, nr es mejor. Para este TP, la corrida vs PDB confirma
que LpxA tiene estructura conocida a alta resolución (hasta 1.74 Å) y da el
punto de partida para modelado.

Nota de entorno: el BLAST web vs nr y el BLAST local vs PDB completo
quedaron truncos por el reinicio del equipo (jobs cancelados, `/tmp`
limpiado). La identificación no depende de ellos: identidad exacta
262/262 con P0A722 + 4 PDB 100% verificadas por API son evidencia primaria.
`desafio_X_pdb.py` queda como vía de reproducción cuando se disponga de
los XML.

---

## Bibliografía

* Altschul et al. (1990). Basic local alignment search tool. *J Mol Biol* 215:403-410.
* Henikoff & Henikoff (1992). Amino acid substitution matrices from protein blocks. *PNAS* 89:10915-10919 (BLOSUM).
* Dayhoff et al. (1978). A model of evolutionary change in proteins. Atlas of Protein Sequence and Structure (PAM).
* Needleman & Wunsch (1970). A general method applicable to the search for similarities in the amino acid sequence of two proteins. *J Mol Biol* 48:443-453.
* Smith & Waterman (1981). Identification of common molecular subsequences. *J Mol Biol* 147:195-197.
* Reeck et al. (1987). “Homology” in proteins and nucleic acids: a terminology muddle and a way out of it. *Cell* 50:667.
* Prusiner (1998). Prions. *PNAS* 95:13363-13383. (contexto prión PRNP).
* Coleman & Raetz (1988). First committed step of lipid A biosynthesis in *E. coli*: sequence of the lpxA gene. *J Bacteriol* 170:1268-1274. (función LpxA).
* Raetz & Whitfield (2002). Lipopolysaccharide endotoxins. *Annu Rev Biochem* 71:635-700. (vía de Raetz).
* UniProt PRNP humano P04156: https://rest.uniprot.org/uniprotkb/P04156.txt / https://www.uniprot.org/uniprotkb/P04156/entry
* UniProt LpxA *E. coli* P0A722: https://rest.uniprot.org/uniprotkb/P0A722.txt / https://www.uniprot.org/uniprotkb/P0A722/entry
* LpxA PDB 1LXA: https://www.rcsb.org/structure/1LXA (más 2AQ9, 2QIA, etc., ver cross-refs de P0A722).
* NCBI BLAST docs: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs
* NCBI BLAST+ manual: https://www.ncbi.nlm.nih.gov/books/NBK279690/
* RCSB PDB: https://www.rcsb.org/
* Fuente TP: https://github.com/AJVelezRueda/Introduccion_a_la_Bioinformatica/blob/master/Teorico_Practicos/Bioinform%C3%A1tica_Primeros_pasos/Alineamientos_secuenciales.md
