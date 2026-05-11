# Interactions Search

Detecta y clasifica interacciones no covalentes entre un ligando y una proteína a partir de archivos PDB. Identifica puentes de hidrógeno, interacciones aromáticas (π-π y T-shaped) e interacciones hidrofóbicas, valida cada contacto por distancia y ángulo, y genera salidas en CSV y scripts TCL para visualización en VMD.

---

## Archivos

| Archivo | Descripción |
|---|---|
| `Interactions_search.py` | Script principal |
| `Interacciones_variables.yml` | Parámetros de distancia, aceptores y dadores por residuo |
| `Geometry.py` | Módulo auxiliar de geometría (referencia, no se importa directamente) |

---

## Dependencias

```
biopython
rdkit
pandas
numpy
pyyaml
```

---

## Uso

### Modo 1 — PDB separados

```bash
python Interactions_search.py -r proteina.pdb -l ligando.pdb -c A
```

### Modo 2 — Batch (múltiples ligandos, un receptor)

```bash
python Interactions_search.py -r proteina.pdb -l lig1.pdb lig2.pdb lig3.pdb -c A
```

Cada par genera su propia carpeta de salida. Los CSV acumulativos (`Interactions_close.csv`, `CM_all.csv`) se appendean automáticamente.

### Modo 3 — PDB complejo

Se tiene un solo PDB con proteína + ligando(s). El script lo separa automáticamente.

```bash
# Un solo HETATM (se selecciona automáticamente)
python Interactions_search.py -x complejo.pdb -c A

# Varios HETATM (se elige con -n)
python Interactions_search.py -x complejo.pdb -c A -n LIG
```

Si hay múltiples grupos HETATM y no se usa `-n`, el script lista los disponibles y sale sin analizar.

### Argumentos

| Argumento | Descripción |
|---|---|
| `-x / --complex` | PDB complejo. Alternativo a `-r`. |
| `-r / --receptor_pdb` | PDB del receptor ya separado. |
| `-l / --ligand_input` | PDB(s) del ligando. Acepta uno o varios (batch). |
| `-c / --chain_receptor` | Cadena de la proteína (ej: `A`). |
| `-n / --lig_name` | Nombre del HETATM cuando hay varios en `--complex`. |

---

## Pipeline de análisis

```
PDB complejo (opcional)
        │
        ▼
[1] split_pdb()
    ├── <stem>_protein.pdb      ← registros ATOM
    └── <stem>_RESNAME.pdb      ← registros HETATM por residuo (sin agua)
        │
        ▼
[2] Limpieza del ligando
    └── remove_bias()           ← elimina átomos CM del PDB del ligando
        │
        ▼
[3] Hot-points del ligando  (RDKit + SMARTS)
    ├── Aceptores de H-bond:  [O;H1], [O;H0], [N;H1], [N;H0], [n], [o], [N+]
    ├── Dadores de H-bond:    [O;H], [N;H2], [N;H], [S;H], [nH]
    └── Anillos aromáticos:   detección por ring_info, filtro por tamaño > 5
        │
        ▼
[4] Sitio activo del receptor  (BioPython)
    └── active_site_residues()
        Residuos cuyo centro de masa está a menos de 12 Å del CM del ligando
        Se excluyen HOH y el propio ligando
        │
        ▼
[5] Puntos de interés del receptor
    └── Coordenadas_interes_receptor()
        ├── Aceptores del receptor  (según tabla en YAML)
        ├── Dadores del receptor    (según tabla en YAML)
        └── Centros de anillos aromáticos (TYR, PHE, TRP)
        │
        ▼
[6] Búsqueda de contactos  (numpy, distancias vectorizadas)
    ├── H-bond:        lig aceptor ↔ rec dador     (umbral: Distances_Hidrogen_Bonds)
    │                  lig dador   ↔ rec aceptor    (umbral: Distances_Hidrogen_Bonds)
    ├── Aromática:     centroide ↔ centroide         (umbral: Distances_Aromatic)
    ├── Hidrofóbica:   C apolar lig ↔ C apolar rec   (umbral: Distances_Hidrofobica)
    ├── Salt bridge:   grupo ± lig ↔ grupo ∓ rec     (umbral: 4.0 Å)
    └── π-catión:      anillo lig ↔ ARG/LYS/HIS rec  (umbral: 5.0 Å)
        │
        ▼
[7] Validación por ángulo
    ├── H-bond:      ángulo D-A···Antecedente  entre 100° y 200°
    ├── Aromática:   ángulo entre planos de los anillos
    │                0°–30°  → π-π (paralela / sandwich)
    │                60°–90° → T-shaped (perpendicular)
    └── Hid. / salt / π-cat:  validadas solo por distancia
        │
        ▼
[8] Salidas
    ├── CSV con todas las interacciones brutas
    ├── CSV filtrado por distancia
    ├── CSV filtrado por distancia + ángulo  (interacciones validadas)
    ├── Resumen en consola  (tabla de interacciones validadas)
    ├── Resumen acumulativo Interactions_close.csv
    ├── Centros de masa     CM_all.csv
    └── Script TCL para VMD (si vmd_output: Yes)
```

---

## Separación de PDB (`split_pdb`)

Cuando se usa `--complex`, la función `split_pdb()` procesa el PDB antes del análisis:

- Separa registros `ATOM` (proteína) de `HETATM` (ligandos/cofactores)
- Agrupa los HETATM por nombre de residuo (`resName`)
- **Excluye agua** automáticamente: HOH, WAT, TIP, TIP3, SOL, DOD
- Distribuye los registros `CONECT` al archivo HETATM correspondiente (según serial de átomo)
- Guarda los archivos en `<stem>_split/`

Ejemplo con un PDB que tiene proteína + ligando LIG + grupo HEM:

```
complejo_split/
├── complejo_protein.pdb   ← toda la proteína
├── complejo_LIG.pdb       ← ligando orgánico + sus CONECT
└── complejo_HEM.pdb       ← grupo hemo
```

---

## Configuración (`Interacciones_variables.yml`)

```yaml
options:
  ligand_plot: 'Yes'    # genera imágenes PNG de aceptores, dadores y anillos del ligando
  vmd_output:  'Yes'    # genera script TCL para visualizar en VMD

distancias:
  Distances_Hidrogen_Bonds: 3.2   # Å — umbral para puentes de H
  Distances_Aromatic:       5.5   # Å — umbral centro-a-centro para aromáticas
  Distances_Hidrofobica:    4.0   # Å — umbral hidrofóbico (alineado con PLIP)
  centroid_distance:        9.0   # Å — radio de búsqueda de sitio activo (referencia)
  Distances_C_Simple:       1.54  # Å — enlace C-C simple (referencia)
  Distances_C_Doble:        2.56  # Å — enlace C=C doble (referencia)

acceptors:             # átomos aceptores por residuo
  ALA: [O]
  TYR: [O, OH]
  ASP: [O, OD1, OD2]
  ...

donors:                # átomos dadores por residuo
  ALA: [N]
  ARG: [N, HNE, HH11, HH12, HH21, HH22, HE, NE, NH1, NH2]
  ...

acceptors_antecedent:  # átomo antecedente del aceptor (para cálculo de ángulo)
  TYR: {OH: CZ}
  ASP: {OD1: CG, OD2: CG}
  ...

special:               # casos especiales (ej: grupo hemo)
  HEM: [FE, 1.59]
```

---

## Salidas

### Archivos CSV

| Archivo | Contenido |
|---|---|
| `<folder>/Interaction_<rec>_<lig>_all.csv` | Todas las interacciones encontradas (sin filtros) |
| `<folder>/Interaction_<rec>_<lig>_threshold.csv` | Filtradas por distancia |
| `<folder>/Interaction_<rec>_<lig>_true.csv` | Validadas por distancia y ángulo |
| `Interactions_close.csv` | Resumen acumulativo por corrida |
| `Interactions_all_count.csv` | Conteo por tipo (acceptor / donor / aromatic) |
| `CM_all.csv` | Centro de masa del ligando por corrida |

Columnas de los CSV de interacciones:

| Columna | Descripción |
|---|---|
| `Pos R` | Número de residuo del receptor |
| `Res` | Nombre del residuo (ej: SER, TYR) |
| `Atom` | Átomo involucrado del receptor |
| `Dist` | Distancia en Å |
| `Lig` | Átomo o anillo del ligando |
| `Type` | Tipo: `acceptor`, `donor`, `aromatic` |
| `Angle` | Ángulo de validación en grados |
| `Interaction` | `Yes` / `No` — si cumple criterios de distancia y ángulo |

### Script VMD (`vmd_<receptor>_<ligando>.tcl`)

Genera una visualización lista para cargar en VMD con:
- Proteína completa en NewCartoon transparente
- Residuos del sitio activo en Licorice
- Ligando en Licorice
- Líneas punteadas para cada interacción validada:
  - **Blanco** — aromáticas
  - **Rojo** — ligando aceptor (receptor dador)
  - **Amarillo** — ligando dador (receptor aceptor)
- Etiqueta con la distancia en Ångströms sobre cada línea

### Imágenes PNG del ligando (si `ligand_plot: Yes`)

| Archivo | Contenido |
|---|---|
| `<lig>_acceptors.png` | Ligando con átomos aceptores resaltados |
| `<lig>_donors.png` | Ligando con átomos dadores resaltados |
| `<lig>_aromatic.png` | Ligando con anillos aromáticos resaltados |

---

## Estructura de carpetas generada

Todo queda dentro de una única carpeta por par `<receptor>_<ligando>/`:

```
<receptor>_<ligando>/
├── <receptor>.pdb             ← copia del PDB del receptor
├── <ligando>.pdb              ← copia del PDB del ligando
├── <ligando>_old.pdb          ← copia pre-limpieza (remove_bias)
├── Interaction_*_all.csv      ← todas las interacciones sin filtro
├── Interaction_*_threshold.csv← filtradas por distancia
├── Interaction_*_true.csv     ← validadas por distancia + ángulo
├── summary.csv                ← conteo por tipo de interacción
├── CM.csv                     ← centro de masa del ligando
├── vmd_*.tcl                  ← script VMD (si vmd_output: Yes)
├── *_acceptors.png            ← ligando con aceptores resaltados
├── *_donors.png               ← ligando con dadores resaltados
└── *_aromatic.png             ← ligando con anillos resaltados
```

En modo batch cada par genera su propia carpeta independiente.

---

## Notas

- El script debe correrse desde el directorio donde están los PDB, o bien usar rutas absolutas.
- Para analizar múltiples ligandos en batch, se puede llamar el script en un loop de shell; los CSV acumulativos (`Interactions_close.csv`, etc.) se van appendeando automáticamente.
- Los residuos no estándar que no figuren en `acceptors` / `donors` del YAML se omiten sin error.
- Los anillos aromáticos del ligando deben tener más de 5 átomos para ser considerados (filtra ciclopentano y similares no aromáticos).
