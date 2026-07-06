# CppLibAmberGromacs — Knowledge Base

**Authors:** Ezequiel Cuenca, Nicolás Loubet
**Language:** C++17

---

## 1. Purpose

Header-only library for reading and analyzing molecular dynamics configurations from **AMBER**, **GROMACS**, and **LAMMPS**. The typical workflow is:

1. Read topology → `TopolInfo`
2. Create a coordinate reader
3. Build `Configuration` objects from coordinate files + topology
4. Compute per-molecule quantities (V4S, interactions per site, defects, order parameters...)

---

## 2. Build Flags

Define before including the main header:

| Flag | Effect |
|---|---|
| `NOHUP` | Progress bar prints line-by-line instead of in-place (for log files) |
| `CONSOLE_WIDTH` | Integer, controls progress bar width (default 150) |
| `USE_FAST_READER` | AMBER only: fast sequential PDB reader. Assumes perfectly ordered PDB where each ATOM...TER block = one molecule in topology order |
| `USE_VECTOR_TOPOLOGY` | Stores per-molecule atom data as `vector` instead of `map`. Reduces RAM for large systems. **Incompatible with LAMMPS** (excludes LAMMPS reader from build) |

```cpp
#define NOHUP
// #define USE_VECTOR_TOPOLOGY
#include <CppLibAmberGromacs.hpp>
```

---

## 3. Type Aliases

| Name | Underlying type |
|---|---|
| `Real` | `float` or `double` (defined globally, typically `float`) |
| `NOT_CLASSIFIED` | `#define NOT_CLASSIFIED -314159265` — sentinel for unclassified molecules |

---

## 4. Class Hierarchy

```
Particle
└── Atom
└── Molecule
    └── Water

CoordinateReader  (abstract)
├── AmberCoordinateReader
├── GromacsCoordinateReader
└── LammpsCoordinateReader

TopologyReader  (abstract)
├── AmberTopologyReader
├── GromacsTopologyReader
└── LammpsTopologyReader

Configuration
├── ConfigurationBulk       (water bulk systems)
├── ConfigurationLipid      (lipid membrane systems)
├── ConfigurationDefects    (defect classification)
└── ConfigurationDiscriminated  (per-atom-name energy breakdown)
```

---

## 5. Reading a System

### 5.1 Factory Pattern

```cpp
// Topology
TopolInfo ti = ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::AMBER)
               ->readTopology("path/to/step5_input.parm7");

// Coordinate reader (reuse across frames)
auto* reader = ReaderFactory::createCoordinateReader(ReaderFactory::ProgramFormat::AMBER);

// Configuration (one per frame)
Configuration conf(reader, "path/to/frame1.pdb", ti);
```

Supported formats: `AMBER`, `GROMACS`, `LAMMPS`.

### 5.2 File Types by Format

| Format | Topology file | Coordinate file |
|---|---|---|
| AMBER | `.parm7` / `.prmtop` | `.pdb` |
| GROMACS | `.top` | `.gro` |
| LAMMPS | `.data` | `.xyz` |

### 5.3 TopolInfo — Key Fields

```cpp
struct TopolInfo {
    int num_molecules;      // total molecules
    int num_solutes;        // non-solvent molecules
    int num_solvents;       // water/solvent molecules
    int total_number_of_atoms;

    // Per-molecule atom data: [mol_idx] -> {type, name, charge, mass}
    // Key is 1-based (int) without USE_VECTOR_TOPOLOGY, 0-based index with it
    vector<map<int,tuple<string,string,Real,Real>>> atom_type_name_charge_mass;
    // or vector<vector<tuple<...>>> with USE_VECTOR_TOPOLOGY

    map<string,pair<Real,Real>> type_LJparam;      // type -> {epsilon [kJ/mol], sigma [Å]}
    map<pair<string,string>,pair<Real,Real>> special_interaction; // (t1,t2) -> {eps, sig}
    map<string,int>    type_Z;                     // type -> atomic number
    map<string,string> name_type;                  // atom_name -> type
    map<string,int>    number_of_each_different_molecule;
    map<string,int>    number_of_atoms_per_different_molecule;
    Vector             default_system_bounds;
};
```

**AMBER charge convention:** charges in the `.parm7` are stored in AMBER units (×18.2223). The reader converts them to elementary charge units automatically.

---

## 6. Configuration

### 6.1 Construction

```cpp
Configuration conf(reader, filename, ti);
```

`reader` can be reused across frames — it is stateless for AMBER/GROMACS. `ti` must remain valid for the lifetime of `conf`.

### 6.2 Core Getters

```cpp
int     conf.getNMolec();            // total number of molecules
Vector  conf.getBounds();            // box dimensions in Å
const Molecule& conf.getMolec(int id);  // 1-indexed
Molecule*       conf.getMolec_ptr(int id);
const Water&    conf.getWater(int id);  // throws if not water
```

### 6.3 Molecule and Atom Access

```cpp
const Molecule& mol = conf.getMolec(m);
int    mol.getNAtoms();
bool   mol.isWater();
int    mol.getClassification();      // NOT_CLASSIFIED if not set
Vector mol.getPosition();            // center of mass (or O for Water)
Real   mol.getCharge();              // total molecular charge
Real   mol.getMass();                // total molecular mass

const Atom& atom = mol.getAtom(a);   // 1-indexed
Vector atom.getPosition();
Real   atom.getCharge();
Real   atom.getMass();
Real   atom.getSigma();
Real   atom.getEpsilon();
int    atom.getZ();                  // atomic number
string atom.getAtomType();
```

For **Water** specifically:
```cpp
const Water& w = *static_cast<Water*>(conf.getMolec_ptr(m));
w.getOxygen();      // atoms[0]
w.getHydrogen_1();  // atoms[1]
w.getHydrogen_2();  // atoms[2]
// TIP4P systems also have MW at atoms[3], part of the same molecule
```

---

## 7. V4S — Interactions Per Site

The central quantity. Assigns each neighbour to one of four tetrahedral sites around a water molecule and accumulates the interaction energy. Returns the 4 values **sorted ascending** (V1S ≤ V2S ≤ V3S ≤ V4S).

### 7.1 Main Overloads

```cpp
// Standard: returns 4 sorted site energies [kJ/mol]
vector<Real> conf.getInteractionsPerSite(int ID, Real R_CUT_OFF = 5.0, int* labels = nullptr);

// With water-water flag: flag_ww=true if ALL 4 sites have ≥1 water interaction ≤ V_CUT_OFF
vector<Real> conf.getInteractionsPerSite(int ID, bool& flag_ww, Real V_CUT_OFF = -12.0, Real R_CUT_OFF = 5.0);

// Returns {full_sorted, water_only_sorted} pair
pair<vector<Real>,vector<Real>> conf.getInteractionsPerSite_waterOnly(
    int ID, Real V_CUT_OFF = -12.0, Real R_CUT_OFF = 5.0);

// Single site value: i_V-th element (1-based), equivalent to getInteractionsPerSite(...)[i_V-1]
Real conf.v_4S(int ID, int i_V = 4, Real R_CUT_OFF = 5.0);
```

### 7.2 Site Assignment Logic

- Tetrahedral sites are computed from the O, H1, H2 positions of the center water using `Geometrics::getPerfectTetrahedron`.
- Each neighbour (water oxygen, solute atom, or ion) is assigned to the **closest** site within `R_CUT_OFF`.
- Water neighbours are pre-screened by `O-O distance > R_CUT_OFF + 1.1 Å` (skip if too far).
- Solute molecules are pre-screened by `molecule_position-to-molecule_position distance`, then atom-by-atom within `R_CUT_OFF + 1.1 Å`.
- Ions (1 atom, |charge| ≥ 1) are always considered without pre-screening.

### 7.3 Water VList (raw sorted potentials, no site assignment)

```cpp
vector<Real> conf.getVList(int ID, Real MAX_V4 = 5.5);  // all W-W potentials, sorted ascending
Real         conf.vI(int ID, int V_index = 4);          // V_index-th smallest (1-based)
```

---

## 8. Typical Usage Patterns

### 8.1 GROMACS Bulk Water Analysis

```cpp
TopolInfo ti = ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::GROMACS)->readTopology("system.top");
auto* reader = ReaderFactory::createCoordinateReader(ReaderFactory::ProgramFormat::GROMACS);

for(int i_conf = 0; i_conf <= N_CONF; i_conf++) {
    if(!fs::exists("em-" + to_string(i_conf) + ".gro")) continue;
    Configuration conf(reader, "em-" + to_string(i_conf) + ".gro", ti);

    for(int m = 1; m <= conf.getNMolec(); m++) {
        if(!conf.getMolec(m).isWater()) continue;
        vector<Real> vis = conf.getInteractionsPerSite(m, /*R_CUT_OFF=*/ 7.0); // vis[0]=V1S ... vis[3]=V4S
        bool is_D = conf.vI(m) > -12.0;
        // ... accumulate
    }
}
```

### 8.2 AMBER Lipid/Water Interfacial Analysis

```cpp
TopolInfo ti = ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::AMBER)->readTopology(path + "step5_input.parm7");
auto* reader = ReaderFactory::createCoordinateReader(ReaderFactory::ProgramFormat::AMBER);

for(int i_conf = 1; i_conf <= N_CONF; i_conf++) {
    Configuration conf_RD(reader, path + "RD/frame" + to_string(i_conf) + ".pdb", ti);
    Configuration conf_IS(reader, path + "IS/frame_min_" + to_string(i_conf) + ".pdb", ti);

    // Center-of-mass z-reference from solute
    Real z_ref = calculateCenterOfMass(conf_RD, ti.num_solutes);

    // Iterate over solvent molecules only
    for(int m = ti.num_solutes + 1; m <= conf_RD.getNMolec(); m++) {
        Real z = fabs(conf_RD.getMolec(m).getPosition().z - z_ref);
        if(z > MAX_BINS || z < MIN_BINS) continue;
        vector<Real> vis = conf_IS.getInteractionsPerSite(m);
        // ... bin by z
    }
}
```

**Note:** `ti.num_solutes` gives the count of non-water molecules. Solvent molecules are indexed `ti.num_solutes + 1 .. num_molecules`.

### 8.3 Molecule Classification (ConfigurationBulk)

```cpp
ConfigurationBulk conf(reader, filename, ti);

// D/T0/T1/T2 scheme (JCP 2023)
conf.classifyMolecules(/*V_index=*/4, /*threshold=*/-12.0);

// D3/D5/TA/TB scheme (PRE 2024)
conf.classifyMolecules_includePentacoordinated(4, -12.0);

// Per-molecule checks (also cache the classification on the Water object)
bool d  = conf.isD(m);
bool d3 = conf.isD3(m);
bool d5 = conf.isD5(m);
bool dx = conf.isDX(m);  // D3 or D5

// Tanaka ζ parameter
Water* w = static_cast<Water*>(conf.getMolec_ptr(m));
Real zeta = conf.Tanaka(w);

// Local Structure Index
Real lsi = conf.LSI(m);
```

Classification constants (both schemes share same integer values):
```cpp
ConfigurationBulk::CLASSIFICATION_D_MOLECULE  = 0
ConfigurationBulk::CLASSIFICATION_T0_MOLECULE = 1
ConfigurationBulk::CLASSIFICATION_T1_MOLECULE = 2
ConfigurationBulk::CLASSIFICATION_T2_MOLECULE = 3
// Pentacoordinated scheme aliases:
ConfigurationBulk::CLASSIFICATION_D3_MOLECULE = 0
ConfigurationBulk::CLASSIFICATION_D5_MOLECULE = 1
ConfigurationBulk::CLASSIFICATION_TA_MOLECULE = 2
ConfigurationBulk::CLASSIFICATION_TB_MOLECULE = 3
```

### 8.4 Defect Analysis (ConfigurationDefects)

```cpp
ConfigurationDefects conf(reader, filename, ti);
DefectInfo info = conf.classifyDefect(m, /*R_CUT_OFF=*/5.0, /*V_CUT_OFF=*/-12.0);

info.is_D3;  // under-coordinated defect
info.is_D5;  // over-coordinated defect
info.is_DJ;  // simultaneously D3 and D5
info.lacking_sites;       // count of empty sites
info.bifurcated_sites;    // count of sites with ≥2 neighbours
info.sum_per_site;        // 4-element sorted vector
info.lacking_site_position;
info.bifurcated_site_position;
info.bifurcated_individual_potentials;  // pair<Real,Real>
info.bifurcated_individual_distances;
info.bifurcated_individual_indices;     // molecule IDs
```

### 8.5 Lipid Order Parameter (ConfigurationLipid)

```cpp
ConfigurationLipid conf(reader, filename, ti);

// Returns [3][22] matrix: [chain_index (SN1/SN2)][carbon_position] -> order parameter
vector<vector<Real>> op = conf.orderParameter(mol_id);
// op[0] = SN1 chain, op[1] = SN2 chain, op[2] = (unused third slot)
```

---

## 9. Neighbours and Distances

```cpp
// All molecule IDs within D_MAX_NEI of ID_CENTER
vector<int> nearby = conf.findNearby(int ID_CENTER, Real D_MAX_NEI);

// Intra-molecule atom neighbours (for ConfigurationLipid chain tracing)
vector<int> near_atoms = mol.findNearbyAtoms(int ID_CENTER, Real D_MAX, const Vector& bounds);

// Particle distances (PBC-aware)
Real d = particleA.distanceTo(particleB, bounds);
Real d = distancePBC(vectorA, vectorB, box);
```

---

## 10. Parallelism (ToolKit)

```cpp
// Function signature required by parallel/serial:
// void processSystem(int id_thread, int n_threads, T entry, Res& output, Args... extra_args)

vector<T>   list   = { ... };
vector<Res> output(list.size());

ToolKit::parallel(processSystem, list, output, /*max_threads=*/26);
ToolKit::serial(processSystem, list, output);

// Wrap everything in a timer:
ToolKit::takeTime([] {
    // ... all work here
});
```

Inside `processSystem`, use `id_thread` and `n_threads` for the progress bar:
```cpp
ToolKit::printPercentageBar(i_conf, N_CONF, id_thread, n_threads);
```

---

## 11. Histogram Utilities (ToolKit)

```cpp
// Returns 0-based bin index for value in [LIMIT_MIN, LIMIT_MAX) split into N_BINS bins
int bin = ToolKit::getBinPosition(Real value, Real LIMIT_MIN, Real LIMIT_MAX,
                                   int N_BINS, bool trim = false);
// trim=true: clamps to [0, N_BINS-1] instead of throwing on out-of-range
```

---

## 12. Sorter

```cpp
Sorter::sort(vector<T>& v);                              // ascending
Sorter::sort(vector<T>& v, Sorter::Order::Descending);

// Co-sort: sort `values` and apply same permutation to `indexes`
Sorter::cosort(vector<T1>& values, vector<T2>& indexes);
Sorter::cosort(vector<T1>& values, vector<T2>& indexes, Sorter::Order::Ascending);
```

---

## 13. Geometrics

```cpp
// Tetrahedral sites (used internally by getInteractionsPerSite)
TetrahedronVertices tv = Geometrics::getPerfectTetrahedron(O_pos, H1_pos, H2_pos, bounds);
vector<Vector> sites = tv.toVector();  // order: H1, H2, L1, L2
// H1/H2 = sites near the hydrogens; L1/L2 = lone-pair sites

// Plane, line, box utilities
Plane p = Geometrics::getPlaneFromPoints(a, b, c, bounds);
Real  d = Geometrics::distanceToPlane(plane, point);
Real  d = Geometrics::distanceToLine(line, point, bounds);
bool inside = Geometrics::isInBox(pos, box);
bool inside = Geometrics::isInBox_2D(pos, box);
```

---

## 14. Physical Constants and Units

All distances in **Å**, energies in **kJ/mol**, charges in **elementary charge (e)**.

```cpp
Real K_COULOMB = 1389.35458;  // kJ/mol per e² per Å  (in Atom.hpp)
Real k_B       = 8.314e-3;    // kJ/(mol·K)           (in ToolKit.hpp)
```

AMBER `.parm7` charges are stored in AMBER units (divide by 18.2223 to get e). The reader does this automatically.

GROMACS `.top` sigmas are in nm; the reader multiplies by 10 to convert to Å.

---

## 15. Known Constraints and Gotchas

- **Molecule indexing is 1-based** everywhere (`getMolec(1)` = first molecule, `getAtom(1)` = first atom).
- **`getInteractionsPerSite` only works on water molecules.** Throws `invalid_argument` if called on a non-water molecule.
- **`USE_VECTOR_TOPOLOGY` is incompatible with LAMMPS.** Enabling it excludes the LAMMPS reader from the build.
- **`classifyMolecules` / `isD` / `isD3`** cache classification on the `Water` object. Calling `isD` after `classifyMolecules` returns the cached value, not a recomputation.
- **AMBER PDB reader:** by default uses a flexible reader that handles multi-residue molecules. Use `USE_FAST_READER` only for clean AMBER-generated PDBs where each ATOM…TER block maps exactly to one topology molecule.
- **TIP4P water:** the MW virtual site (index 3, 0-based, within each WAT molecule) has charge −1.1128 e and zero LJ parameters. It is part of the same `Water` object as OW/HW1/HW2 and is included automatically in Coulomb calculations inside `potentialWith`.
- **`NOT_CLASSIFIED` = −314159265**: sentinel used as default classification. Check `mol.isClassified()` before reading `mol.getClassification()`.
- **Center-of-mass vs oxygen position:** `Molecule::getPosition()` returns the geometric center of mass (average of atom positions, not mass-weighted). For `Water`, the position is explicitly set to `atoms[0].getPosition()` (the oxygen). For general molecules, it is the unweighted centroid.
