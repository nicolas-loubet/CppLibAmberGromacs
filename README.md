# CppLibAmberGromacs
Library for handling AMBER, GROMACS and LAMMPS molecular dynamics files.

### Version: 2.1.8

## Authors

- Ezequiel Cuenca
- Nicolás Loubet

## Folder Structure

- `files`: Contains example files for AMBER, GROMACS and LAMMPS.

- `cppambergromacs`: Contains the .hpp files and related components for studying molecular dynamics. Refer to the figure in the `files` folder to understand the class hierarchy and relationships.

- `tests`: Contains .cpp files used for testing specific functionalities.

## Optional commands:

- NOHUP: Changes the visualization of progress bar to a normal text with line-break.

- CONSOLE_WIDTH: Controls the max length of the progress bar

- USE_FAST_READER: Uses a faster, lighter PDB coordinate reader that reads molecules sequentially. Assumes the PDB is perfectly ordered and each ATOM...TER block maps exactly to one molecule in topology order. Use only if your PDB was generated cleanly by AMBER. Default (without this flag) uses a slower but more flexible reader that handles PDBs where molecules are split across multiple residues.

- USE_VECTOR_TOPOLOGY: Stores per-molecule atom data as `vector` instead of `map` during topology reading. Dramatically reduces RAM usage for large systems. Incompatible with LAMMPS (enabling this flag excludes the LAMMPS reader from the build).
