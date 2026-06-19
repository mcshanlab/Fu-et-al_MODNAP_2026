# MODNAP Benchmark

## Introduction

The MODNAP (Modified Nucleic Acids & Nucleic Acid–Protein Complexes) benchmark set is an open-source resource for evaluating structure prediction and molecular modeling of modified nucleic acid and nucleic acid-protein complexes. It is also a resource of modified nucleotides and nucleosides present in structures deposited to the Protein Data Bank. MODNAP was generated from atomic coordinates deposited in the RCSB Protein Data Bank (PDB): https://www.rcsb.org/.

Key functions of this repository:

- MODNAP dataset: High-quality  nucleic acid and nucleic acid-protein structure benchmarking datasets containing modified nucleotides.

- Raw data and code for reproducing the manuscript results


## Dataset availability:


### Raw Data and Scripts Provided

- directory `chemminetools` contains raw data related to the ChemMine Tools analysis

- directory `list-of-modified-nucleotides-in-the-PDB` contains information on the modified nucleotide and nucleotides reported in the manuscript

- directory `MODNAP` contains MODNAP entries DPB files, AlphaFold 3 model files, information on the features of MODNAP, and manuscript raw data / analysis scripts for benchmarking effforts reported in the manuscript

- directory `posebusters` contains raw data and analysis scripts reported for the PoseBusters analysis reported in the manuscript

- directory `trees` contains Netwick tree files to generated dendrograms of modified nucleotides / nucleosides

- directory `example_AlphaFold3_run_CCD-77Y_PDB-6EO6` contains scripts for an example AlphaFold 3 on a SLURM high performance cluster


### MODNAP datasets

- The MODNAP benchmark dataset is located in the `MODNAP/MODNAP` directory.

- AlphaFold 3 models of MODNAP entires are available in the `MODNAP/MODNAP_AF3_models` directory.


## Prerequisites

Before starting, ensure you have the following installed on your system (depending on what you what to test) or available on online servers:

- ChemMine Tools: https://chemminetools.ucr.edu/ 
- PoseBusters: https://github.com/maabuu/posebusters
- AlphaFold3:  https://github.com/google-deepmind/alphafold3
- OpenStructure: https://openstructure.org/install
- PyMOL (command line version): https://www.pymol.org/
- Python 3 and associated dependencies: https://www.python.org/downloads/
