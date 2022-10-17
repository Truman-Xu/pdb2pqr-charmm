[![Documentation Status](https://readthedocs.org/projects/pdb2pqr/badge/?version=latest)](https://pdb2pqr.readthedocs.io/en/latest/?badge=latest)

# PDB2PQR-CHARMM

This package contains the PDB2PQR software.  For more information, please see

* Home page:  http://www.poissonboltzmann.org/
* Documentation: http://pdb2pqr.readthedocs.io


This is a fork of PDB2PQR that address the patch problems when generating CHARMM compatible PDB files. To use the CHARMM specific protonation routine, the protein's chains and ligands have to be separated. This will be automated in the future pyCHARMM releases.
Currently, there are three major modifications from the original repo:
1. Added the patch for neutral lysine and its parameters.
2. Added script and `charmm_biomol` object that automagically protonates and generates CHARMM compatible PDB files.
3. Revert the patch names (PRES) back to their residue names (RESI) before generating pdb strings, as CHARMM/pyCHARMM **does not** read the patch name directly from pdb files. Rather, the patch command is needed for CHARMM to recognize those patches and to generate structures accordingly. The residue(s) need to be patched are also printed out in the PDB header or accessible from the `charmm_biomol` object.

## Installation
You will need all the dependencies for PDB2PQR plus pandas.
Clone this repo, and cd into it.
```
git clone https://github.com/Truman-Xu/pdb2pqr-charmm.git
cd pdb2pqr-charmm
python setup.py install
```
## Usage
To invoke the CHARMM specific routine for protein protonation,
```
python -m pdb2pqr.charmm_protonator [-h] [-x PREFIX] [-o OUTPATH] [--ph PH]
                                    [--save-df] [--report-propka]
                                    I

positional arguments:
  I                     Input filepath of .pdb file

options:
  -h, --help            show this help message and exit
  -x PREFIX, --prefix PREFIX
                        Prefix of the output .pdb files default to
                        "in_file_prefix_protonated"
  -o OUTPATH, --outpath OUTPATH
                        Output path of the processed .pdb file. default = "./"
  --ph PH               The pH value under which the protein's protonation
                        state is determined. Default = 7.0
  --save-df             Flag to save propka results dataframe
  --report-propka       Flag to report propka results
```
Example of using the `charmm_biomol` object is in the jupyter notebook `CHARMM-Prot.ipynb`
## Known Issues
The terminal patches PDB2PQR applies have additional hydrogen and oxygens. They will not be recognized by CHARMM



