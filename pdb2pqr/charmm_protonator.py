from io import StringIO
from collections import namedtuple, OrderedDict

import pandas as pd
import os
import propka
import propka.input as pk_in
import propka.output as pk_out
from propka.parameters import Parameters
from propka.molecular_container import MolecularContainer

import pdb2pqr
from pdb2pqr import io
from pdb2pqr import hydrogens, debump, forcefield
from pdb2pqr.utilities import noninteger_charge
from pdb2pqr.main import is_repairable

propka_param = os.path.join(propka.__path__[0],'propka.cfg')

def run_propka(propka_param, biomolecule, options):
    """Run a PROPKA calculation.
    :param biomolecule:  biomolecule object
    :type biomolecule:  Biomolecule
    :return:  (DataFrame-convertible table of assigned pKa values,
               pKa information from PROPKA)
    :rtype:  (list, str)
    """

    lines = io.print_biomolecule_atoms(
        atomlist=biomolecule.atoms, chainflag=True, pdbfile=True
    )

    with StringIO() as fpdb:
        fpdb.writelines(lines)
        parameters = pk_in.read_parameter_file(propka_param, Parameters())
        molecule = MolecularContainer(parameters, options)
        # needs a mock name with .pdb extension to work with stream data, hence the "input.pdb"
        molecule = pk_in.read_molecule_file("input.pdb", molecule, fpdb)
    
    molecule.calculate_pka()

    lines = []

    lines.append(
        pk_out.get_folding_profile_section(
            molecule,
            conformation="AVR",
            reference="neutral",
            window=[0.0, 14.0, 1.0],
        )
    )
    lines.append(
        pk_out.get_charge_profile_section(molecule, conformation="AVR")
    )
    lines.append(pk_out.get_the_line())
    pka_str = "\n".join(lines)

    # Summarize in pKas in DataFrame for later use
    conformation = molecule.conformations["AVR"]
    rows = []
    for group in conformation.groups:
        row_dict = OrderedDict()
        atom = group.atom
        row_dict["res_num"] = atom.res_num
        row_dict["ins_code"] = atom.icode
        row_dict["res_name"] = atom.res_name
        row_dict["chain_id"] = atom.chain_id
        row_dict["group_label"] = group.label
        row_dict["group_type"] = getattr(group, "type", None)
        row_dict["pKa"] = group.pka_value
        row_dict["model_pKa"] = group.model_pka
        row_dict["buried"] = group.buried
        if group.coupled_titrating_group:
            row_dict["coupled_group"] = group.coupled_titrating_group.label
        else:
            row_dict["coupled_group"] = None
        rows.append(row_dict)
    return rows, pka_str

class charmm_biomol(pdb2pqr.biomolecule.Biomolecule):
    ## TODO: get rid of terminal patches from pdb2pqr from upstream
    def __init__(self, pdb_file):
        # PRES (patched residue ids) from CHARMM 
        # top_all36_prot.rtf
        self.patch_names = [
             'NTER','GLYP','PROP','ACE','ACED',
             'ACP','ACPD','NNEU','NGNE','CTER',
             'CNEU','PCTE','CT1','CT2','CT3',
             'ASPP','GLUP','LSN','CYSD','SERD',
             'DISU','HS2','LIG1','LINK'
        ]
        # Definitions contain reference structures for amino acids and nucleic acids
        definition = io.get_definitions()
        pdblist, is_cif = io.get_molecule(pdb_file)
        super().__init__(pdblist, definition)
        self.set_termini(neutraln=False, neutralc=False)
        self.update_bonds()
        self.forcefield_ = forcefield.Forcefield(
            'CHARMM', definition, userff = None, usernames = None
        )
        # Default options for propka
        options_def = namedtuple('options',['keep_protons', 'protonate_all', 'chains', 'titrate_only', 'display_coupled_residues'])
        self.propka_options = options_def(False, False, None, None, False)
        hydrogen_handler = hydrogens.create_handler()
        debumper = debump.Debump(self)
        self.hydrogen_routines = hydrogens.HydrogenRoutines(
            debumper, hydrogen_handler
        )
        self.is_protonated = False
        
    def prep_for_protonation(self):

        if is_repairable(self, has_ligand = False):
            print(
                f"Attempting to repair {self.num_missing_heavy:d} "
                "missing atoms in biomolecule."
            )
            self.repair_heavy()

        self.update_ss_bridges()
        self.remove_hydrogens()
        
    def protonate(self, ph, report = False):
        self.ph = ph
        pka_ordered_dict, self.pka_str = run_propka(propka_param, self, self.propka_options)
        if report:
            print(f"PROPKA information:\n{self.pka_str}")
            
        self.apply_pka_values(
            self.forcefield_.name,
            ph,
            {
                f"{row['res_name']} {row['res_num']} {row['chain_id']}": row[
                    "pKa"
                ]
                for row in pka_ordered_dict
                if row["group_label"].startswith(row["res_name"])
            },
        )
        
        self.add_hydrogens()
        self.is_protonated = True
        self.pka_df = pd.DataFrame(pka_ordered_dict)
    
    def optimize_hydrogens(self):
        
        self.hydrogen_routines.set_optimizeable_hydrogens()
        self.hold_residues(None)
        self.hydrogen_routines.initialize_full_optimization()

        self.hydrogen_routines.optimize_hydrogens()
        self.hydrogen_routines.cleanup()
        
    def apply_ffname(self):
        # set_state() is where the magic happens (patch names applied)
        self.set_states()
        # convert patch names to the forcefield specific ones 
        name_scheme = self.forcefield_
        self.apply_name_scheme(name_scheme)
        
        # and add attributes of radius and charge on atoms
        self.matched_atoms, self.missing_atoms = self.apply_force_field(self.forcefield_)
        # modifications to set the patch name back to resi names for CHARMM
        self.patched_res = set()
        self.pdb_strings = []
        
        for res in self.residues:
            # if patch attribute exists for the residue and it's not peptide
            if hasattr(res,'patches') and res.patches != ['PEPTIDE']:
                for atom in res.atoms:
                    if atom.res_name in self.patch_names:
                        # record which residues the patches have been applied
                        self.patched_res.add((res.res_seq, res.name, atom.res_name))
                        # revert the name back to standard CHARMM resi names
                        atom.res_name = res.name
                    elif atom.res_name == 'TER':
                        if not (atom.name.startswith('HT') or atom.name.startswith('OT')):
                            atom.res_name = res.name

    def check_charges(self):

        total_charge = 0
        for residue in self.residues:
            charge = residue.charge
            charge_err = noninteger_charge(charge)
            if charge_err:
                print(
                    f"Residue {residue} has non-integer charge:  {charge_err}"
                )
            total_charge += charge
        charge_err = noninteger_charge(total_charge)

        if charge_err:
            raise ValueError(charge_err)
        
    def auto_protonate(self, ph):
        self.prep_for_protonation()
        self.protonate(ph)
        self.optimize_hydrogens()
        self.apply_ffname()
        self.check_charges()
        
    def get_pdb_header(self):
        reslist, charge = self.charge
        if self.is_protonated:
            header = io.print_pqr_header(
                self.pdblist,
                self.missing_atoms,
                reslist,
                charge,
                'CHARMM',
                'propka',
                self.ph,
                'CHARMM',
                include_old_header=False,
            )
        else:
            header = io.get_old_header(self.pdb_list)
            print("No protonation performed. Original header used.")
            
        if self.patched_res:
            header+= 'REMARK   7 FOLLOWING PATCHE(S) NEED TO BE APPLIED IN CHARMM\n'
            for res in self.patched_res:
                header+= 'REMARK   7 RESID: {} RESI: {} PRES: {}\n'.format(*res)
        return header

    def get_pdb_str(self):
        pdb_strings = []
        for atom in self.matched_atoms:
            pdb_strings.append(atom.get_pdb_string())
        return pdb_strings
    
    def save_pdb(self, pdb_filename):
        header = self.get_pdb_header()
        with open(pdb_filename, 'w') as f:
            f.write(header)
            for l in self.get_pdb_str():
                # ignore the TERminal hydrogens and oxigens added by pdb2pqr
                if l[17:20] != 'TER':
                    f.write(l+'\n')
            f.write('TER\nEND\n')

if __name__ == '__main__':
    import argparse
    import os
    import sys
    parser = argparse.ArgumentParser(
        prog="CHARMM protein protonator",
        description='Protonate a protein with pdb2pqr and propka \
            and generate PDB file that \
            is compatible with CHARMM residue and atom naming scheme.'
        )

    parser.add_argument('inpath',metavar='I', type=str,
                        help='Input filepath of .pdb file')
    parser.add_argument('-x','--prefix', type=str,
                        help='Prefix of the output .pdb files\
                        default to "in_file_prefix_protonated"')
    parser.add_argument('-o','--outpath', type=str, default='./',
                        help='Output path of the processed\
                            .pdb file. default = "./"')
    parser.add_argument('--ph', type=float, default=7.0,
                        help="The pH value under which the protein's \
                            protonation state is determined. \
                            Default = 7.0") 
    parser.add_argument('--save-df', action='store_true',
                        help="Flag to save propka results dataframe")
    parser.add_argument('--report-propka', action='store_true',
                        help="Flag to report propka results")

    args = parser.parse_args()
    in_path = os.path.abspath(args.inpath)
    out_path = os.path.abspath(args.outpath)
    os.makedirs(out_path, exist_ok=True)
    if args.prefix is None:
        # Get the file stem from infile path
        if sys.platform == "win32":
        # This applies to both 32- and 64-bit Windows machines
            out_prefix = in_path.split('\\')[-1].rstrip('.pdb') 
        else:
            out_prefix = in_path.split('/')[-1].rstrip('.pdb')
    else:
        out_prefix = args.prefix

    out_pdb = os.path.join(out_path, out_prefix+'_protonated.pdb')
    biomol = charmm_biomol(in_path)
    biomol.auto_protonate(ph = args.ph)
    biomol.save_pdb(out_pdb)
    if args.save_df:
        biomol.pka_df.to_csv(out_path+'.csv')
    if args.report_propka:
        print(biomol.pka_str)
