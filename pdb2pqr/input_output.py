"""Functions related to reading and writing data."""
import logging
import io
from collections import Counter
from pathlib import Path
from sys import path as sys_path
import requests
from . import psize
from . import inputgen
from . import cif
from . import pdb
from . import definitions as defns
from .config import FORCE_FIELDS, TITLE_FORMAT_STRING, VERSION
from .config import FILTER_WARNINGS_LIMIT, FILTER_WARNINGS
from .config import AA_DEF_PATH, NA_DEF_PATH, PATCH_DEF_PATH


_LOGGER = logging.getLogger(__name__)


class DuplicateFilter(logging.Filter):
    """Filter duplicate messages."""
    def __init__(self):
        super().__init__()
        self.warn_count = Counter()

    def filter(self, record):
        """Filter current record."""
        if record.levelname == "WARNING":
            for fwarn in FILTER_WARNINGS:
                if record.getMessage().startswith(fwarn):
                    self.warn_count.update([fwarn])
                    if self.warn_count[fwarn] > FILTER_WARNINGS_LIMIT:
                        return False
                    elif self.warn_count[fwarn] == FILTER_WARNINGS_LIMIT:
                        _LOGGER.warning(
                            "Suppressing further '%s' messages", fwarn)
                        return False
                    else:
                        return True
        return True


def print_protein_atoms(atomlist, chainflag=False, pdbfile=False):
    """Get text lines for specified atoms

    Args:
        atomlist:  The list of atoms to include (list)
        chainflag:  Flag whether to print chainid or not
    Returns:
        text:  list of (stringed) atoms (list)
    """
    text = []
    currentchain_id = None
    for atom in atomlist:
        # Print the "TER" records between chains
        if currentchain_id is None:
            currentchain_id = atom.chain_id
        elif atom.chain_id != currentchain_id:
            currentchain_id = atom.chain_id
            text.append("TER\n")

        if pdbfile is True:
            text.append("%s\n" % atom.get_pdb_string())
        else:
            text.append("%s\n" % atom.get_pqr_string(chainflag=chainflag))
    text.append("TER\nEND")
    return text


def get_old_header(pdblist):
    """Get old header from list of PDBs.

    Args:
        pdblist:  list of PDBs
    Returns:
        Old header as string.
    """
    old_header = io.StringIO()
    header_types = (
        pdb.HEADER, pdb.TITLE, pdb.COMPND, pdb.SOURCE, pdb.KEYWDS, pdb.EXPDTA,
        pdb.AUTHOR, pdb.REVDAT, pdb.JRNL, pdb.REMARK, pdb.SPRSDE, pdb.NUMMDL)
    for pdb_obj in pdblist:
        if not isinstance(pdb_obj, header_types):
            break
        old_header.write(str(pdb_obj))
        old_header.write('\n')
    return old_header.getvalue()


def print_pqr_header(
        pdblist, atomlist, reslist, charge, force_field, ph_calc_method, ph,
        ffout, include_old_header=False):
    """Print the header for the PQR file

    Args:
        pdblist:  list of lines from original PDB with header
        atomlist: A list of atoms that were unable to have charges assigned
                  (list)
        reslist:  A list of residues with non-integral charges (list)
        charge:  The total charge on the protein (float)
        ff:  The forcefield name (string)
        ph:  pH value, if any. (float)
        ffout:  ff used for naming scheme (string)
    Returns
        header:  The header for the PQR file (string)
    """
    if force_field is None:
        force_field = 'User force field'
    else:
        force_field = force_field.upper()
    head = "\n".join(
        [
            "REMARK   1 PQR file generated by PDB2PQR",
            "REMARK   1 %s" % TITLE_FORMAT_STRING.format(version=VERSION),
            "REMARK   1",
            "REMARK   1 Forcefield Used: %s" % force_field, ""])
    if ffout is not None:
        head += "REMARK   1 Naming Scheme Used: %s\n" % ffout
    head += head + "REMARK   1\n"
    if ph_calc_method is not None:
        head += "\n".join(
            [
                (
                    "REMARK   1 pKas calculated by %s and assigned "
                    "using pH %.2f\n"
                ) % (ph_calc_method, ph), "REMARK   1", ""])
    if len(atomlist) != 0:
        head += "\n".join(
            [
                "REMARK   5 WARNING: PDB2PQR was unable to assign charges",
                "REMARK   5 to the following atoms (omitted below):", ""])
        for atom in atomlist:
            head += "REMARK   5    %i %s in %s %i\n" % (
                atom.serial, atom.name, atom.residue.name,
                atom.residue.res_seq)
        head += "\n".join(
            [
                "REMARK   5 This is usually due to the fact that this residue "
                "is not",
                "REMARK   5 an amino acid or nucleic acid; or, there are no "
                "parameters",
                "REMARK   5 available for the specific protonation state of "
                "this",
                "REMARK   5 residue in the selected forcefield.",
                "REMARK   5", ""])
    if len(reslist) != 0:
        head += "\n".join([
            "REMARK   5 WARNING: Non-integral net charges were found in",
            "REMARK   5 the following residues:", ""])
        for residue in reslist:
            head += (
                "REMARK   5    %s - Residue Charge: %.4f\n"
                % (residue, residue.charge))
        head += "REMARK   5\n"
    head += "\n".join([
        "REMARK   6 Total charge on this protein: %.4f e\n" % charge,
        "REMARK   6", ""])
    if include_old_header:
        head += "\n".join(["REMARK   7 Original PDB header follows",
                           "REMARK   7", ""])
        head += get_old_header(pdblist)
    return head


def print_pqr_header_cif(atomlist, reslist, charge, force_field,
                         ph_calc_method, ph, ffout, include_old_header=False):
    """Print the header for the PQR file in cif format.

    Args:
        atomlist: A list of atoms that were unable to have charges assigned
                  (list)
        reslist:  A list of residues with non-integral charges (list)
        charge:  The total charge on the protein (float)
        force_field:  The forcefield name (string)
        ph:  pH value, if any. (float)
        ffout:  ff used for naming scheme (string)
    Returns
        header:  The header for the PQR file (string)
    """
    if force_field is None:
        force_field = "User force field"
    else:
        force_field = force_field.upper()

    header = "\n".join([
        "#", "loop_", "_pdbx_database_remark.id",
        "_pdbx_database_remark.text", "1", ";",
        "PQR file generated by PDB2PQR",
        TITLE_FORMAT_STRING.format(version=VERSION), "",
        "Forcefield used: %s\n" % force_field, ""])
    if ffout is not None:
        header += "Naming scheme used: %s\n" % ffout
    header += "\n"
    if ph_calc_method is not None:
        header += (
            "pKas calculated by %s and assigned using pH %.2f\n"
            % (ph_calc_method, ph))
    header += "\n".join([";", "2", ";", ""])
    if len(atomlist) > 0:
        header += "\n".join(["Warning: PDB2PQR was unable to assign charges",
                             "to the following atoms (omitted below):", ""])
        for atom in atomlist:
            header += "    %i %s in %s %i\n" % (
                atom.serial, atom.name, atom.residue.name,
                atom.residue.res_seq)
        header += "\n".join([
            "This is usually due to the fact that this residue is not",
            "an amino acid or nucleic acid; or, there are no parameters",
            "available for the specific protonation state of this",
            "residue in the selected forcefield.", ""])
    if len(reslist) > 0:
        header += "\n".join([
            "Warning: Non-integral net charges were found in",
            "the following residues:", ""])
        for residue in reslist:
            header += "    %s - Residue Charge: %.4f\n" % (
                residue, residue.charge)
    header += "\n".join([
        ";", "3", ";",
        "Total charge on this protein: %.4f e" % charge, ";", ""])
    if include_old_header:
        _LOGGER.warning("Including original CIF header not implemented.")
    header += "\n".join([
        "#", "loop_", "_atom_site.group_PDB", "_atom_site.id",
        "_atom_site.label_atom_id", "_atom_site.label_comp_id",
        "_atom_site.label_seq_id", "_atom_site.Cartn_x",
        "_atom_site.Cartn_y", "_atom_site.Cartn_z",
        "_atom_site.pqr_partial_charge", "_atom_site.pqr_radius", ""])
    return header


def dump_apbs(output_pqr, output_path):
    """Generate and dump APBS input files related to output_pqr.

    Args:
        output_pqr:  path to PQR file used to generate APBS input file
        output_path:  path for APBS input file output
    """
    method = "mg-auto"
    size = psize.Psize()
    size.parse_input(output_pqr)
    size.run_psize(output_pqr)
    input_ = inputgen.Input(output_pqr, size, method, 0, potdx=True)
    input_.print_input_files(output_path)


def test_for_file(name, type_):
    """Test for the existence of a file with a few name permutations.

    Args:
        name:  name of file
        type_:  type of file
    Returns:
        path to file or None
    """
    if name is None:
        return ''
    test_names = [name, name.upper(), name.lower()]
    test_suffixes = ["", ".%s" % type_.upper(), ".%s" % type_.lower()]
    test_dirs = [
        Path(p).joinpath("pdb2pqr", "dat") for p in sys_path + [Path.cwd()]]
    if name.lower() in FORCE_FIELDS:
        name = name.upper()
    for test_dir in test_dirs:
        test_dir = Path(test_dir)
        for test_name in test_names:
            for test_suffix in test_suffixes:
                test_path = test_dir / (test_name + test_suffix)
                if test_path.is_file():
                    _LOGGER.debug("Found %s file %s", type_, test_path)
                    return test_path
    _LOGGER.warning("Unable to find %s file for %s", type_, name)
    return ""


def test_names_file(name):
    """Test for the *.names file that contains the XML mapping.

    Args:
        name:  The name of the forcefield (string)
    Returns
        path:  The path to the file (string)
    """
    return test_for_file(name, "NAMES")


def test_dat_file(name):
    """Test for the existence of the forcefield file with a few name
    permutations.

    Args:
        name of forcefield
    Returns:
        filename or empty string
    """
    return test_for_file(name, "DAT")


def test_xml_file(name):
    """Test for the existence of the forcefield file with a few name
    permutations.

    Args:
        name of the xml file
    Returns:
        filename or empty string
    """
    return test_for_file(name, "xml")


def get_pdb_file(name):
    """Obtain a PDB file.

    First check the path given on the command line - if that file is not
    available, obtain the file from the PDB webserver at
    http://www.rcsb.org/pdb/

    TODO - this should be a context manager (to close the open file)

    Args:
        name:  Name of PDB file to obtain (string)

    Returns
        file:  File object containing PDB file (file object)
    """
    path = Path(name)
    if path.is_file():
        return open(path, "rt", encoding="utf-8")
    else:
        url_path = "https://files.rcsb.org/download/" + path.stem + ".pdb"
        _LOGGER.debug("Attempting to fetch PDB from %s", url_path)
        resp = requests.get(url_path)
        if resp.status_code != 200:
            errstr = "Got code %d while retrieving %s" % (
                resp.status_code, url_path)
            raise IOError(errstr)
        return io.StringIO(resp.text)


def get_molecule(input_path):
    """Get molecular structure information.

    Args:
        input_path:  structure file PDB ID or path
    Returns:
        list of molecule records (lines)
        Boolean indicating whether entry is CIF
    Raises:
        RuntimeError:  problems with structure file
    """
    path = Path(input_path)
    input_file = get_pdb_file(input_path)
    is_cif = False

    if path.suffix.lower() == ".cif":
        pdblist, errlist = cif.read_cif(input_file)
        is_cif = True
    else:
        pdblist, errlist = pdb.read_pdb(input_file)

    if len(pdblist) == 0 and len(errlist) == 0:
        raise RuntimeError("Unable to find file %s!" % path)

    if len(errlist) != 0:
        if is_cif:
            _LOGGER.warning("Warning: %s is a non-standard CIF file.\n", path)
        else:
            _LOGGER.warning("Warning: %s is a non-standard PDB file.\n", path)
        _LOGGER.error(errlist)

    return pdblist, is_cif


def get_definitions(aa_path=AA_DEF_PATH, na_path=NA_DEF_PATH,
                    patch_path=PATCH_DEF_PATH):
    """Get topology definition files.

    Args:
        aa_path:  likely location of amino acid topology definitions
        na_path:  likely location of nucleic acid topology definitions
        patch_path:  likely location of patch topology definitions

    Returns:
        Definitions object.
    """
    aa_path_ = test_xml_file(aa_path)
    if not aa_path_:
        raise FileNotFoundError("Unable to locate %s" % aa_path)
    na_path_ = test_xml_file(na_path)
    if not na_path_:
        raise FileNotFoundError("Unable to locate %s" % na_path)
    patch_path_ = test_xml_file(patch_path)
    if not patch_path_:
        raise FileNotFoundError("Unable to locate %s" % patch_path)
    with open(aa_path_, "rt") as aa_file:
        with open(na_path_, "rt") as na_file:
            with open(patch_path_, "rt") as patch_file:
                definitions = defns.Definition(
                    aa_file=aa_file, na_file=na_file, patch_file=patch_file)
    return definitions

def setup_logger(output_pqr, level='DEBUG'):
    """ Setup the logger to output the log file to the same dir as pqr output"""
    # Get the output logging location
    output_pth = Path(output_pqr)
    log_file = Path(output_pth.parent, output_pth.stem + '.log')
    _LOGGER.info('Logs stored: %s', log_file)

    logging.basicConfig(
        filename=log_file,
        level=getattr(logging, level)
        )
