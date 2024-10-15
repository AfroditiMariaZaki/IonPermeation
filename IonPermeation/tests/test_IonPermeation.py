import pytest
import MDAnalysis as mda
from ..ion_permeation import IonPermeation

@pytest.fixture(scope='module')
def setup_universe():
    """
    Create a sample Universe with a test PDB and trajectory file
    """
    topfile = "../data/test.pdb"
    trjfile = "../data/test.xtc"
    u = mda.Universe(topfile, trjfile)

    return u

@pytest.fixture(scope='module')
def initialization(setup_universe: mda.Universe):
    """
    Initialize the IonPermeation module
    """
    u = setup_universe()
    ions = "resname SOD"
    in_sel = "protein and resid 271 272 652 653 654 and name CA"
    out_sel = "protein and resid 690 694 308 312 and name CA"

    # Initialize IonPermeation
    IP = IonPermeation(u, ions, in_sel, out_sel)

    assert IP.u == u
    assert IP.atomgroup == ions
    assert IP.in_sel is not None
    assert IP.out_sel is not None

