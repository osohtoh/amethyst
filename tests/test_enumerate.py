import pytest

from amethyst.amethyst import enumerate
from amethyst.io import parse_file_input
from rdkit.Chem.AllChem import CanonSmiles, MolFromSmiles

path = "tests\\newline_input.txt"
core = ""


def test_enumerate_file():
    # enumerate(core, subs_path=path)
    pass


def test_enumerate_mol():
    core = MolFromSmiles("[1*]c1cnccc1[2*]")
    subs = [["[1*]CC", "[1*]CCC"], ["[2*]NC", "[2*]NCC"]]
    out = [
        CanonSmiles(x)
        for x in ["CCc1cnccc1NC", "CCNc1ccncc1CC", "CCCc1cnccc1NC", "CCCc1cnccc1NCC"]
    ]
    generated = enumerate(core=core, multiple_rs=True, subs_mol=subs, output_smi=True)
    assert len(generated) == len(out)
    assert generated == out
