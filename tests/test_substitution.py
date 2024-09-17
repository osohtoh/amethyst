import pytest
from rdkit.Chem.AllChem import CanonSmiles, MolFromSmiles, MolToSmiles

from amethyst.io import parse_mol_input
from amethyst.substitution import general_sub

core = MolFromSmiles("[*:2]Cc1cccc([*:1])c1")
r = [["[*:1]NC", "[*:1]NCc1ccccc1"], ["[*:2]c1ccccc1", "[*:2]OC"]]


def test_general_sub_relabelled():
    subs = parse_mol_input(r)
    result = general_sub(core, subs)
    assert len(result) == 4
    assert [CanonSmiles(MolToSmiles(x)) for x in result] == [
        CanonSmiles(x)
        for x in [
            "CNc1cccc(Cc2ccccc2)c1",
            "CNc1cccc(COC)c1",
            "c1ccc(CNc2cccc(Cc3ccccc3)c2)cc1",
            "COCc1cccc(NCc2ccccc2)c1",
        ]
    ]


def test_general_sub_unlabelled():
    subs = parse_mol_input(r)
    core_unlabelled = MolFromSmiles("[2*]Cc1cccc([1*])c1")
    result = general_sub(core_unlabelled, subs, False)
    assert len(result) == 4
    assert [CanonSmiles(MolToSmiles(x)) for x in result] == [
        CanonSmiles(x)
        for x in [
            "CNc1cccc(Cc2ccccc2)c1",
            "CNc1cccc(COC)c1",
            "c1ccc(CNc2cccc(Cc3ccccc3)c2)cc1",
            "COCc1cccc(NCc2ccccc2)c1",
        ]
    ]
