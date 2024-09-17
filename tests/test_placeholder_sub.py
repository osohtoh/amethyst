from typing import List

import pytest
from rdkit.Chem.AllChem import CanonSmiles, Mol, MolFromSmiles, MolToSmiles

from amethyst.substitution import placeholder_atom_sub

# from loguru import logger
# logger.add(sink='substitution.log', level=10)


core_smi_N: str = "CCN"  # Ethylamine
placeholder_atom_N: str = "N"  # Nitrogen atom
r_groups_N: List[str] = ["NCC", "NC(O)C"]
core_mol_N: Mol = MolFromSmiles(core_smi_N)

core_smi_Lu: str = "CC[Lu]"
placeholder_atom_Lu: str = "[Lu]"
r_groups_Lu: List[str] = ["NCC", "NC(O)C"]
core_mol_Lu: Mol = MolFromSmiles(core_smi_Lu)


def test_placeholder_substitution_front():
    sub_mol_N = placeholder_atom_sub(core_mol_N, placeholder_atom_N, r_groups_N)
    assert MolToSmiles(sub_mol_N[0]) == CanonSmiles("CCNCC")
    assert MolToSmiles(sub_mol_N[1]) == CanonSmiles("CCNC(O)C")
    sub_mol_Lu = placeholder_atom_sub(core_mol_Lu, placeholder_atom_Lu, r_groups_Lu)
    assert MolToSmiles(sub_mol_Lu[0]) == CanonSmiles("CCNCC")
    assert MolToSmiles(sub_mol_Lu[1]) == CanonSmiles("CCNC(O)C")


core_smi_inner_N: str = "CNC"
core_mol_inner_N: Mol = MolFromSmiles(core_smi_inner_N)
core_smi_inner_Lu: str = "C[Lu]C"
core_mol_inner_Lu: Mol = MolFromSmiles(core_smi_inner_Lu)

# FIXME - Rewrite proper SMILES, these are wrong (should the insertion be in brackets)

# def test_placeholder_substitution_inner():
#     sub_mol_inner_N = placeholder_atom_sub(core_mol_inner_N, placeholder_atom_N, r_groups_N)
#     assert MolToSmiles(sub_mol_inner_N[0]) == CanonSmiles("CNCCC")
#     assert MolToSmiles(sub_mol_inner_N[1]) == CanonSmiles("CNC(O)CC")
