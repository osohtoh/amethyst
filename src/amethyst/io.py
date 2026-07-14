import os
import re
from dataclasses import dataclass
from typing import Optional, Union

from loguru import logger
from rdkit.Chem.rdmolfiles import MolFromSmiles, MolToSmiles
from rdkit.Chem.rdchem import Mol

from amethyst.utils import mols_to_str


@dataclass
class Substituents:
    """Dataclass used for storing a listof R-groups with their respective R#.

    Attributes:
        r_num (int): R#
        subs (list[Mol]): listof substituents in Mol object

    """

    r_num: int
    subs: list[Mol]


def parse_file_input(
    filepath: str,
    r_num: Optional[int] = None,
    delimiter: Optional[str] = None,
    multiple_rs: Optional[bool] = False,
) -> list[Substituents]:
    """Parses provided file to a Substituents dataclass. Accepts R-groups marked as either isotope labels or atom maps. Newline separated file can only be for one R#.

    Args:
        filepath (str): Path to file with R-groups. Required.
        r_num (int): Number of R-group attached to the Substituents dataclass. Can be passed as R# (case insensitive) in a filename.
        delimiter (Optional[str], optional): Delimiter for one-line files. Defaults to newline.
        multiple_rs (Optional[bool], optional): Flag determining if each line in a file should be considered different R-group. First line is R1 and so on. Defaults to False.

    Raises:
        FileNotFoundError: Raised if supplied path isn't a file.
        ValueError: Raised if no r_num was provided.

    Returns:
        Substituents: Parsed file saved to dataclass.
    """
    if os.path.isfile(filepath):
        logger.debug("File path is good.")
        pass
    else:
        logger.error(f"{filepath} is not a file!")
        raise FileNotFoundError(f"{filepath} is not a file!")

    if r_num is None:
        if multiple_rs:
            r_num = 1
        else:
            path_split = re.split(r"(\\\\)|(/)|(\\)", filepath)
            m = re.match("[rR][0-9]+", path_split[-1])
            if m is not None:
                r_num = int(m.group(0))
            else:
                raise ValueError("R# is missing.")
    logger.debug(f"File given for R{r_num}.")

    subs_list: list[Mol] = []
    r = f"[*:{r_num}]"

    subs: list[Substituents] = []

    r_group_regex = r"\[[0-9]+\*\]|\[[^\]]+\:[0-9]+\]"

    with open(filepath, "r") as file:
        if delimiter is None:
            for i in file:
                m = re.sub(r_group_regex, r, i)
                mol = MolFromSmiles(m)
                if mol is None:
                    raise ValueError(f"Could not parse SMILES: {m.strip()!r}")
                subs_list.append(mol)
                logger.debug(f"SMILES added: {m}")
            if not subs_list:
                raise ValueError(f"No substituents found in {filepath}")
            subs.append(Substituents(r_num, subs_list))
        else:
            lines: list[str] = file.readlines()
            logger.debug(lines)
            if multiple_rs:
                # TODO - Write tests
                for line in lines:
                    r = f"[*:{r_num}]"
                    m = re.sub(r_group_regex, r, line)
                    split_lines: list[str] = m.split(delimiter)
                    mols = []
                    for smi in split_lines:
                        mol = MolFromSmiles(smi)
                        if mol is None:
                            raise ValueError(f"Could not parse SMILES: {smi.strip()!r}")
                        mols.append(mol)
                    if not mols:
                        raise ValueError(f"No substituents found for R{r_num}")
                    subs.append(Substituents(r_num, mols))
                    logger.debug(f"R{r_num} SMILES: {split_lines}")
                    r_num = r_num + 1
            else:
                for line in lines:
                    m = re.sub(r_group_regex, r, line)
                    split_lines: list[str] = m.split(delimiter)
                    mols = []
                    for smi in split_lines:
                        mol = MolFromSmiles(smi)
                        if mol is None:
                            raise ValueError(f"Could not parse SMILES: {smi.strip()!r}")
                        mols.append(mol)
                    if not mols:
                        raise ValueError(f"No substituents found in {filepath}")
                    subs.append(Substituents(r_num, mols))
                    logger.debug(f"SMILES added: {split_lines}")

    logger.debug(f"Final sub list: {mols_to_str(subs_list)}")

    return subs


def parse_mol_input(mols: list[list[Union[Mol, str]]]) -> list[Substituents]:
    """Parses listof molecules to a Substituents class. R# is handled via the listindex (n+1) e.g., first listof Mol's in the listpassed will have R1 number and so on.

    Args:
        mols (list[list[Mol]]): listcontaining another listof R-groups.

    Raises:
        ValueError: Raised when input isn't Mol or str.

    Returns:
        list[Substituents]: Returns subs parsed into a listof Substituents dataclass.
    """
    r_group_regex = r"\[[0-9]+\*\]|\[[^\]]+\:[0-9]+\]"

    r_num = 1
    substituents_list = []
    for i in mols:
        if not i:
            raise ValueError(f"No substituents provided for R{r_num}")
        if type(i[0]) is Mol:
            logger.debug("Mol input detected.")
            mols_smi = [MolToSmiles(x) for x in i]
        elif isinstance(i[0], str):
            logger.debug("String input detected.")
            mols_smi = mols[r_num - 1]
        else:
            raise ValueError("Wrong input type.")

        smis_relabelled = []
        for j in mols_smi:
            r = f"[*:{r_num}]"
            logger.debug(f"Inner: {j}")
            if isinstance(j, str):
                m = re.sub(r_group_regex, r, j)
            elif type(j) is Mol:
                smi = MolToSmiles(j)
                m = re.sub(r_group_regex, r, smi)
            else:
                raise ValueError("Wrong input type.")
            logger.debug(f"Relabelled SMILES: {m}")
            mol = MolFromSmiles(m)
            if mol is None:
                raise ValueError(f"Could not parse SMILES: {m!r}")
            smis_relabelled.append(m)
            logger.debug(f"SMILES relabelled: {smis_relabelled}")

        mols_relabelled = [MolFromSmiles(x) for x in smis_relabelled]
        substituents_list.append(Substituents(r_num, mols_relabelled))
        logger.debug(f"R{r_num} SMILES: {mols_to_str(mols_relabelled)}")

        r_num = r_num + 1

    return substituents_list
