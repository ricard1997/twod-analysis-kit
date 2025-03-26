"""
Location of data files
======================

Use as ::

    from twodanalysis.data.files import *

"""

__all__ = [
    "MDANALYSIS_LOGO",  # example file of MDAnalysis logo
]

import importlib.resources
import pooch
data_directory = importlib.resources.files("twodanalysis") / "data"




MDANALYSIS_LOGO = data_directory / "mda.txt"



#test

MEMBRANE_XTC = pooch.retrieve(
    url="https://zenodo.org/records/14834046/files/md_membrane_nowater.xtc",
    known_hash=None,
)
MEMBRANE_TPR = pooch.retrieve(
    url="https://zenodo.org/records/14834046/files/md_membrane_nowater.tpr",
        known_hash=None,
)


MD_NOWATER_TPR = pooch.retrieve(
    url="https://zenodo.org/records/14834046/files/md_biopolymer_nowater.tpr",
        known_hash=None,
)
MD_NOWATER_PDB = pooch.retrieve(
    url="https://zenodo.org/records/14834046/files/md_biopolymer_nowater.pdb",
        known_hash=None,
)
MD_TRAJ = pooch.retrieve(
    url="https://zenodo.org/records/14834046/files/md_biopolymer_nowater.xtc",
        known_hash=None,
)


