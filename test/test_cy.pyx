"""
"""
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import cmf

from src import iso_Storages as strg

# Define new isotope storages
S1 = strg.iso_storage(1, 2, 1)
S2 = strg.iso_storage(1,2,1,conc_iso_liquid={"2H": 0.0, "18O": 0.0})

print('done YaY!')

