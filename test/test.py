"""
"""
# -*- coding: utf-8 -*-


#from src import iso_Storages_cy as strg
import example_cy, example_py
import time
# from src import iso_Storages_cy as strg


tic = time.perf_counter()
sum_py = example_py.test(5)
toc = time.perf_counter()

time_py = toc - tic

tic = time.perf_counter()
sum_cy = example_cy.test(5)
toc = time.perf_counter()

time_cy = toc - tic

print('Cython is {} faster then Python'.format(time_cy / time_py))



# Define new isotope stor   ages
#S1 = strg.iso_storage(1, 2, 1)
#S2 = strg.iso_storage(1,2,1,conc_iso_liquid={"2H": 0.0, "18O": 0.0})


print('done YaY!')

