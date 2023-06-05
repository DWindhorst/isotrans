""""
"""
# -*- coding: utf-8 -*-

from iso_Storages import iso_storage as Storage
import iso_flux_connections as flx
import matplotlib.pyplot as plt

# Define a new Storages

S1 = Storage(1, 2, 1)
S2 = Storage(1,2,1,conc_iso_liquid={"2H": 0.0, "18O": 0.0})

# Set a connection between storages
c = flx.linear_connection(S1,S2,0.2,0.05,'2H', 100)
result = c.run()

# Visualize
fig = plt.figure(figsize=(10, 6 ), dpi=100)
ax = fig.subplots()

lgd = ['Liquid', 'vapor']
i = 0
for conc in result:

    ax.plot(conc[0], label='Storage 1_' + lgd[i])
    ax.plot(conc[1], label='Storage 2_' + lgd[i])

    ax.set_xlabel('time [days]')
    ax.set_ylabel('conc [kg/m3] ')
    ax.legend()
    i +=1
#ax.legend(('Liquid', 'Vapor'))
plt.show()



#S1.get_storage_i('2H')
#print(S1)


#C = flx.linear_connection(S1, S2, 10, 0, '2H')

""""

adl = flx.liquid_advection(q_l=0.5)
adv = flx.vapor_advection(q_v=0.3)
dil = flx.liquid_diffusion()
div = flx.vapor_diffusion()

"""
print('done YaY!')
