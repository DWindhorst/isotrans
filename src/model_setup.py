
'''
Created on 15.05.2023
@author: poudel-b
'''

# -*- coding: utf-8 -*-

# import tkinter
# import matplotlib
# matplotlib.use('TkAgg')

import matplotlib.pyplot as plt
import numpy as np
import solve_iso_transport as iso_solve


def main():

    # TODO: diffusion time scale

    # initialize
    # units in kg - m - days

    Ci = 0  # initial isotope concentration

    # layer information
    soil_depth = 20  # m
    layers = 100  # No of discrete layers
    dz = soil_depth / layers  # m

    # Time information
    time_final = 200  # days
    dt = 1/86400 * 60 * 60  # 1 hour
    time_steps = int(time_final / dt)

    # Flux information
    ql = [0.01]  # m / day  [ql / 86400  per sec]
    qv = [0.0]  # m/day


    # soil properties
    theta = 0.5  # initial volumetric water content
    theta_0 = 0.01  # volumetric (residual) water content at high suctions m3/m3
    theta_sat = 0.35  # saturate water content m3/m3
    porosity = 0.35  # soil porosity m3/m3
    tortuosity = 0 #0.67  # tortuosity of the soil m / m   #TODO: tortuosity = 0 for diffusion = 0

    soil = [theta, theta_0, theta_sat, porosity, tortuosity]
    layer = [layers, dz]
    time = [time_steps, dt]
    Q = zip(ql, qv)

    # Run multiple simulations and store the results
    C_t = []
    for q_l, q_v in Q:

        C_t.append(iso_solve.run_1D_model(Ci, soil, layer, time, [q_l, q_v]))

    #print(C_t)


    # VISUALIZATION

    # Concentration vs time
    fig1, ax = plt.subplots(figsize=(10, 5), dpi=300)
    for C, q_l, q_v in zip(C_t, ql, qv):
        for l in range(0, layers, 10):

            ax.plot(np.arange(0, time_steps) * dt, C[l], label = '{}m'.format(l * dz))
                     #label='ql = {}, qv = {}, dl = {}, dv = {},'.format(q_l, qv_down, 0, 0))

        ax.set_xlabel('time [days]')
        ax.set_ylabel('conc')
        ax.set_title('conc at breakthrough (ql = {}, qv = {})'.format(q_l, q_v))
        #plt.yscale('log')
        ax.legend()

    #timestamp1 = datetime.now().strftime("%Y%m%d-%H%M%S")
    #fig1.savefig( 'D:\\Isotope transport\\Scripts\\output\\{}.png'.format(timestamp1))
    plt.show()

    # Concentration vs depth

    fig2, ax1 = plt.subplots(figsize=(10, 5), dpi=300)
    for C, q_l, q_v in zip(C_t, ql, qv):
        for t in range(0, time_steps, 48):


            ax1.plot(C[:, t], -np.arange(0, layers) * dz,
                     label='day {}'.format(t / 24))

        ax1.set_xlabel('conc')
        ax1.set_ylabel('depth [m]')
        ax1.set_title('conc profiles (ql = {}, qv = {}) '.format(q_l, q_v))
        #plt.yscale('log')
        ax1.legend()

    #timestamp2 = datetime.now().strftime("%Y%m%d-%H%M%S")
    #fig2.savefig( 'D:\\Isotope transport\\Scripts\\output\\{}.png'.format(timestamp2))
    plt.show()


if __name__ == '__main__':
    main()

