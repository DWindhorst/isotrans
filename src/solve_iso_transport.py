'''
Created on 15.05.2023

@author: poudel-b
'''
# -*- coding: utf-8 -*-

from math import exp
from scipy.sparse import lil_matrix
import numpy as np
from scipy.sparse.linalg import spsolve

# import os

# TODO: diffusion time scale

# Function to run 1D vertical transport


def run_1D_model(Ci, soil, layer, time, Q): # n_layers, dz, n_timesteps, dt, ql, qv):

    theta, theta_0, theta_sat, porosity, tortuosity = soil
    n_layers, dz = layer
    n_timesteps, dt = time
    ql, qv = Q

    Ci_layers = [Ci] * n_layers

    Ci_t = []
    Ci_t.append(Ci_layers)

    for t in range(1, n_timesteps):  # TODO review n_time steps once again

        A1_ij = []
        A2_ij = []
        A3_ij = []
        Bij = []

        for layer in range(n_layers):

            def dl_i(tortuosity, theta, T=283.15, Isotopologue='2H', formulation="Melayah", ignore_dl_i=False):

                """
               Calculates the liquid diffusivity for isotope i in the soil liquide phase (m**2/s)
                """

                try:
                    if formulation == "Melayah":
                        if "18O" == Isotopologue:
                            a_i = 0.9669
                        elif "2H" == Isotopologue:
                            a_i = 0.9833
                        else:
                            raise NotImplementedError

                        dl_i = a_i * 1.10 ** -9 * exp(-(535400 / T ** 2) + (1393.3 / T) + 2.1876)

                    # SLI: cable_sli_solve.f90::L2381-2390
                    elif formulation == "Cuntz":
                        if "18O" == Isotopologue:
                            a_i = 1 / 1.026
                        elif "2H" == Isotopologue:
                            a_i = 1 / 1.013
                        else:
                            raise NotImplementedError

                        dl_i = a_i * 1e-7 * exp(-577 / (T - 145))

                    else:
                        raise NotImplementedError

                    return dl_i * theta * tortuosity

                except ValueError as err:
                    print(err)
                    raise NotImplementedError

            def d_l_eff(ql, tortuosity, theta, dispersivity=0):

                """
                Calculates the total (effective) liquid diffusivity for isotope i in the soil liquide phase (m**2/s)

                """
                d_l = dl_i(tortuosity, theta)
                d_L_eff = d_l + dispersivity * abs(ql)

                return d_L_eff

            def dv_i_eff(theta, theta_0, porosity, tortuosity, Isotopologue='2H'):
                """
                Calculates the total (effective) vapour diffusivity of isotope in soil air space (m**2/s)
                """
                try:
                    nD = 0.67 + 0.33 * exp(1 - theta / theta_0)

                except ZeroDivisionError as err:
                    nD = 0.67

                try:
                    dv_eff = (porosity - theta) * tortuosity * dv_air() \
                             * (dv_i(Isotopologue) / dv_air()) ** nD

                    return dv_eff

                except ValueError as err:
                    print(err)
                    raise NotImplementedError

            def dv_i(Isotopologue='2H', ignore_dv_i=False):
                """
                Calculates the vapour diffusivity of Isotopologue in free air (m**2/s)
                """
                try:
                    # SLI: cable_sli_solve.f90::L2366 - L2367
                    if "18O" == Isotopologue:
                        b_i = 1.0285  # H218O diffusivity in air (Merlivat 1978) // SLI: 1/b_i = alphak_vdiff
                    elif "2H" == Isotopologue:
                        b_i = 1.0251  # HDO diffusivity in air (Merlivat 1978) // SLI: 1/b_i = alphak_vdiff
                    else:
                        raise NotImplementedError

                    dv = dv_air()

                    if ignore_dv_i == True:
                        return dv
                    else:
                        # SLI: cable_sli_solve.f90::L2373 - L2374
                        dv_i = dv / b_i  # SLI: dv = Dvs ; dv_i = Divs ; 1/b_i = alphak_vdiff

                        return dv_i

                except ValueError as err:
                    print(err)
                    raise NotImplementedError

            def dv_air(T=283.15, Pa=10 ** 5):  # TODO: Check the formulation and Units also highlight TODO
                """
                Calculates the vapour diffusivity of water in free air (m**2/s)
                """
                try:
                    dva = 2.17e-5  # vapour diffusivity of water in air at 0 degC [m2/s]
                    dv = dva * 1e5 / Pa * (T / 273.16) ** 1.88

                    return dv

                except ValueError as err:
                    print(err)
                    raise NotImplementedError

            def D_lv_eff(liquid_diffusivity, vapor_diffusivity, beta):

                return liquid_diffusivity + vapor_diffusivity * beta

            def beta(T=283.15, solute='2H', ignore_alpha_i=False, density_h2o_vapour=0.059,
                     density_h2o_liquid=1000.00):
                try:
                    # SLI: cable_sli_solve.f90::L2275 - L2305
                    if ignore_alpha_i == False:
                        if "18O" == solute:
                            alpha_i = exp(-(24844 / T ** 2 + (-76.248) / T + 0.052612))
                        elif "2H" == solute:
                            alpha_i = exp(-(1137 / T ** 2 + (-0.4156) / T - 0.0020667))
                        else:
                            raise NotImplementedError
                        return alpha_i
                    else:
                        return 1.0

                # beta_i = alpha_i according to SLI
                except ValueError as err:
                    print(err)
                    raise NotImplementedError

            vapor_diffusivity = dv_i_eff(theta, theta_0, porosity, tortuosity)
            liquid_diffusivity = d_l_eff(ql, tortuosity, theta)

            Dlv_upper = D_lv_eff(liquid_diffusivity, vapor_diffusivity, beta())
            Dlv_current = D_lv_eff(liquid_diffusivity, vapor_diffusivity, beta())
            Dlv_lower = D_lv_eff(liquid_diffusivity, vapor_diffusivity, beta())

            Dlv_up = (Dlv_upper + Dlv_current) / 2
            Dlv_down = (Dlv_current + Dlv_lower) / 2

            Dv_upper = dv_i_eff(theta, theta_0, porosity, tortuosity)
            Dv_current = dv_i_eff(theta, theta_0, porosity, tortuosity)
            Dv_lower = dv_i_eff(theta, theta_0, porosity, tortuosity)

            current_beta_upper = beta()
            current_beta = beta()
            current_beta_lower = beta()

            ql_up = ql
            ql_down = ql
            qv_up = qv
            qv_down = qv

            Qlv_up = ql_up + qv_up * (current_beta_upper + current_beta) / 2 \
                     - (Dv_upper + Dv_current) / 2 * (current_beta - current_beta_upper) / dz

            Qlv_down = ql_down + qv_down * (current_beta_lower + current_beta) / 2 \
                       - (Dv_lower + Dv_current) / 2 * (current_beta_lower - current_beta) / dz

            theta_eff = theta + (porosity - theta) * current_beta

            delta_z = dz
            delta_t = dt

            if layer == 0:  # upper boundary

                A1 = 0
                A1_ij.append(A1)

                A2 = 1
                A2_ij.append(A2)

                A3 = 0
                A3_ij.append(A3)

                B = 1
                Bij.append(B)

            if layer != 0 and layer != n_layers - 1:  # intermediate layers

                A1 = - Dlv_up / dz - Qlv_up / 2
                A1_ij.append(A1)

                A2 = delta_z / delta_t * theta_eff + Dlv_down / dz + Qlv_down / 2 \
                     + Dlv_up / dz - Qlv_up / 2
                A2_ij.append(A2)

                A3 = - Dlv_down / dz + Qlv_down / 2
                A3_ij.append(A3)

                B = delta_z / delta_t * theta_eff * Ci_layers[layer]
                Bij.append(B)

            if layer == n_layers - 1:  # lower boundary

                A1 = 0
                A1_ij.append(A1)

                A2 = 1
                A2_ij.append(A2)

                A3 = 0
                A3_ij.append(A3)

                B = 0
                Bij.append(B)

        A = lil_matrix((n_layers, n_layers))
        A.setdiag(A1_ij[1::], k=-1)     # TODO: Check for the diagonal-1 matrix starts from - A[1,0]
        A.setdiag(A2_ij, k=0)
        A.setdiag(A3_ij, k=1)

        B = np.asarray(Bij)
        A = A.tocsr()
        Ci_layers = spsolve(A, B).tolist()
        Ci_t.append(Ci_layers)

    C = np.array(Ci_t).transpose()
    return C




