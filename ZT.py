#!/usr/bin/env python
#
# ZT.py
#
# script to evaluate thermoelectric properties
# and ZT from transmittance calculated by OpenMX code.
#
#
# Copyright (c) 2019 Yuto Tanaka
#

"""
--- how to use ---

$ python ZT.py tran_file kl_file (--Tmin=100 --Tmax=1000 --dT=100 --cmesh=500 --cpmax=1.0)

"""

import argparse
import numpy as np

usage = "usage: %prog [options]"
parser = argparse.ArgumentParser(usage=usage)
parser.add_argument("tran_file", help="transmittance data file")
parser.add_argument("kp_file", help="thermal conductance data file")
parser.add_argument("--Tmin", action="store", default="100",
                    help="print the minimum temperature you want to calculate.")
parser.add_argument("--Tmax", action="store", default="1000",
                    help="print the maximum temperature you want to calculate.")
parser.add_argument("--dT", action="store", default="100",
                    help="print the width of temperature. Default is 100 K.")
parser.add_argument("--cmesh", action="store", default="500",
                    help="print the number of mesh for chemical potential. Default is 500.")
parser.add_argument("--cpmax", action="store", default="1.0",
                    help="print the maximum of chemical potential. Default is 1.0 eV.")

# Constant data
kb = 8.6173303e-5  # Boltzman constant (ev/K)
G0 = 3.87404194169e-5  # q^2/h  (C^2/J/s)


def fermi_function(x):
    fermi = 1 / (1 + np.exp(x))
    return fermi


def diff_fermi_function(beta, E):
    x = beta * E
    dfermi = beta * np.exp(x) * fermi_function(x) * fermi_function(x)
    return dfermi


def asist_function(E, tran, dfermi):
    num_data = np.shape(E)[0]
    L0 = 0.0
    L1 = 0.0
    L2 = 0.0
    gx = tran * dfermi
    for i in range(num_data-1):
        dE = E[i+1] - E[i]
        L0 += 0.5 * (gx[i] + gx[i+1]) * dE
        L1 += 0.5 * (E[i]*gx[i] + E[i+1]*gx[i+1]) * dE
        L2 += 0.5 * (E[i]*E[i]*gx[i] + E[i+1]*E[i+1]*gx[i+1]) * dE

    return L0, L1, L2


# calculate thermoelectric properties due to electron
def thermoelectric(L0, L1, L2, T):
    G = L0  # electric conductance [G0]
    S = -L1/(L0*T)  # Seebeck coefficient [V/K]
    P = G0 * G * S * S  # power factor [W/K^2]
    K_e = G0 * (L2 - L1 * L1 / L0) / T  # electron thermal conductance [W/K]

    return G, S, P, K_e


def phonon_thermal_concutance(kp_data, T):
    idx_array = np.where(kp_data == T)[0]
    # if the data at the temperature T exits in the kp_data
    if idx_array != []:
        kappa_p = kp_data[idx_array[0]][1]
    else:
        idx1 = np.where(kp_data >= T)[0][0]
        idx2 = np.where(kp_data >= T)[0][1]
        kappa_p = (kp_data[idx2][1] - kp_data[idx1][1]) / \
            (kp_data[idx2][0] - kp_data[idx1][0]) * \
            (T - kp_data[idx1][0]) + kp_data[idx1][1]

    return kappa_p


def calc_thermoelectric_properties(T_min, T_max, T_width, tran_file, kp_file, num_cmesh, max_chem):

    # the number of temperature grid T_point
    if T_width != 0:
        T_point = int((T_max - T_min) / T_width + 1)
    else:
        print("calculate only at %3.1f K." % (T_min))
        T_point = 1

    ZT_data = np.zeros([T_point, 2])  # array od ZT data
    tep_data = np.zeros([num_cmesh, 6])  # array of thermoelectric properties

    # phonon thermal conductance data
    kp_org = np.loadtxt(kp_file)
    kp_data = (kp_org.T[:2]).T

    # transmission data
    data = np.loadtxt(tran_file, delimiter='  ')
    tran = data.T[1]

    # loop temperature
    for t in range(0, T_point):
        ZT_max = 0.0  # initialize maximum of ZT
        T = T_min + T_width * t  # temperature
        if T == 0:
            ZT_data[t][0] = T
            ZT_data[t][1] = ZT_max
            continue

        beta = 1 / (kb*T)

        # find phonon thermal conductance at t K
        kappa_p = phonon_thermal_concutance(kp_data, T)

        # specify the name of the data file
        tep_file = 'tep' + str(int(T)) + '.data'

        # loop chemical potential
        for i in range(num_cmesh):
            # chemical potential mu
            mu = i / float(num_cmesh) * max_chem * 2 - max_chem
            E = data.T[0] - mu

            # differential Ferimi distribution function
            dfermi = diff_fermi_function(beta, E)

            # assist function L_n
            L0, L1, L2 = asist_function(E, tran, dfermi)

            # calculate thermoelectric properties
            ele_cond, seebeck, power, kappa_e = thermoelectric(L0, L1, L2, T)

            # estimate the figure of merit ZT
            ZT = T * power / (kappa_e + kappa_p)

            # compare ZT with the maximum
            if ZT_max < ZT:
                ZT_max = ZT

            # store the thermoelectric properties
            tep_data[i] = [mu, ele_cond, seebeck, power, kappa_e, ZT]

        # store temperature and ZT_max in ZT_data array
        ZT_data[t][0] = T
        ZT_data[t][1] = ZT_max

        # save tep_data
        np.savetxt(tep_file, tep_data, delimiter='  ')

    return ZT_data, ZT_max


def main():

    # Set args
    options = parser.parse_args()
    if options.tran_file:
        tran_file = options.tran_file
    else:
        print("tran file is not selected.")
        exit(1)

    if options.kp_file:
        kp_file = options.kp_file
    else:
        print("thermal conductance file is not selected.")
        exit(1)

    if options.Tmin:
        T_min = float(options.Tmin)
        print("The minimum tempreature : %3.1f K" % (T_min))
        if T_min < 0:
            print("Specified temperature is negative.")
            exit(1)

    if options.Tmax:
        T_max = float(options.Tmax)
        print("The maximum tempreature : %3.1f K" % (T_max))
        if T_min > T_max:
            print("Tmin is larger than Tmax. Check arguments.")
            exit(1)

    if options.dT:
        T_width = float(options.dT)
        if T_max == T_min:
            T_width = 0
        print("The tempreature width   : %3.1f K" % (T_width))

    if options.cmesh:
        num_cmesh = int(options.cmesh)
        if num_cmesh <= 0:
            print("cmaeh should be positive interger.")
            exit(1)

    if options.cpmax:
        max_chem = float(options.cpmax)
        print(
            "The chemical potential range (eV) : [%3.1f, %3.1f]" % (-max_chem, max_chem))

    print("Start calculation...")

    # calculate the thermoelectric properties
    ZT_data, ZT_max = calc_thermoelectric_properties(
        T_min, T_max, T_width, tran_file, kp_file, num_cmesh, max_chem)
    np.savetxt("ZT.data", ZT_data, delimiter='  ')

    print("ZTmax = %5.4f" % (ZT_max))
    print("Finish")


if __name__ == "__main__":
    main()
