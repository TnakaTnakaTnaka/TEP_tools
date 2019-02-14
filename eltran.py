#!/usr/bin/env python
#
# tran.py
#
# script to calculate the transmittance by the Brillouin
# zone integration of k-resolved transmittance.
#
#
# Copyright (c) 2018 Yuto Tanaka
#


"""
--- how to use ---

$ python tran.py


"""


import os
import sys
import fnmatch
import numpy as np


def data_column():
    """find k-resolved transmittance file"""
    flag = 0
    for file_name in os.listdir(os.getcwd()):
        if fnmatch.fnmatch(file_name, '*.tran0_0'):
            prefix = file_name.split(".")[0]
            data = np.loadtxt(file_name)
            energy = data.T[3]
            num_data = int(np.shape(data)[0])
            num_column = int(np.shape(data)[1])
            flag = 1

    # If the *.tran0_0 file is not exist, exit this program.
    if flag == 0:
        print("*.tran0_0 file is not exist.")
        sys.exit()

    return prefix, num_data, num_column, energy


def transmission(num_data, num_column):
    """integrate transmission in 1st BZ"""
    tran = np.zeros(num_data)
    kpoint = 0
    for file_name in os.listdir('.'):
        if fnmatch.fnmatch(file_name, '*.tran*_*'):
            kpoint += 1  # kpoint countor
            data = np.loadtxt(file_name)
            if num_column == 7:  # SOC
                tran += data.T[5]

            elif num_column == 9:  # Without SOC
                tran += data.T[5] + data.T[7]

            else:
                print("The number of column is not proper.")
                exit()

    tran /= kpoint
    print("Loaded", kpoint, "file(s).")

    return tran


def main():
    """Set transmission array"""
    prefix, num_data, num_column, energy = data_column()
    tran = np.zeros([num_data, 2])

    tran.T[0] = energy
    tran.T[1] = transmission(num_data, num_column)

    # Save transmission file.
    tran_file = prefix + ".tran"
    np.savetxt(tran_file, tran, delimiter='  ')


if __name__ == "__main__":
    main()
