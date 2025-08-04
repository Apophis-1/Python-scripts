#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

FC = [[ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. ],
      [ 0., 2.875, -3.966, 2.496, 2.652, -.31, -.511, 0., 0., 0. ],
      [ 0., 1.967, -2.943, 3.755, 1.206, -.727, -.026, 0., 0., 0. ],
      [ 0., 3.862, -2.486, 1.527, 1.247, -.076, -.183, 0., 0., 0. ],
      [ 0., -2.688, -7.843, 2.471, 2.758, -.233, -.747, 0., 0., 0. ],
      [ 0., -2.416, -5.118, 1.869, -.4, .064, .05, 0., 0., 0. ],
      [ 0., 3.017, 2.446, -.306, 0., 0., 0., 0., 0., 0. ],
      [ 0., 2.635, 2.986, -.488, 0., 0., 0., 0., 0., 0. ],
      [ 0., 3.826, 1.619, -.099, 0., 0., 0., 0., 0., 0. ],
      [ 0., -2.204, 1.831, .149, 0., 0., 0., 0., 0., 0. ],
      [ 0., -1.676, 1.075, .404, 0., 0., 0., 0., 0., 0. ],
      [ 0., 4.026, 4.277, -1., 25.48, 36.93, -2.792, -3.226, -5.317, 1.648 ],
      [ 0., 2.582, .829, -1., 9.375, .333, .543, -1.376, -.049, .036 ],
      [ 0., 10.05, 8.204, -1., 151.7, 254.5, -11.46, -13.16, -31.68, 2.408 ],
      [ 0., 5.303, 5.918, -1., 16.58, -4.292, -72.89, -7.881, -13.76, 3.206 ],
      [ 0., -14.6, 3.125, 1., 251., 15.63, 72.24, 18.2, 12.21, .781 ],
      [ 0., 3.997, -1., 25.83, -3.268, 0., 0., 0., 0., 0. ],
      [ 0., 3.059, -1., 14.76, -2.049, 0., 0., 0., 0., 0. ],
      [ 0., 8.177, -1., 105.5, -10.1, 0., 0., 0., 0., 0. ],
      [ 0., -6.144, 1., 52.54, 6.711, 0., 0., 0., 0., 0. ],
      [ 0., -1.33, 1., 5.919, 2.475, 0., 0., 0., 0., 0. ]]

def calc_Lmin(M, X, FC_row):
    logM = np.log10(M)
    F = FC[FC_row]
    return F[1] + F[2]*X + (F[3] + F[4]*X)*logM + (F[5] + F[6]*X)*logM**2

def calc_LHe(M, FC_row):
    logM = np.log10(M)
    F = FC[FC_row]
    return F[1] + F[2]*logM + F[3]*logM**2

def calc_Mmax(logL, X, FC_row):
    F = FC[FC_row]
    f = F[4] + F[5]*X + F[6]*X**2 + (F[7] + F[8]*X)*logL
    return (F[1] + F[2]*X + F[3]*np.sqrt(f)) / (1 + F[9]*X)

def calc_MHe(logL, FC_row):
    F = FC[FC_row]
    return F[1] + F[2]*np.sqrt(F[3] + F[4]*logL)

def main():
    warnings = []
    choice = input("Choose: 1 for Luminosity Calculator, 2 for Mass Calculator: ").strip()

    if choice == "1":
        print("\nInputs:")
        try:
            M = float(input("  Mass (Msun): "))
            X = float(input("  X (mass fraction): "))
        except ValueError:
            print("\nError:\n  All inputs must be valid numbers.")
            return
        if M <= 0 or X < 0:
            print("\nError:\n  Yea, nice try :) Zero or negative input(s).")
            return
        if X > 1:
            print("\nError:\n  Yea, nice try :) X > 1")
            return

        print("\nOutputs:")
        if X == 0:
            within_range = False
            if 8 <= M <= 250:
                print(f"  log(L_He/Lsun):  {calc_LHe(M, 6):.5f}, slope: inf [8–250] (Sequence no. 6)")
                within_range = True
            if 0.6 <= M <= 100:
                print(f"  log(L_He/Lsun):  {calc_LHe(M, 7):.5f}, slope: inf [0.6–100] (Sequence no. 7)")
                within_range = True
            if 60 <= M <= 1000:
                print(f"  log(L_He/Lsun):  {calc_LHe(M, 8):.5f}, slope: inf [60–1000] (Sequence no. 8)")
                within_range = True
            if not within_range:
                if M < 0.6:
                    print(f"  log(L_He/Lsun):  {calc_LHe(M, 7):.5f}, slope: inf [0.6–100] (Sequence no. 7)")
                elif M > 1000:
                    print(f"  log(L_He/Lsun):  {calc_LHe(M, 8):.5f}, slope: inf [60–1000] (Sequence no. 8)")
        else:
            within_range = False
            if 12 <= M <= 250:
                print(f"  log(L_min/Lsun): {calc_Lmin(M, X, 1):.5f}, slope: 0 [12–250] (Sequence no. 1)")
                within_range = True
            if 2 <= M <= 100:
                print(f"  log(L_min/Lsun): {calc_Lmin(M, X, 2):.5f}, slope: 0 [2–100] (Sequence no. 2)")
                within_range = True
            if 60 <= M <= 4000:
                print(f"  log(L_min/Lsun): {calc_Lmin(M, X, 3):.5f}, slope: 0 [60–4000] (Sequence no. 3)")
                within_range = True
            if not within_range:
                if M < 2:
                    print(f"  log(L_min/Lsun): {calc_Lmin(M, X, 2):.5f}, slope: 0 [2–100] (Sequence no. 2)")
                elif M > 4000:
                    print(f"  log(L_min/Lsun): {calc_Lmin(M, X, 3):.5f}, slope: 0 [60–4000] (Sequence no. 3)")

            within_range = False
            if 8 <= M <= 250:
                print(f"  log(L_He/Lsun):  {calc_LHe(M, 6):.5f}, slope: inf [8–250] (Sequence no. 6)")
                within_range = True
            if 0.6 <= M <= 100:
                print(f"  log(L_He/Lsun):  {calc_LHe(M, 7):.5f}, slope: inf [0.6–100] (Sequence no. 7)")
                within_range = True
            if 60 <= M <= 1000:
                print(f"  log(L_He/Lsun):  {calc_LHe(M, 8):.5f}, slope: inf [60–1000] (Sequence no. 8)")
                within_range = True
            if not within_range:
                if M < 0.6:
                    print(f"  log(L_He/Lsun):  {calc_LHe(M, 7):.5f}, slope: inf [0.6–100] (Sequence no. 7)")
                elif M > 1000:
                    print(f"  log(L_He/Lsun):  {calc_LHe(M, 8):.5f}, slope: inf [60–1000] (Sequence no. 8)")
                
        if M < 2 or M > 4000:
            warnings.append("Input mass is outside the grid range for L_min (2 ≤ M/Msun ≤ 4000)")
        if M < 0.6 or M > 1000:
            warnings.append("Input mass is outside the grid range for L_He (0.6 ≤ M/Msun ≤ 1000)")
        if X > 0.7:
            warnings.append("Input X is outside the grid range (0 ≤ X ≤ 0.7)")


    elif choice == "2":
        print("\nInputs:")
        try:
            L = float(input("  Luminosity (logL): "))
            X = float(input("  X (mass fraction): "))
        except ValueError:
            print("\nError:\n  All inputs must be valid numbers.")
            return
        if X < 0:
            print("\nError:\n  Yea, nice try :) Zero or negative input(s).")
            return
        if X > 1:
            print("\nError:\n  Yea, nice try :) X > 1")
            return

        print("\nOutputs:")
        if X == 0:  
            within_range = False
            if 8 <= 10**calc_MHe(L, 16) <= 250:
                print(f"  log(M_He/Msun): {10**calc_MHe(L, 16):.5f}, slope: inf [8–250] (Sequence no. 16)")
                within_range = True 
            if 0.6 <= 10**calc_MHe(L, 17) <= 100:
                print(f"  log(M_He/Msun): {10**calc_MHe(L, 17):.5f}, slope: inf [0.6–100] (Sequence no. 17)")
                within_range = True
            if 60 <= 10**calc_MHe(L, 18) <= 1000:
                print(f"  log(M_He/Msun): {10**calc_MHe(L, 18):.5f}, slope: inf [60–1000] (Sequence no. 18)")
                within_range = True
            if not within_range:
                if 10**calc_MHe(L, 17) < 0.6:
                    print(f"  log(M_He/Msun): {10**calc_MHe(L, 17):.5f}, slope: inf [0.6–100] (Sequence no. 17)")
                elif 10**calc_MHe(L, 18) > 1000:
                    print(f"  log(M_He/Msun): {10**calc_MHe(L, 18):.5f}, slope: inf [60–1000] (Sequence no. 18)")                   
        else:
            within_range = False
            if 12 <= 10**calc_Mmax(L, X, 11) <= 250:
                print(f"  log(M_max/Msun): {10**calc_Mmax(L, X, 11):.5f}, slope: 0 [12–250] (Sequence no. 11)")
                within_range = True
            if 2 <= 10**calc_Mmax(L, X, 12) <= 100:
                print(f"  log(M_max/Msun): {10**calc_Mmax(L, X, 12):.5f}, slope: 0 [2–100] (Sequence no. 12)")
                within_range = True
            if 60 <= 10**calc_Mmax(L, X, 13) <= 4000:
                print(f"  log(M_max/Msun): {10**calc_Mmax(L, X, 13):.5f}, slope: 0 [60–4000] (Sequence no. 13)")
                within_range = True
            if not within_range:
                if 10**calc_Mmax(L, X, 12) < 2:
                    print(f"  log(M_max/Msun): {10**calc_Mmax(L, X, 12):.5f}, slope: 0 [2–100] (Sequence no. 12)")
                elif 10**calc_Mmax(L, X, 13) > 4000:
                    print(f"  log(M_min/Msun): {10**calc_Mmax(L, X, 13):.5f}, slope: 0 [60–4000] (Sequence no. 13)")
                    
            within_range = False
            if 8 <= 10**calc_MHe(L, 16) <= 250:
                print(f"  log(M_He/Msun): {10**calc_MHe(L, 16):.5f}, slope: inf [8–250] (Sequence no. 16)")
                within_range = True 
            if 0.6 <= 10**calc_MHe(L, 17) <= 100:
                print(f"  log(M_He/Msun): {10**calc_MHe(L, 17):.5f}, slope: inf [0.6–100] (Sequence no. 17)")
                within_range = True
            if 60 <= 10**calc_MHe(L, 18) <= 1000:
                print(f"  log(M_He/Msun): {10**calc_MHe(L, 18):.5f}, slope: inf [60–1000] (Sequence no. 18)")
                within_range = True
            if not within_range:
                if 10**calc_MHe(L, 17) < 0.6:
                    print(f"  log(M_He/Msun): {10**calc_MHe(L, 17):.5f}, slope: inf [0.6–100] (Sequence no. 17)")
                elif 10**calc_MHe(L, 18) > 1000:
                    print(f"  log(M_He/Msun): {10**calc_MHe(L, 18):.5f}, slope: inf [60–1000] (Sequence no. 18)")  
                    
        if 10**calc_Mmax(L, X, 12) < 2 or 10**calc_Mmax(L, X, 13) > 4000:
            warnings.append("Output M_max is outside the grid range (2 ≤ M/Msun ≤ 4000)")
        if 10**calc_MHe(L, 17) < 0.6 or 10**calc_MHe(L, 18) > 1000:
            warnings.append("Output M_He is outside the grid range (0.6 ≤ M/Msun ≤ 1000)")
        if X > 0.7:
            warnings.append("Input X is outside the grid range (0 ≤ X ≤ 0.7)")
    else:
        print("Invalid choice.")
        return

    if warnings:
        print("\nWarning(s):")
        for w in warnings:
            print(" ", w)

if __name__ == "__main__":
    main()
