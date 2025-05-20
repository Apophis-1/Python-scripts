#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ---------------------------------  How to Use ----------------------------------
# When the script is executed, the user is prompted to choose a calculator: enter 1 for the 
# Luminosity Calculator or 2 for the Mass Calculator. Based on the selected option, the user 
# must input either the stellar mass or luminosity, along with the hydrogen (X) and metal (Z) 
# abundances as mass fractions. The script then outputs the minimum, maximum, and pure-helium 
# values of either mass or luminosity, depending on the chosen calculator.

# Errors are displayed if the inputs are not valid numbers, or if the mass is zero or negative, 
# or if X or Z is negative. X = 0 and Z = 0 are allowed.

# A set of warnings is printed based on the parameter range of the synthetic model grid. If the 
# inputs fall outside the grid’s tested parameter range, a general warning is shown. If the inputs 
# are significantly beyond the grid range such that the minimum or maximum value of M or L is not 
# truly a minimum or maximum, then a warning is issued indicating that the ML fits may be 
# unreliable. 

# The model grid was computed for Z = 0.008 and Z = 0.004, which approximately correspond
# to LMC-like (0.4 Zsun) and SMC-like (0.2 Zsun) metallicities, where Zsun = 0.02. For 
# any Z value other than 0.008 or 0.004, interpolation or extrapolation is performed, and a 
# corresponding warning is provided.





import numpy as np

Z1, Z2 = 0.008, 0.004

L_min_Z1 = [2.053491, 3.790927, -0.802070, -2.976704, 0.965973, 0.185089, 0.369268, -0.374144, 0.105449, 0.005]
L_max_Z1 = [3.751088, 2.209607, -0.453056, -0.520778, 0.245808, -0.016714, -1.329120, 1.228870, -0.262928, 0.005]
L_min_Z2 = [2.125432, 3.689468, -0.763519, -2.900812, 0.934060, 0.173159, 0.308744, -0.307890, 0.090761, 0.005]
L_max_Z2 = [3.733297, 2.198926, -0.424813, -0.552451, 0.309716, -0.060483, -1.305613, 1.228668, -0.286156, 0.005]

s_LMC = [0.698967, -0.025170, 0.003576, 5.017684, -1.125765, 1.362459, -2.995227, 1.177010, -0.692827]
s_SMC = [0.709244, 0.007519, -0.020923, 4.636537, -1.631714, 1.744423, -2.711640, 2.088682, -1.369483]

def calc_L(M, X, params):
    logm = np.log10(M)
    return sum(params[i] * logm**i for i in range(3)) + \
           sum(params[i] * logm**(i-3) for i in range(3,6)) * X + \
           sum(params[i] * logm**(i-6) for i in range(6,9)) * np.exp(-X / params[9])

def get_L_values(M, X, Z):
    factor = (Z - Z1) / (Z2 - Z1)
    L_min = calc_L(M, X, L_min_Z1) + factor * (calc_L(M, X, L_min_Z2) - calc_L(M, X, L_min_Z1))
    L_max = calc_L(M, X, L_max_Z1) + factor * (calc_L(M, X, L_max_Z2) - calc_L(M, X, L_max_Z1))
    L_he  = calc_L(M, 0, L_max_Z1) + factor * (calc_L(M, 0, L_max_Z2) - calc_L(M, 0, L_max_Z1))
    return L_min, L_max, L_he

def get_slope(M, X, Z):
    logm = np.log10(M)
    def s_formula(coeffs):
        f1, f2, f3, f4, f5, f6, f7, f8, f9 = coeffs
        return (f1 + f2 * logm + f3 * logm**2) + \
               (f4 + f5 * logm + f6 * logm**2) * X + \
               (f7 + f8 * logm + f9 * logm**2) * X**2
    s1 = s_formula(s_LMC)
    s2 = s_formula(s_SMC)
    factor = (Z - Z1) / (Z2 - Z1)
    return s1 + factor * (s2 - s1)

def root_find_mass_for_Z(L_val, x_val, m_low, m_high, target, Z_fixed):
    def luminosity_diff(m):
        L_min, L_max, _ = get_L_values(m, x_val, Z_fixed)
        return (L_min if target == "L_min" else L_max) - L_val
    return bisection(luminosity_diff, m_low, m_high)

def root_find_mass(L_val, x_val, m_low, m_high, target, Z):
    m_at_Z1 = root_find_mass_for_Z(L_val, x_val, m_low, m_high, target, Z1)
    m_at_Z2 = root_find_mass_for_Z(L_val, x_val, m_low, m_high, target, Z2)
    if m_at_Z1 is None or m_at_Z2 is None:
        return None
    factor = (Z - Z1) / (Z2 - Z1)
    return m_at_Z1 + factor * (m_at_Z2 - m_at_Z1)

def bisection(f, a, b, tol=1e-6, max_iter=100):
    fa, fb = f(a), f(b)
    if fa * fb > 0:
        return None
    for _ in range(max_iter):
        c = (a + b) / 2
        fc = f(c)
        if abs(fc) < tol or (b - a) / 2 < tol:
            return round(c, 5)
        if fa * fc < 0:
            b, fb = c, fc
        else:
            a, fa = c, fc
    return round((a + b) / 2, 5)

def main():
    warnings = []
    choice = input("Choose: 1 for Luminosity Calculator, 2 for Mass Calculator: ").strip()

    if choice == "1":
        print("\nInputs:")
        try:
            M = float(input("  Mass (Msun): "))
            X = float(input("  X (mass fraction): "))
            Z = float(input("  Z (mass fraction): "))
        except ValueError:
            print("\nError:\n  All inputs must be valid numbers.")
            return
        if M <= 0 or X < 0 or Z < 0:
            print("\nError:\n  Yea, nice try :) Zero or negative input(s).")
            return
        if X + Z > 1:
            print("\nError:\n  Yea, nice try :) X + Z > 1")
            return
        
        s = get_slope(M, X, Z)

        print("\nOutputs:")
        if X == 0:
            L_he = get_L_values(M, 0, Z)[1]
            print(f"  log(L_He/Lsun):  {L_he:.5f}, slope: inf")
        else:
            L_min, L_max, L_he = get_L_values(M, X, Z)
            print(f"  log(L_min/Lsun): {L_min:.5f}, slope: 0")
            print(f"  log(L_max/Lsun): {L_max:.5f}, slope: {s:.2f}")
            print(f"  log(L_He/Lsun):  {L_he:.5f}, slope: inf")

        if L_min > L_max or L_min > L_he or L_max < L_he:
            print("\nWarning(s):\n One or more inputs are well beyond the grid range. \n The fit calculations may not be reliable.")
            return
        if M < 1 or M > 18:
            warnings.append("Input mass is outside the grid range for L_max (1 ≤ M ≤ 18)")
        if M < 1 or M > 40:
            warnings.append("Input mass is outside the grid range for L_min and L_He (1 ≤ M ≤ 40)")

        if X > 0.7:
            warnings.append("Input X is outside grid range (0 ≤ X ≤ 0.7)")
        if Z != Z1 and Z != Z2:
            if min(Z1,Z2) < Z < max(Z1,Z2):
                warnings.append("Luminosity and slope values are interpolated between Z = 0.008 and 0.004")
            else:
                warnings.append("Luminosity and slope values are extrapolated beyond Z = 0.008 and 0.004")

    elif choice == "2":
        print("\nInputs:")
        try:
            L = float(input("  Luminosity (logL): "))
            X = float(input("  X (mass fraction): "))
            Z = float(input("  Z (mass fraction): "))
        except ValueError:
            print("\nError:\n  All inputs must be valid numbers.")
            return
        if X < 0 or Z < 0:
            print("\nError:\n  Yea, nice try :) Zero or negative input(s).")
            return
        if X + Z > 1:
            print("\nError:\n  Yea, nice try :) X + Z > 1")
            return

        print("\nOutputs:")
        if X == 0:
            m_he = root_find_mass(L, 0, 0.01, 100, "L_max", Z)
            if m_he is not None and (m_he < 1 or m_he > 40):
                warnings.append("Output M_He is outside grid range (1 ≤ M ≤ 40)")
            print(f"  M_He:  {m_he:.5f}, slope: inf")
        else:
            m_max = root_find_mass(L, X, 0.01, 100, "L_min", Z)
            m_min = root_find_mass(L, X, 0.01, 50, "L_max", Z)
            m_he = root_find_mass(L, 0, 0.01, 100, "L_max", Z)
            s = get_slope(m_min, X, Z) if m_min is not None else None
            if None in (m_min, m_max, m_he):
                print("\nError: \n One or more inputs are well beyond the grid range. \n The fit calculations failed.")
                return
            else:
                print(f"  M_min/Msun: {m_min:.5f}, slope: {s:.2f}")
                print(f"  M_max/Msun: {m_max:.5f}, slope: 0")
                print(f"  M_He/Msun:  {m_he:.5f}, slope: inf")

        if m_min is not None and (m_min > m_max or m_min > m_he or m_max < m_he):
            print("\nWarning(s):\n One or more inputs are well beyond the grid range. \n The fit calculations may not be reliable.")
            return

        if m_min is not None and (m_min < 1 or m_min > 18):
            warnings.append("Output M_min is outside grid range (1 ≤ M ≤ 18)")
        if m_max is not None and (m_max < 1 or m_max > 40):
            warnings.append("Output M_max is outside grid range (1 ≤ M ≤ 40)")
        if m_he is not None and (m_he < 1 or m_he > 40):
            warnings.append("Output M_He is outside grid range (1 ≤ M ≤ 40)")
        if X > 0.7:
            warnings.append("Input X is outside grid range (0 ≤ X ≤ 0.7)")
        if Z != Z1 and Z != Z2:
            if min(Z1,Z2) < Z < max(Z1,Z2):
                warnings.append("Mass and slope values are interpolated between Z = 0.008 and 0.004")
            else:
                warnings.append("Mass and slope values are extrapolated beyond Z = 0.008 and 0.004")

    else:
        print("Invalid choice.")
        return

    if warnings:
        print("\nWarning(s):")
        for w in warnings:
            print(" ", w)

if __name__ == "__main__":
    main()
