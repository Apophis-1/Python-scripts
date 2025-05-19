#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ----------- How to Use -------------

# Upon running the Python script, the user will be prompted to choose a calculator: Enter 1 for the Luminosity Calculator or 2 for the Mass Calculator.
# Based on the selected option, the user will be asked to input either the stellar mass or luminosity, along with the H and Z abundances (as mass fractions).
# The script will output the minimum, maximum, and pure-helium values of mass or luminosity, depending on the chosen calculator.
# Warnings will be displayed if any of the input or output values exceed the synthetic model grid parameter ranges.
# The synthetic model grid was run for Z = 0.008 and 0.004, roughly corresponding to LMC-like (0.4Zsun) and SMC-like (0.2Zsun) metallicity where Zsun = 0.02.
# For all values of input Z other than Z = 0.008 and 0.004, interpolation or extrapolation is performed.

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

Z1, Z2 = 0.008, 0.004

L_min_Z1 = [2.053491, 3.790927, -0.802070, -2.976704, 0.965973, 0.185089, 0.369268, -0.374144, 0.105449, 0.005]
L_max_Z1 = [3.751088, 2.209607, -0.453056, -0.520778, 0.245808, -0.016714, -1.329120, 1.228870, -0.262928, 0.005]
L_min_Z2 = [2.125432, 3.689468, -0.763519, -2.900812, 0.934060, 0.173159, 0.308744, -0.307890, 0.090761, 0.005]
L_max_Z2 = [3.733297, 2.198926, -0.424813, -0.552451, 0.309716, -0.060483, -1.305613, 1.228668, -0.286156, 0.005]

s_LMC = [0.698967, -0.025170, 0.003576, 5.017684, -1.125765, 1.362459, -2.995227, 1.177010, -0.692827]
s_SMC = [0.709244, 0.007519, -0.020923, 4.636537, -1.631714, 1.744423, -2.711640, 2.088682, -1.369483]

def calc_L(m, x, params):
    logm = np.log10(m)
    return sum(params[i] * logm**i for i in range(3)) + \
           sum(params[i] * logm**(i-3) for i in range(3,6)) * x + \
           sum(params[i] * logm**(i-6) for i in range(6,9)) * np.exp(-x / params[9])

def get_L_values(m, x, Z):
    factor = (Z - Z1) / (Z2 - Z1)

    L_min = calc_L(m, x, L_min_Z1) + factor * (calc_L(m, x, L_min_Z2) - calc_L(m, x, L_min_Z1))
    L_max = calc_L(m, x, L_max_Z1) + factor * (calc_L(m, x, L_max_Z2) - calc_L(m, x, L_max_Z1))
    L_he  = calc_L(m, 0, L_max_Z1) + factor * (calc_L(m, 0, L_max_Z2) - calc_L(m, 0, L_max_Z1))
    return L_min, L_max, L_he

def get_slope(m, x, Z):
    logm = np.log10(m)

    def s_formula(coeffs):
        f1, f2, f3, f4, f5, f6, f7, f8, f9 = coeffs
        return (f1 + f2 * logm + f3 * logm**2) + \
               (f4 + f5 * logm + f6 * logm**2) * x + \
               (f7 + f8 * logm + f9 * logm**2) * x**2
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
        m = float(input("  Mass (Msun): "))
        if m <= 0:
            print("\nError:\n  Input mass must be positive.")
            return
        x = float(input("  X (mass fraction): "))
        if x < 0:
            print("\nError:\n  Input X must be positive.")
            return
        Z = float(input("  Z (mass fraction): "))
        if Z < 0:
            print("\nError:\n  Input Z must be positive.")
            return

        if Z != Z1 and Z != Z2:
            if min(Z1,Z2) < Z < max(Z1,Z2):
                warnings.append("Luminosity values are interpolated between Z = 0.008 and Z = 0.004")
            else:
                warnings.append("Luminosity values are extrapolated beyond Z = 0.008 and Z = 0.004")

        if x < 0 or x > 0.7:
            warnings.append("Input X is outside grid range (0 ≤ X ≤ 0.7)")
        if x + Z > 1:
            print("\nError:\n  X + Z > 1")
            return

        if m < 1 or m > 18:
            warnings.append("Input mass is outside the grid range for L_max (1 ≤ M ≤ 18)")
        if m < 1 or m > 40:
            warnings.append("Input mass is outside the grid range for L_min and L_He (1 ≤ M ≤ 40)")

        s = get_slope(m, x, Z)

        print("\nOutputs:")
        if x == 0:
            L_he = get_L_values(m, 0, Z)[1]
            print(f"  L_He:  {L_he:.5f}, slope: inf")
        else:
            L_min, L_max, L_he = get_L_values(m, x, Z)
            print(f"  L_min: {L_min:.5f}, slope: 0")
            print(f"  L_max: {L_max:.5f}, slope: {s:.2f}")
            print(f"  L_He:  {L_he:.5f}, slope: inf")



    elif choice == "2":
        print("\nInputs:")
        L = float(input("  Luminosity (logL): "))
        x = float(input("  X (mass fraction): "))
        if x < 0:
            print("\nError:\n  Input X must be positive.")
            return
        Z = float(input("  Z (mass fraction): "))
        if Z < 0:
            print("\nError:\n  Input Z must be positive.")
            return

        if Z != Z1 and Z != Z2:
            if min(Z1,Z2) < Z < max(Z1,Z2):
                warnings.append("Mass values are interpolated between Z = 0.008 and Z = 0.004")
            else:
                warnings.append("Mass values are extrapolated beyond Z = 0.008 and 0.004")

        if x < 0 or x > 0.7:
            warnings.append("Input X is outside grid range (0 ≤ X ≤ 0.7)")
        if x + Z > 1:
            print("\nError:\n  X + Z > 1")
            return

        print("\nOutputs:")
        if x == 0:
            m_he = root_find_mass(L, 0, 0.5, 20, "L_max", Z)
            if m_he is not None and (m_he < 1 or m_he > 40):
                warnings.append("Output M_He is outside grid range (1 ≤ M ≤ 40)")
            print(f"  M_He:  {m_he:.5f}, slope: inf")
        else:
            m_max = root_find_mass(L, x, 0.5, 100, "L_min", Z)
            m_min = root_find_mass(L, x, 0.5, 50, "L_max", Z)
            m_he = root_find_mass(L, 0, 0.5, 100, "L_max", Z)

            s = get_slope(m_min, x, Z) if m_min is not None else None

            if m_min is not None and (m_min < 1 or m_min > 18):
                warnings.append("Output M_min is outside grid range (1 ≤ M ≤ 18)")
            if m_max is not None and (m_max < 1 or m_max > 40):
                warnings.append("Output M_max is outside grid range (1 ≤ M ≤ 40)")
            if m_he is not None and (m_he < 1 or m_he > 40):
                warnings.append("Output M_He is outside grid range (1 ≤ M ≤ 40)")

            print(f"  M_min: {m_min:.5f}, slope: {s:.2f}")
            print(f"  M_max: {m_max:.5f}, slope: 0")
            print(f"  M_He:  {m_he:.5f}, slope: inf")

    else:
        print("Invalid choice.")
        return

    if warnings:
        print("\nWarning(s):")
        for w in warnings:
            print(" ", w)

if __name__ == "__main__":
    main()
