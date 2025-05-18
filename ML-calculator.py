#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ----------- How to Use -------------

# Upon running the Python script, the user will be prompted to choose a calculator: Enter 1 for the Luminosity Calculator or 2 for the Mass Calculator.
# Based on the selected option, the user will be asked to input either the stellar mass or luminosity, along with the H and Z abundances (as mass fractions).
# The script will output the minimum, maximum, and pure-helium values of mass or luminosity, depending on the chosen calculator.
# Warnings will be displayed if any of the input or output values exceed the tested grid parameter ranges.


import numpy as np

Z1, Z2 = 0.008, 0.004

L_min_Z1 = [2.053491, 3.790927, -0.802070, -2.976704, 0.965973, 0.185089, 0.369268, -0.374144, 0.105449, 0.005]
L_max_Z1 = [3.751088, 2.209607, -0.453056, -0.520778, 0.245808, -0.016714, -1.329120, 1.228870, -0.262928, 0.005]
L_min_Z2 = [2.125432, 3.689468, -0.763519, -2.900812, 0.934060, 0.173159, 0.308744, -0.307890, 0.090761, 0.005]
L_max_Z2 = [3.733297, 2.198926, -0.424813, -0.552451, 0.309716, -0.060483, -1.305613, 1.228668, -0.286156, 0.005]


def interpolate_params(Z, p1, p2):
    return [p1[i] + (p2[i] - p1[i]) * (Z - Z1) / (Z2 - Z1) for i in range(len(p1))]

def calc_L(m, x, params):
    logm = np.log10(m)
    return sum(p * logm**i for i, p in enumerate(params[:3])) + \
           sum(p * logm**i for i, p in enumerate(params[3:6])) * x + \
           sum(p * logm**i for i, p in enumerate(params[6:9])) * np.exp(-x / params[9])

def get_min_max_L(m, x, Z):
    p_min = interpolate_params(Z, L_min_Z1, L_min_Z2)
    p_max = interpolate_params(Z, L_max_Z1, L_max_Z2)
    return calc_L(m, x, p_min), calc_L(m, x, p_max)

def bisection(f, a, b, tol=1e-6, max_iter=100):
    fa = f(a)
    fb = f(b)
    if fa * fb > 0:
        return None
    for _ in range(max_iter):
        c = (a + b) / 2.0
        fc = f(c)
        if abs(fc) < tol or (b - a) / 2.0 < tol:
            return round(c, 5)
        if fa * fc < 0:
            b, fb = c, fc
        else:
            a, fa = c, fc
    return round((a + b) / 2.0, 5)

def root_find_mass(L_val, x_val, m_low, m_high, target, Z):
    def luminosity_diff(m):
        L_min, L_max = get_min_max_L(m, x_val, Z)
        return (L_min if target == "L_min" else L_max) - L_val
    return bisection(luminosity_diff, m_low, m_high)

def main():
    warnings = []
    choice = input("Choose: 1 for Luminosity Calculator, 2 for Mass Calculator: ").strip()

    if choice == "1":
        print("\nInputs:")
        m = float(input("  M: "))
        x = float(input("  X: "))
        Z = float(input("  Z: "))

        if Z != Z1 and Z != Z2:
            if Z1 < Z < Z2 or Z2 < Z < Z1:
                warnings.append("Values are interpolated between Z = 0.008 and Z = 0.004")
            else:
                warnings.append("Values are extrapolated beyond Z = 0.008 and Z = 0.004")

        if x < 0 or x > 0.7:
            warnings.append("Warning: X is outside tested range (0 <= X <= 0.7)")
        if x + Z > 1:
            print("\nError:\n  X + Z > 1")
            return

        if m < 1 or m > 18:
            warnings.append("Input Mass is outside the tested range for Lmax (1 < M < 18)")
        if m < 1 or m > 40:
            warnings.append("Input Mass is outside the tested range for Lmin and LHe (1 < M < 40)")

        print("\nOutputs:")
        if x == 0:
            L_he = get_min_max_L(m, 0, Z)[1]
            print(f"\n  Pure He luminosity: {L_he:.5f}")
        else:
            L_min, L_max = get_min_max_L(m, x, Z)
            L_he = get_min_max_L(m, 0, Z)[1]
            print(f"  L_min: {L_min:.5f}")
            print(f"  L_max: {L_max:.5f}")
            print(f"  Pure He luminosity: {L_he:.5f}")

    elif choice == "2":
        print("\nInputs:")
        L_val = float(input("  L: "))
        x_val = float(input("  X: "))
        Z = float(input("  Z: "))

        if Z != Z1 and Z != Z2:
            if Z1 < Z < Z2 or Z2 < Z < Z1:
                warnings.append("Values are interpolated between Z = 0.008 and Z = 0.004")
            else:
                warnings.append("Values are extrapolated beyond Z = 0.008 and Z = 0.004")

        if x_val < 0 or x_val > 0.7:
            warnings.append("Warning: X is outside tested range (0 <= X <= 0.7)")
        if x_val + Z > 1:
            print("\nError:\n  X + Z > 1")
            return

        print("\nOutputs:")
        if x_val == 0:
            m_he = root_find_mass(L_val, 0, 0.5, 20, "L_max", Z)
            if m_he is not None and (m_he < 1 or m_he > 40):
                warnings.append("Output pure He mass is outside tested range (1 < M < 40)")
            print(f"\n  Pure He mass: {m_he}")
        else:
            m_max = root_find_mass(L_val, x_val, 0.5, 100, "L_min", Z)
            m_min = root_find_mass(L_val, x_val, 0.5, 50, "L_max", Z)
            m_he = root_find_mass(L_val, 0, 0.5, 100, "L_max", Z)

            if m_min is not None and (m_min < 1 or m_min > 18):
                warnings.append("Output M_min is outside tested range (1 < M < 18)")
            if m_max is not None and (m_max < 1 or m_max > 40):
                warnings.append("Output M_max is outside tested range (1 < M < 40)")
            if m_he is not None and (m_he < 1 or m_he > 40):
                warnings.append("Output pure He mass is outside tested range (1 < M < 40)")

            print(f"  M_min: {m_min}")
            print(f"  M_max: {m_max}")
            print(f"  Pure He mass: {m_he}")
    else:
        print("Invalid choice.")
        return

    if warnings:
        print("\nWarning(s):")
        for w in warnings:
            print(" ", w)

if __name__ == "__main__":
    main()

