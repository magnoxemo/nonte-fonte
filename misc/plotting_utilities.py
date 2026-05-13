"""
plotting_utilities.py
This script is generate by Chatgpt
Generalized plotting utilities for nonte-fonte (local_fet) C++ library results.
"""

#!/usr/bin/env bash
import numpy as np
import matplotlib.pyplot as plt
import argparse
from numpy.polynomial.legendre import legval


# ----------------------------
# Domain mapping
# ----------------------------
def map_to_reference(x, a, b):
    return (2.0 * x - (a + b)) / (b - a)


# ----------------------------
# Legendre reconstruction
# ----------------------------
def reconstruct_1d(x, coeffs, domain):
    x0, x1 = domain
    xi = map_to_reference(x, x0, x1)
    return legval(xi, coeffs)


def reconstruct_2d(x, y, coeffs, domain):
    x0, x1, y0, y1 = domain
    xi = map_to_reference(x, x0, x1)
    eta = map_to_reference(y, y0, y1)

    Nx, Ny = coeffs.shape

    val = 0.0
    for i in range(Nx):
        Pi = legval(xi, [0]*i + [1])
        for j in range(Ny):
            Pj = legval(eta, [0]*j + [1])
            val += coeffs[i, j] * Pi * Pj

    return val


# ----------------------------
# Parsing FET output
# ----------------------------
def parse_fet_file(filename, is_2d=True, nx=None, ny=None):
    """
    Expected format:

    Domain = x0, x1, y0, y1   (or x0, x1)
    coeff rows...

    Returns list of (domain, coeffs)
    """
    blocks = []

    with open(filename, "r") as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if line.startswith("Domain"):
            parts = line.replace("Domain =", "").split(",")

            if is_2d:
                x0, x1, y0, y1 = map(float, parts)
                domain = (x0, x1, y0, y1)
            else:
                x0, x1 = map(float, parts)
                domain = (x0, x1)

            i += 1
            coeffs = []

            while i < len(lines) and not lines[i].startswith("Domain"):
                row = lines[i].strip().replace(",", "")
                if row:
                    coeffs.extend([float(v) for v in row.split()])
                i += 1

            coeffs = np.array(coeffs)

            if is_2d:
                if nx is None or ny is None:
                    raise ValueError("nx, ny must be provided for 2D case")
                coeffs = coeffs.reshape(nx, ny)

            blocks.append((domain, coeffs))
        else:
            i += 1
    return blocks

def plot_1d(blocks, num_points=300):
    plt.figure()

    for domain, coeffs in blocks:
        x0, x1 = domain
        x = np.linspace(x0, x1, num_points)

        y = np.array([reconstruct_1d(xi, coeffs, (x0, x1)) for xi in x])

        plt.plot(x, y, label=f"[{x0:.2f},{x1:.2f}]")

    plt.title("FET Legendre Reconstruction (1D)")
    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


def plot_2d(blocks, num_points=60):
    for domain, coeffs in blocks:
        x0, x1, y0, y1 = domain

        x = np.linspace(x0, x1, num_points)
        y = np.linspace(y0, y1, num_points)

        X, Y = np.meshgrid(x, y)
        Z = np.zeros_like(X)

        for i in range(num_points):
            for j in range(num_points):
                Z[j, i] = reconstruct_2d(X[j, i], Y[j, i], coeffs, domain)

        plt.figure()
        plt.contourf(X, Y, Z, levels=30)
        plt.colorbar()
        plt.title(f"FET Cell [{x0},{x1}] x [{y0},{y1}]")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.tight_layout()
        plt.show()


# ----------------------------
# CLI
# ----------------------------
def main():
    parser = argparse.ArgumentParser(description="FET Legendre reconstruction plotter")

    parser.add_argument("file", type=str, help="FET output text file")
    parser.add_argument("--dim", type=int, default=2, choices=[1, 2], help="Dimension: 1 or 2")
    parser.add_argument("--nx", type=int, default=None, help="Legendre order x (2D only)")
    parser.add_argument("--ny", type=int, default=None, help="Legendre order y (2D only)")
    parser.add_argument("--points", type=int, default=80, help="Resolution")

    args = parser.parse_args()

    if args.dim == 1:
        blocks = parse_fet_file(args.file, is_2d=False)
        plot_1d(blocks, num_points=args.points)

    else:
        if args.nx is None or args.ny is None:
            raise ValueError("For 2D, you must specify --nx and --ny")

        blocks = parse_fet_file(args.file, is_2d=True, nx=args.nx, ny=args.ny)
        plot_2d(blocks, num_points=args.points)


if __name__ == "__main__":
    main()