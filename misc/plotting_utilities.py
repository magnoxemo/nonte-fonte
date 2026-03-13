"""
plotting_utilities.py
This script is generate by Chatgpt
Generalized plotting utilities for nonte-fonte (local_fet) C++ library results.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from typing import Optional, List, Callable
import os


def plot_1d_results(csv_file: str,
                    title: str = "1D FET Results",
                    figsize: tuple = (10, 6),
                    save_name: Optional[str] = None):
    """
    Plot 1D results from CSV file.

    Expected columns: x, true_pdf, [tally_0, tally_1, ...]
    """
    df = pd.read_csv(csv_file)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, sharex=True)

    # Plot PDFs
    ax1.plot(df['x'], df['true_pdf'], 'k-', linewidth=2, label='True PDF', alpha=0.7)

    # Plot tallies
    tally_cols = [col for col in df.columns if col.startswith('tally_')]
    colors = plt.cm.tab10(np.linspace(0, 1, len(tally_cols)))

    for i, col in enumerate(tally_cols):
        # Only plot non-zero values
        mask = df[col] != 0.0
        if mask.any():
            ax1.plot(df[mask]['x'], df[mask][col], '--',
                     color=colors[i], linewidth=2, label=col.replace('_', ' ').title())

    ax1.set_ylabel('PDF', fontsize=12)
    ax1.set_title(title, fontsize=14, fontweight='bold')
    ax1.legend(loc='best', fontsize=10)
    ax1.grid(True, alpha=0.3)

    # Plot errors
    for i, col in enumerate(tally_cols):
        mask = df[col] != 0.0
        if mask.any():
            error = df[mask]['true_pdf'] - df[mask][col]
            ax2.plot(df[mask]['x'], error, '-', color=colors[i],
                     linewidth=1.5, label=col.replace('_', ' ').title())

    ax2.axhline(0, color='k', linestyle='-', linewidth=1, alpha=0.5)
    ax2.set_xlabel('x', fontsize=12)
    ax2.set_ylabel('Error', fontsize=12)
    ax2.set_title('Reconstruction Error', fontsize=12)
    ax2.legend(loc='best', fontsize=10)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_name:
        plt.savefig(save_name, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {save_name}")

    plt.show()


def plot_2d_results(csv_file: str,
                    title: str = "2D FET Results",
                    figsize: tuple = (14, 5),
                    save_name: Optional[str] = None):
    """
    Plot 2D results from CSV file.

    Expected columns: x, y, true_pdf, [legendre_pdf, ...]
    """
    df = pd.read_csv(csv_file)

    # Pivot for contour plots
    x_unique = np.sort(df['x'].unique())
    y_unique = np.sort(df['y'].unique())

    X, Y = np.meshgrid(x_unique, y_unique)
    Z_true = df.pivot(index='y', columns='x', values='true_pdf').values

    # Find reconstruction columns
    recon_cols = [col for col in df.columns if col.endswith('_pdf') and col != 'true_pdf']

    n_plots = 1 + len(recon_cols)
    fig, axes = plt.subplots(1, n_plots, figsize=figsize)

    if n_plots == 1:
        axes = [axes]

    # True PDF
    im = axes[0].contourf(X, Y, Z_true, levels=20, cmap='viridis')
    axes[0].set_title('True PDF', fontsize=12, fontweight='bold')
    axes[0].set_xlabel('x')
    axes[0].set_ylabel('y')
    plt.colorbar(im, ax=axes[0])

    # Reconstructions
    for i, col in enumerate(recon_cols, 1):
        Z_recon = df.pivot(index='y', columns='x', values=col).values
        im = axes[i].contourf(X, Y, Z_recon, levels=20, cmap='viridis')
        axes[i].set_title(col.replace('_', ' ').title(), fontsize=12, fontweight='bold')
        axes[i].set_xlabel('x')
        axes[i].set_ylabel('y')
        plt.colorbar(im, ax=axes[i])

    plt.suptitle(title, fontsize=14, fontweight='bold')
    plt.tight_layout()

    if save_name:
        plt.savefig(save_name, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {save_name}")

    plt.show()


def compute_statistics(csv_file: str, tally_columns: Optional[List[str]] = None):
    """
    Compute and print statistics for tally results.

    Computes L2 error, max error, RMS error for each tally.
    """
    df = pd.read_csv(csv_file)

    if tally_columns is None:
        tally_columns = [col for col in df.columns if col.startswith('tally_')]

    print("\n" + "="*60)
    print("RECONSTRUCTION STATISTICS")
    print("="*60)

    for col in tally_columns:
        # Only compute for non-zero values
        mask = df[col] != 0.0
        if not mask.any():
            continue

        error = df[mask]['true_pdf'] - df[mask][col]

        # Compute metrics
        l2_error = np.sqrt(np.mean(error**2))
        max_error = np.max(np.abs(error))
        mean_error = np.mean(error)

        print(f"\n{col.replace('_', ' ').title()}:")
        print(f"  L2 Error:      {l2_error:.6e}")
        print(f"  Max Error:     {max_error:.6e}")
        print(f"  Mean Error:    {mean_error:.6e}")
        print(f"  RMS Error:     {l2_error:.6e}")

    print("\n" + "="*60)


def plot_histogram_results(csv_file: str,
                           figsize: tuple = (12, 5),
                           save_name: Optional[str] = None):
    """
    Plot histogram and FET comparison.

    Expected columns: x, true_pdf, histogram, [tally1, tally2, ...]
    """
    df = pd.read_csv(csv_file)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    # Left: PDFs
    ax1.plot(df['x'], df['true_pdf'], 'k-', linewidth=2, label='True PDF', alpha=0.7)

    if 'histogram' in df.columns:
        ax1.step(df['x'], df['histogram'], 'r-', linewidth=1.5,
                 where='mid', label='Histogram', alpha=0.7)

    tally_cols = [col for col in df.columns
                  if col.startswith('tally') and col != 'histogram']
    colors = plt.cm.tab10(np.linspace(0, 1, len(tally_cols)))

    for i, col in enumerate(tally_cols):
        mask = df[col] != 0.0
        if mask.any():
            ax1.plot(df[mask]['x'], df[mask][col], '--',
                     color=colors[i], linewidth=2, label=col.replace('_', ' ').title())

    ax1.set_xlabel('x', fontsize=12)
    ax1.set_ylabel('PDF', fontsize=12)
    ax1.set_title('PDF Comparison', fontsize=14, fontweight='bold')
    ax1.legend(loc='best')
    ax1.grid(True, alpha=0.3)

    # Right: Errors
    if 'histogram' in df.columns:
        hist_error = df['true_pdf'] - df['histogram']
        ax2.plot(df['x'], hist_error, 'r-', linewidth=1.5,
                 label='Histogram Error', alpha=0.7)

    for i, col in enumerate(tally_cols):
        mask = df[col] != 0.0
        if mask.any():
            error = df[mask]['true_pdf'] - df[mask][col]
            ax2.plot(df[mask]['x'], error, '-', color=colors[i],
                     linewidth=1.5, label=f"{col} Error")

    ax2.axhline(0, color='k', linestyle='-', linewidth=1, alpha=0.5)
    ax2.set_xlabel('x', fontsize=12)
    ax2.set_ylabel('Error', fontsize=12)
    ax2.set_title('Reconstruction Error', fontsize=14, fontweight='bold')
    ax2.legend(loc='best')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_name:
        plt.savefig(save_name, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {save_name}")

    plt.show()


def plot_convergence(csv_files: List[str],
                     labels: List[str],
                     figsize: tuple = (10, 6),
                     save_name: Optional[str] = None):
    """
    Plot convergence study from multiple CSV files.

    Each CSV should have x, true_pdf, and tally columns.
    """
    fig, ax = plt.subplots(figsize=figsize)

    errors = []

    for csv_file, label in zip(csv_files, labels):
        df = pd.read_csv(csv_file)
        tally_col = [col for col in df.columns if col.startswith('tally_')][0]

        mask = df[tally_col] != 0.0
        if mask.any():
            error = df[mask]['true_pdf'] - df[mask][tally_col]
            l2_error = np.sqrt(np.mean(error**2))
            errors.append(l2_error)

    ax.semilogy(labels, errors, 'o-', linewidth=2, markersize=8)
    ax.set_xlabel('Configuration', fontsize=12)
    ax.set_ylabel('L2 Error', fontsize=12)
    ax.set_title('Convergence Study', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_name:
        plt.savefig(save_name, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {save_name}")

    plt.show()


# Quick use functions
def quick_plot_1d(csv_file: str):
    """Quick 1D plot with auto-naming"""
    base_name = os.path.splitext(csv_file)[0]
    plot_1d_results(csv_file, save_name=f"{base_name}_plot.png")
    compute_statistics(csv_file)


def quick_plot_2d(csv_file: str):
    """Quick 2D plot with auto-naming"""
    base_name = os.path.splitext(csv_file)[0]
    plot_2d_results(csv_file, save_name=f"{base_name}_plot.png")


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: python plotting_utilities.py <csv_file>")
        print("  Automatically detects 1D or 2D and plots accordingly")
        sys.exit(1)

    csv_file = sys.argv[1]

    # Auto-detect 1D or 2D
    df = pd.read_csv(csv_file)

    if 'y' in df.columns:
        print("Detected 2D data")
        quick_plot_2d(csv_file)
    else:
        print("Detected 1D data")
        quick_plot_1d(csv_file)
