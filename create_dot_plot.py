"""
sequencing/main.py

Summary
-------
This small script converts a DNA sequence (string) into a square binary matrix
where element (i,j) == 1 when the characters at positions i and j match, else 0.

Provided utilities
- strToMatrix(s): returns a square list-of-lists with 1/0 entries for equality
- cleanliness(r, c, arr): returns the sum of values in a 5x5 neighborhood
- cleanliness_x_fraction(r, c, arr): returns fraction of nonzero cells in a 7x7
    window centered at (r,c) that belong to a thick "X" shape (see docstring)
- plot_binary_matrix(arr, out=.., show=..): plot/save a black/white image

Usage
-----
Run the module directly; by default it saves `matrix.png` in the working folder.
Requires: numpy, matplotlib

Notes
-----
The X-shape contains cells in the 7x7 window whose absolute row/col offsets
from the center differ by at most 1 (a thick diagonal/X). The fraction returned
is (# nonzero cells in that X) / (total X cells in the available window).
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

def same(let1, let2):
    return 1 if let1 == let2 else 0

def strToMatrix(s: str):
    """Return a square 2D list where entry (i,j) == 1 if s[i]==s[j], else 0."""
    arr = [[same(cLet, rLet) for cLet in s] for rLet in s]
    return arr


def cleaningFunction(r, c, arr, radius=3):
    """Return fraction of cells in a (2*radius+1)x(2*radius+1) window centered at
    (r,c) that are part of an "X" shape and are nonzero in `arr`.

    The X shape is defined as positions whose row and column offsets from the
    center are within 1 of each other: abs(abs(dr) - abs(dc)) <= 1. This
    gives a 1-cell thick diagonal (both diagonals) with a small thickness.

    - r, c: center cell indices
    - arr: 2D array-like (list of lists or numpy array)
    - radius: half-width of the square window (default 3 -> 7x7 window)

    Returns a float in [0.0, 1.0]. If no X-shape cells fall inside the array
    bounds, returns 0.0.
    """
    a = np.array(arr)
    rows, cols = a.shape
    rmin = max(r - radius, 0)
    rmax = min(r + radius, rows - 1)
    cmin = max(c - radius, 0)
    cmax = min(c + radius, cols - 1)

    total_in_x = 0
    nonzero_in_x = 0
    for i in range(rmin, rmax + 1):
        for j in range(cmin, cmax + 1):
            dr = abs(i - r)
            dc = abs(j - c)
            # X-shape membership: row and column offsets within 1 of each other
            if abs(dr - dc) <= 1:
                total_in_x += 1
                if a[i, j] != 0:
                    nonzero_in_x += 1

    if total_in_x == 0:
        return 0.0
    return float(nonzero_in_x) / float(total_in_x)


def plot_binary_matrix(arr, save_path=None, show=False, figsize=6):
    """Plot a 2D array-like (zeros and nonzeros). Nonzero -> black, zero -> white.

    - arr: 2D array-like (list of lists or numpy array)
    - save_path: if provided, save the figure to this path
    - show: if True, call plt.show()
    """
    def col(r, c):
        return 0 if arr[r][c] == 0 else 1 if cleaningFunction(r, c, arr) > 1/3 else 2
    cleanArr = [[col(r, c) for c in range(len(arr[0]))] for r in range(len(arr))]
    mask = np.array(cleanArr, dtype=int)
    cmap = ListedColormap(['white', 'black', 'red'])
    fig, ax = plt.subplots(figsize=(figsize, figsize))
    ax.imshow(mask, cmap=cmap, interpolation='nearest')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('Binary array: black=nonzero, white=zero')
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved image to {save_path}")
    if show:
        plt.show()
    plt.close(fig)


if __name__ == '__main__':
    # built-in default sequence (used when --seq is not provided)
    s_default = "TTGATTTTTGGACAGGAATACTATTCTGAGAAGATGCTCCGTATTTATTCTGGCCTCTTTATATTAAATAAAAATCAAAATAAAACTAAAATAAAATGAGGCATTTACATTGAAATCATGTTTTTACTCCCCCCATTTAAATGAGGTGTGG"

    parser = argparse.ArgumentParser(description='Plot binary matrix from string equality.')
    parser.add_argument('--seq', '-s', default=s_default,
                        help='Sequence string to use (overrides built-in default)')
    parser.add_argument('--out', '-o', default='matrix.png', help='Output image path (PNG)')
    parser.add_argument('--show', action='store_true', help='Show the image window')
    args = parser.parse_args()
    if set(args.seq) - set('ACGT'):
        parser.error('Error: sequence contains invalid characters; only A,C,G,T allowed.')
    
    # use a prefix so we get a manageable square matrix for plotting
    arr = strToMatrix(args.seq)

    # print a small summary
    print(f"Making matrix from segment length={len(args.seq)} -> matrix shape {len(arr)}x{len(arr[0])}")
    plot_binary_matrix(arr, save_path=args.out, show=args.show)

