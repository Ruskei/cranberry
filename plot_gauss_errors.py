#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot Gauss-law error metrics.")
    parser.add_argument(
        "metrics",
        nargs="?",
        type=Path,
        default=Path("out/gauss-law/metrics.npy"),
    )
    parser.add_argument("--save", type=Path, help="Save the figure instead of showing it.")
    args = parser.parse_args()

    data = np.load(args.metrics, allow_pickle=False)
    if data.ndim != 2 or data.shape[1] != 8:
        raise ValueError(f"expected an (N, 8) metrics array, got {data.shape}")
    if data.shape[0] == 0:
        raise ValueError("metrics array contains no samples")

    time = data[:, 0]
    fig, (relative, absolute) = plt.subplots(2, 1, sharex=True, figsize=(9, 7))

    relative.plot(time, data[:, 5], label="Relative L2")
    relative.plot(time, data[:, 7], label="Relative Linf")
    relative.set_ylabel("Error (%)")
    relative.grid(alpha=0.3)
    relative.legend()

    absolute.plot(time, np.abs(data[:, 3]), label="|Sum residual|")
    absolute.plot(time, data[:, 4], label="Residual L2")
    absolute.plot(time, data[:, 6], label="Residual Linf")
    absolute.set_xlabel("Simulation time")
    absolute.set_ylabel("Absolute error")
    absolute.grid(alpha=0.3)
    absolute.legend()

    fig.tight_layout()
    if args.save:
        fig.savefig(args.save, dpi=160)
    else:
        plt.show()


if __name__ == "__main__":
    main()
