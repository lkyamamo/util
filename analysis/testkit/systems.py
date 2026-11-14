"""Synthetic systems whose analysis results are known analytically.

Each generator states the answer it was built to have, so a test can compare an
analysis script's output to that answer instead of to a stored reference.
All randomness goes through a numpy Generator the caller seeds.
"""

import numpy as np


def ideal_gas(n_atoms, box_length, rng):
    """Uniform random positions in a cubic box [0, L)^3.

    Known answer: g(r) = 1 at every r (up to finite-N corrections of order 1/N),
    and S(q) = 1 for every q != 0.
    """
    return rng.uniform(0.0, box_length, size=(n_atoms, 3))


def random_walk(n_atoms, n_frames, dt, diffusion, rng):
    """Unwrapped trajectories of independent Gaussian random walks.

    Each step is N(0, 2*D*dt) per axis, so the known answer is
    MSD(t) = 6*D*t exactly in expectation at every lag t.
    Returns positions with shape (n_frames, n_atoms, 3), starting at the origin.
    """
    steps = rng.normal(0.0, np.sqrt(2.0 * diffusion * dt), size=(n_frames - 1, n_atoms, 3))
    return np.concatenate([np.zeros((1, n_atoms, 3)), np.cumsum(steps, axis=0)])


def harmonic_oscillators(n_atoms, n_frames, dt, omega, amplitude, rng):
    """Independent 3D harmonic oscillators at one angular frequency omega.

    x(t) = A cos(omega t + phi) per axis, with a random phase phi per atom and axis.
    Known answer: the velocity autocorrelation is cos(omega t) times a constant, so
    the vibrational density of states has a single peak at omega (frequency omega/2pi).
    Returns (positions, velocities), each with shape (n_frames, n_atoms, 3).
    """
    t = np.arange(n_frames)[:, None, None] * dt
    phase = rng.uniform(0.0, 2.0 * np.pi, size=(1, n_atoms, 3))
    positions = amplitude * np.cos(omega * t + phase)
    velocities = -amplitude * omega * np.sin(omega * t + phase)
    return positions, velocities


def dipole_series(n_frames, mean, fluctuation, rng):
    """Total-dipole time series M(t) with a set fluctuation.

    Each component is Gaussian with variance fluctuation/3 around mean, so the
    known answer is <|M|^2> - |<M>|^2 = fluctuation, the quantity the
    fluctuation formula for the dielectric constant is built on.
    Returns an array with shape (n_frames, 3).
    """
    mean = np.broadcast_to(np.asarray(mean, dtype=float), (3,))
    return mean + rng.normal(0.0, np.sqrt(fluctuation / 3.0), size=(n_frames, 3))


def lattice_from_data(path):
    """Read a perfect crystal from a LAMMPS data file (e.g. starting-structures/).

    Returns (positions, types, cell lengths). Known answer: the RDF peaks sit at
    the crystal's neighbor-shell distances, which follow from the cell.
    """
    from ase.io import read

    atoms = read(str(path), format="lammps-data")
    return atoms.get_positions(), atoms.arrays["type"].copy(), atoms.cell.lengths()
