"""The testkit's generators must have the answers they claim, or every test
built on them inherits the error. Each test checks one generator's stated
known answer, with a tolerance set by its sampling error."""

from pathlib import Path

import numpy as np
import pytest

from testkit import (
    dipole_series,
    harmonic_oscillators,
    ideal_gas,
    lattice_from_data,
    random_walk,
    read_dump,
    write_data,
    write_dump,
)

REPO = Path(__file__).resolve().parents[2]


def test_random_walk_msd_is_6Dt():
    rng = np.random.default_rng(1)
    n_atoms, n_frames, dt, D = 4000, 50, 0.01, 0.3
    pos = random_walk(n_atoms, n_frames, dt, D, rng)
    lag = np.arange(n_frames) * dt
    msd = ((pos - pos[0]) ** 2).sum(axis=2).mean(axis=1)
    # MSD at lag t averages n_atoms samples of a chi-square with 3 dof scaled by 2Dt:
    # relative standard error sqrt(2/(3 n_atoms)) = 1.3%; allow 5 standard errors.
    rel_se = np.sqrt(2.0 / (3.0 * n_atoms))
    assert np.allclose(msd[1:], 6 * D * lag[1:], rtol=5 * rel_se)
    assert np.array_equal(pos[0], np.zeros((n_atoms, 3)))


def test_ideal_gas_is_uniform():
    rng = np.random.default_rng(2)
    n, L = 200_000, 10.0
    pos = ideal_gas(n, L, rng)
    assert pos.min() >= 0.0 and pos.max() < L
    counts, _ = np.histogram(pos[:, 0], bins=20, range=(0, L))
    expected = n / 20
    # Poisson counts: each bin within 5 sigma of n/20.
    assert np.all(np.abs(counts - expected) < 5 * np.sqrt(expected))


def test_harmonic_oscillator_spectrum_peaks_at_omega():
    rng = np.random.default_rng(3)
    dt, omega, n_frames = 0.001, 2 * np.pi * 25.0, 4000
    pos, vel = harmonic_oscillators(50, n_frames, dt, omega, 0.1, rng)
    # Invariant: velocity is the time derivative of position. A central difference
    # of A cos(wt) is off by at most A w (w dt)^2 / 6 (its truncation error).
    truncation = 0.1 * omega * (omega * dt) ** 2 / 6
    assert np.allclose(np.gradient(pos, dt, axis=0)[1:-1], vel[1:-1], rtol=0, atol=1.01 * truncation)
    spectrum = (np.abs(np.fft.rfft(vel, axis=0)) ** 2).sum(axis=(1, 2))
    freqs = np.fft.rfftfreq(n_frames, dt)
    # 25 Hz falls exactly on a bin (frequency resolution 1/(n dt) = 0.25 Hz).
    assert freqs[np.argmax(spectrum)] == pytest.approx(25.0, abs=0.25)


def test_dipole_series_fluctuation():
    rng = np.random.default_rng(4)
    n, mean, fluct = 200_000, (1.0, -2.0, 0.5), 30.0
    m = dipole_series(n, mean, fluct, rng)
    measured = (m**2).sum(axis=1).mean() - (m.mean(axis=0) ** 2).sum()
    # Variance of a sum of 3 Gaussian variances: relative SE sqrt(2/(3n)); allow 5 SE.
    assert measured == pytest.approx(fluct, rel=5 * np.sqrt(2 / (3 * n)))
    assert m.mean(axis=0) == pytest.approx(mean, abs=5 * np.sqrt(fluct / 3 / n))


def test_dump_round_trip(tmp_path):
    rng = np.random.default_rng(5)
    frames = [
        {"id": np.arange(1, 6), "element": ["O", "H", "H", "O", "H"],
         "x": rng.uniform(-3, 13, 5), "y": rng.uniform(0, 10, 5), "z": rng.uniform(0, 10, 5)}
        for _ in range(3)
    ]
    path = tmp_path / "t.lammpstrj"
    write_dump(path, frames, ((0, 0, 0), (10, 10, 10)), ("id", "element", "x", "y", "z"),
               timesteps=[0, 10, 20], wrap=True)
    back = read_dump(path)
    assert [s for s, _, _ in back] == [0, 10, 20]
    for (_, box, data), frame in zip(back, frames):
        assert np.allclose(box, [[0, 10]] * 3)
        assert list(data["element"]) == frame["element"]
        assert np.allclose(data["x"], np.mod(frame["x"], 10.0), atol=1e-6)
        assert data["x"].min() >= 0 and data["x"].max() < 10


def test_data_file_read_by_ase(tmp_path):
    """Independent implementation: ASE must read back what write_data wrote."""
    from ase.io import read

    pos = np.array([[0.5, 0.5, 0.5], [1.5, 2.0, 2.5], [3.0, 1.0, 0.25]])
    path = tmp_path / "t.data"
    write_data(path, pos, [1, 2, 2], {1: 15.9994, 2: 1.008}, ((0, 0, 0), (4, 4, 4)))
    atoms = read(str(path), format="lammps-data", atom_style="atomic")
    assert np.allclose(atoms.get_positions(), pos)
    assert np.allclose(atoms.cell.lengths(), 4.0)
    assert list(atoms.arrays["type"]) == [1, 2, 2]


def test_alpha_cristobalite_lattice():
    pos, types, lengths = lattice_from_data(REPO / "starting-structures/Alpha-cristobalite.data")
    assert len(pos) == 12 and (types == 1).sum() == 4 and (types == 2).sum() == 8
    si, o = pos[types == 1], pos[types == 2]
    d = si[:, None, :] - o[None, :, :]
    d -= lengths * np.round(d / lengths)
    nearest = np.linalg.norm(d, axis=2).min()
    # Si-O bond in alpha-cristobalite: 1.60-1.61 A (Downs & Palmer, Am. Mineral. 79, 9 (1994)).
    assert 1.59 < nearest < 1.63
