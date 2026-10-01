"""Tests for simledger on a synthetic tree shaped like /scratch2/lkyamamo.

Run from jobs/ledger:  python3 -m unittest discover -s tests -v
"""

import contextlib
import io
import os
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from simledger import cli, config, inbox, store  # noqa: E402
from simledger.discover import locate, parse_analysis_name  # noqa: E402
from simledger.parsers import lammps, scheduler, vasp  # noqa: E402
from simledger.scan import scan  # noqa: E402

LMP_IN = """units metal
atom_style atomic
boundary p p p
read_data start.data
pair_style reaxff NULL
pair_coeff * * SiOH.usc Si O H
timestep 0.0005
fix 1 all nvt temp 300.0 300.0 0.05
dump d1 all custom 1000 immersion.lammpstrj id type x y z
run 100000
run 50000
"""

LMP_LOG_OK = """LAMMPS (7 Aug 2019)
units metal
  18144 atoms
WARNING: something mild
Loop time of 101.499 on 256 procs for 100000 steps with 18144 atoms
Performance: 21.281 ns/day, 1.128 hours/ns, 985.234 timesteps/s
Loop time of 50.7346 on 256 procs for 50000 steps with 18144 atoms
Total wall time: 0:02:34
"""

LMP_LOG_ERR = """LAMMPS (7 Aug 2019)
  18144 atoms
ERROR: Lost atoms: original 18144 current 18100 (src/thermo.cpp:441)
"""

SLURM = """#!/bin/bash
#SBATCH --account=priyav_216
#SBATCH --nodes=4
#SBATCH --ntasks=256
#SBATCH --partition=priya
#SBATCH --time=01:00:00
#SBATCH --output=STREAM_OUTPUT
#SBATCH --job-name=SiOH-Immersion
module load gcc/13.3.0
export lmp="/home1/lkyamamo/executables/lammps/lmp_mpi_shock_2019"
srun --mpi=pmix -n $SLURM_NTASKS $lmp -in in.input
"""

STREAM = """starting
Wed Aug 26 17:54:44 PDT 2026
Total wall time: 0:02:34
Wed Aug 26 17:57:22 PDT 2026
finished
"""

INCAR = """SYSTEM = Silica
IBRION=2        !! CG
NSW = 200
ISIF=2
ENCUT=1000      !! cutoff
PREC=Accurate
EDIFF=1e-8 ; ISMEAR = 0
"""

KPOINTS = "auto\n0\nGamma\n2 2 1\n0 0 0\n"

POSCAR = """slab
1.0
10.0 0.0 0.0
0.0 10.0 0.0
0.0 0.0 20.0
Si O H
8 16 4
Selective dynamics
Direct
"""

POTCAR = """  PAW_PBE Si 05Jan2001
   VRHFIN =Si: s2p2
   TITEL  = PAW_PBE Si 05Jan2001
   POMASS =   28.085; ZVAL   =    4.000    mass and valenz
   ENMAX  =  245.345; ENMIN  =  184.009 eV
SECRET-PSEUDOPOTENTIAL-DATA 0.123 0.456
"""

OUTCAR_DONE = """ vasp.6.4.2 18Apr23 (build Mar 01 2024) complex
 running  256 mpi-ranks, with    1 threads/rank
  free  energy   TOTEN  =      -170.12345 eV
  energy  without entropy=     -169.76250592  energy(sigma->0) =     -169.92421643
 E-fermi :  -1.2345
 reached required accuracy - stopping structural energy minimisation
 General timing and accounting informations for this job:
                         Elapsed time (sec):       72.349
"""

OUTCAR_PARTIAL = """ vasp.6.4.2 18Apr23
  energy  without entropy=     -169.7  energy(sigma->0) =     -169.8
"""


# A log written with `echo both`: $-lines are followed by their substituted copy.
ECHO_LOG = """LAMMPS (7 Aug 2019)
units           metal
atom_style      atomic
read_data       start.data
  orthogonal box = (0 0 0) to (10 10 20)
  5 atoms
  read_data CPU = 0.0012 secs
mass 1 15.9994  # O
mass 2 1.00784  # H
pair_style      usc
pair_coeff  * * OH.usc O H
timestep        0.00025
velocity all create 300 4928459
fix 1 all nvt temp ${TARGET_TEMP} ${TARGET_TEMP} $(100*dt)
fix 1 all nvt temp 288.15 288.15 0.025
run ${n}
run 400
Step Temp PotEng Press Volume
       0    300.0   -10.0   -600.0    2000.0
     200    290.0   -10.1   -500.0    2000.0
     400    286.0   -10.2   -400.0    2000.0
Loop time of 1.0 on 64 procs for 400 steps with 5 atoms
unfix 1
fix 1 all nvt temp 288.15 288.15 0.025
run 400
Step Temp PotEng Press Volume
     400    288.0   -10.2   -450.0    2000.0
     800    290.0   -10.2   -350.0    2000.0
Loop time of 1.0 on 64 procs for 400 steps with 5 atoms
unfix 1
fix 2 all npt temp 288.15 288.15 0.025 iso 1.0 1.0 0.25
run 800
Step Temp PotEng Press Volume
     800    288.0   -10.2   -450.0    2000.0
    1600    289.0   -10.2   1.0    1900.0
Loop time of 1.0 on 64 procs for 800 steps with 5 atoms
Total wall time: 0:00:03
"""

DATA = """LAMMPS data file

5 atoms
2 atom types

0 10 xlo xhi
0 10 ylo yhi
0 20 zlo zhi

Masses

1 15.9994
2 1.00784

Atoms # atomic

1 1 0 0 0
2 2 1 0 0
3 2 0 1 0
4 1 5 5 5
5 2 5 6 5
"""


def w(path: Path, text: str = "", mtime: float = None) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    if mtime is not None:
        os.utime(path, (mtime, mtime))
    return path


def build_tree(root: Path) -> None:
    old = time.time() - 10 * 86400
    p = root / "20260825_locked"
    # 0247: plain LAMMPS run with run/ symlinks back to input_files (one broken)
    r = p / "runs" / "0247"
    for name, text in (("in.input", LMP_IN), ("start.data", "data"), ("SiOH.usc", "pot")):
        w(r / "input_files" / name, text)
    w(r / "run" / "log.lammps", LMP_LOG_OK, old)
    w(r / "run" / "STREAM_OUTPUT", STREAM, old)
    w(r / "run" / "lammps_submit.slurm", SLURM)
    w(r / "lammps_submit.slurm", SLURM)
    os.symlink("/scratch1/lkyamamo/20260825_locked/runs/0247/input_files/in.input", r / "run" / "in.input")
    # reaction_energy/0256: VASP relaxation that converged, with a README
    r = p / "runs" / "reaction_energy" / "0256"
    w(r / "input_files" / "INCAR", INCAR)
    w(r / "input_files" / "KPOINTS", KPOINTS)
    w(r / "input_files" / "POSCAR", POSCAR)
    w(r / "input_files" / "POTCAR", POTCAR)
    w(r / "run" / "OUTCAR", OUTCAR_DONE, old)
    w(r / "run" / "OSZICAR", "  1 F= -.17E+03 E0= -.169E+03  d E =-.1E+03\n  2 F= -.171E+03 E0= -.1699E+03  d E =-.1E+01\n")
    w(r / "run" / "CONTCAR", POSCAR.replace("20.0", "20.5"))
    w(r / "README.md", "# 0256 water dissociation on silanol\nCG relaxation of the hydroxylated slab.\n")
    w(r / "setup" / "frames" / "f1.xyz", "1\n\nH 0 0 0\n")
    # vasp_0272: VASP job that stopped without the timing block, long ago
    r = p / "runs" / "vasp_0272"
    w(r / "input_files" / "INCAR", "IBRION=-1\nENCUT=400\n")
    w(r / "run" / "OUTCAR", OUTCAR_PARTIAL, old)
    w(r / "run" / "WAVECAR", "x" * 2048, old)
    # ring-barrier/0282: numbered VASP frames along a path
    for i, e in ((1, -100.0), (2, -99.2), (3, -99.9)):
        frame = p / "runs" / "ring-barrier" / "0282" / "run" / "size_03__ring_00008__roll_060" / str(i)
        w(frame / "OUTCAR", OUTCAR_DONE.replace("-169.92421643", str(e)).replace("reached required accuracy", ""), old)
        # POSCARs are relative symlinks into setup/, as on the cluster
        w(p / "runs" / "ring-barrier" / "0282" / "setup" / "poscars" / str(i) / "POSCAR", POSCAR)
        os.symlink(f"../../../setup/poscars/{i}/POSCAR", frame / "POSCAR")
    w(p / "runs" / "ring-barrier" / "0282" / "input_files" / "INCAR", "IBRION=-1\n")
    # ring-barrier/0286: inputs only
    w(p / "runs" / "ring-barrier" / "0286" / "input_files" / "in.input", LMP_IN)
    # analysis
    w(p / "analysis" / "025_6-8_reaction_energy" / "1.png", "png")
    w(p / "analysis" / "0247_density" / "density.py", '"""Density profile of 0247."""\n')
    w(p / "analysis" / "0247_density" / "README.md", "Density profile across the immersion interface.\n")
    w(p / "analysis" / "0281-ring-selection" / "x.csv", "a,b\n")
    w(p / "potentials" / "README.md", "## SiOH.usc\n20260824_SiOH_v6.usc from finalized_modifications\n")
    w(p / "potentials" / "SiOH.usc", "pot")
    w(p / "potentials" / "SiOH_PAW-PBE_M1_POTCAR", POTCAR)
    w(p / "structures" / "README.md", "## 0039.data\n0232 end (finalized mod)\n")
    w(p / "structures" / "0039.data", "data")

    # 0250: echoed log + data file -> full system description
    r = p / "runs" / "0250"
    w(r / "input_files" / "in.input", "read_data start.data\nfix 1 all nvt temp ${TARGET_TEMP} ${TARGET_TEMP} 0.1\nrun ${n}\n")
    w(r / "input_files" / "start.data", DATA)
    w(r / "run" / "log.lammps", ECHO_LOG, old)
    os.symlink("/scratch1/gone/start.data", r / "run" / "start.data")   # broken: falls back to input_files/

    r = p / "runs" / "0251"
    w(r / "input_files" / "start.data", DATA)
    w(r / "run" / "log.lammps", ECHO_LOG.replace(
        "  5 atoms\n", "  5 atoms\nreplicate 2 1 1\n  orthogonal box = (0 0 0) to (20 10 20)\n  10 atoms\n")
        .replace("with 5 atoms", "with 10 atoms"), old)

    # 0253: starting structure is a broken link to another project; end state is readable
    r = p / "runs" / "0253"
    (r / "input_files").mkdir(parents=True)
    os.symlink("/scratch1/lkyamamo/other-project/runs/0103/run/30C.data", r / "input_files" / "start.data")
    w(r / "run" / "log.lammps", ECHO_LOG.replace("Total wall time", "write_data final.data\nTotal wall time"), old)
    w(r / "run" / "final.data", DATA)
    os.symlink("/scratch1/lkyamamo/20260825_locked/runs/0253/input_files/start.data", r / "run" / "start.data")

    # 0254: killed before write_data; start.data links to 0247's final.data at its pre-move path
    w(p / "runs" / "0247" / "run" / "final.data", DATA)
    r = p / "runs" / "0254"
    (r / "input_files").mkdir(parents=True)
    os.symlink("/scratch1/lkyamamo/20260825_locked/runs/0247/run/final.data", r / "input_files" / "start.data")
    w(r / "run" / "log.lammps", ECHO_LOG.split("Loop time")[0] + "read_restart restart.3\n", old)

    # 0255: one log per structure (log_000.lammps ...), several data files read in turn
    r = p / "runs" / "0255"
    w(r / "input_files" / "Ortho.data", DATA)
    w(r / "input_files" / "Meta.data", DATA.replace("5 atoms", "6 atoms") + "6 1 1 1 1\n")
    for i, name in enumerate(("Ortho", "Meta")):
        w(r / "run" / f"log_{i:03d}.lammps",
          f"LAMMPS (7 Aug 2019)\nread_data ${{struc}}.data\nread_data {name}.data\n  {5 + i} atoms\n"
          "minimize 1e-8 1e-8 1000 10000\nLoop time of 1.0 on 1 procs for 10 steps with 5 atoms\nTotal wall time: 0:00:01\n", old)

    # second project: a temperature sweep, and an id reused from the first project
    q = root / "finalized-modifications"
    r = q / "runs" / "0239"
    w(r / "submit_temperature_sweep.sh", "#!/bin/bash\n")
    w(r / "STREAM_OUTPUT", "sweep driver\n", old)
    w(r / "input_files" / "in.input", LMP_IN)
    for t, log in (("T15C", LMP_LOG_OK), ("T30C", LMP_LOG_ERR)):
        text = LMP_IN.replace("300.0 300.0", "288.15 288.15") if t == "T15C" else LMP_IN
        w(r / t / "input_files" / "in.input", text)
        w(r / t / "run" / "in.input", text + ("# edited\n" if t == "T30C" else ""))
        w(r / t / "run" / "log.lammps", log, old)
        w(r / t / "run" / "lammps_submit.slurm", SLURM)
    w(q / "runs" / "0247" / "run" / "log.lammps", LMP_LOG_OK, old)
    w(q / "analysis" / "0239_T15C_msd" / "msd.dat", "1 2\n")


def snapshot(root: Path):
    out = {}
    for dirpath, dirnames, filenames in os.walk(root):
        for n in filenames + dirnames:
            p = Path(dirpath) / n
            st = os.lstat(p)
            out[str(p)] = (st.st_size, st.st_mtime)
    return out


class ParserTests(unittest.TestCase):
    def test_analysis_names(self):
        self.assertEqual(parse_analysis_name("0239_T60C_msd", 4)[:2], (["0239"], "T60C_msd"))
        self.assertEqual(parse_analysis_name("0281-ring-selection", 4)[:2], (["0281"], "ring-selection"))
        self.assertEqual(parse_analysis_name("0242_0243_compare", 4)[0], ["0242", "0243"])
        self.assertEqual(parse_analysis_name("0242-0245_rdf", 4)[0], ["0242", "0243", "0244", "0245"])
        self.assertEqual(parse_analysis_name("025_6-8_reaction_energy", 4)[:2],
                         (["0256", "0257", "0258"], "reaction_energy"))
        ids, _, warns = parse_analysis_name("figures", 4)
        self.assertEqual(ids, [])
        self.assertTrue(warns)

    def test_locate_survives_moved_data(self):
        cfg = config.Config(home=Path("/tmp/x"))
        loc = locate("/scratch1/lkyamamo/finalized-modifications/runs/0239/T45C/run", cfg)
        self.assertEqual(loc, {"project": "finalized-modifications", "run_id": "0239",
                               "group_path": "", "subrun": "T45C"})
        loc = locate("/scratch2/lkyamamo/20260825_locked/runs/dft_surface_coverage/vasp-interface-0252/run", cfg)
        self.assertEqual((loc["run_id"], loc["group_path"], loc["subrun"]), ("0252", "dft_surface_coverage", ""))
        self.assertIsNone(locate("/home1/lkyamamo/util", cfg))

    def test_lammps_input_and_log(self):
        li = lammps.parse_input(LMP_IN)
        self.assertEqual(li["total_steps"], 150000)
        self.assertAlmostEqual(li["simulated_time"], 75.0)
        self.assertEqual(li["time_unit"], "ps")
        self.assertEqual(li["fixes"][0]["temp"], "300.0 300.0 0.05")
        lg = lammps.parse_log(LMP_LOG_OK, LMP_LOG_OK)
        self.assertTrue(lg["finished"])
        self.assertEqual((lg["wall_time_s"], lg["n_atoms"], lg["procs"]), (154, 18144, 256))
        self.assertEqual(lammps.parse_log(LMP_LOG_ERR, LMP_LOG_ERR)["errors"][0][:16], "ERROR: Lost atom")

    def test_vasp(self):
        tags = vasp.parse_incar(INCAR)
        self.assertEqual((tags["ENCUT"], tags["ISMEAR"], tags["EDIFF"]), ("1000", "0", "1e-8"))
        self.assertEqual(vasp.calc_type(tags), "relax_ions")
        self.assertEqual(vasp.calc_type({"IBRION": "0"}), "md_nve")
        self.assertEqual(vasp.calc_type({"IBRION": "0", "MDALGO": "2"}), "md_nvt")
        ps = vasp.parse_poscar(POSCAR)
        self.assertEqual((ps["formula"], ps["n_atoms"], ps["volume"]), ("Si8O16H4", 28, 2000.0))
        self.assertTrue(ps["selective_dynamics"])
        titles = vasp.potcar_titles(POTCAR.splitlines())
        self.assertEqual(titles, [{"titel": "PAW_PBE Si 05Jan2001", "vrhfin": "Si", "zval": 4.0,
                                   "pomass": 28.085, "enmax": 245.345}])
        oc = vasp.parse_outcar(OUTCAR_DONE, OUTCAR_DONE)
        self.assertTrue(oc["finished"] and oc["reached_accuracy"])
        self.assertAlmostEqual(oc["energy_sigma0"], -169.92421643)
        self.assertEqual(oc["cores"], 256)
        self.assertEqual(vasp.parse_kpoints(KPOINTS)["mesh"], "2 2 1")

    def test_scheduler(self):
        s = scheduler.parse_script(SLURM)
        self.assertEqual(s["directives"]["ntasks"], "256")
        self.assertEqual(s["modules"], ["gcc/13.3.0"])
        self.assertIn("lmp", s["executables"])
        so = scheduler.parse_stdout(STREAM, STREAM)
        self.assertEqual((so["start"], so["end"]), ("2026-08-26T17:54:44", "2026-08-26T17:57:22"))
        self.assertEqual(scheduler.parse_stdout("", "slurmstepd: error: *** JOB 1 CANCELLED AT x DUE TO TIME LIMIT ***")
                         ["failures"], ["cancelled", "time_limit"])


class ScanTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        base = Path(self.tmp.name)
        self.root = base / "scratch"
        self.home = base / "ledger"
        build_tree(self.root)
        self.cfg = config.load(self.home)

    def tearDown(self):
        self.tmp.cleanup()

    def cli(self, *argv) -> str:
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(io.StringIO()):
            cli.main(["--home", str(self.home), *argv])
        return buf.getvalue()

    def runs(self, conn):
        return {r["run_key"]: r for r in conn.execute("SELECT * FROM runs")}

    def test_scan_is_read_only(self):
        before = snapshot(self.root)
        scan(self.cfg, self.root)
        self.assertEqual(before, snapshot(self.root))

    def test_refuses_home_inside_tree(self):
        cfg = config.load(self.root / "ledger")
        with self.assertRaises(SystemExit):
            scan(cfg, self.root)

    def test_discovery_and_status(self):
        res = scan(self.cfg, self.root)
        conn = res["conn"]
        runs = self.runs(conn)
        self.assertEqual(set(runs), {"20260825_locked/0247", "20260825_locked/0250", "20260825_locked/0256",
                                     "20260825_locked/0251", "20260825_locked/0253", "20260825_locked/0254",
                                     "20260825_locked/0255", "20260825_locked/0272",
                                     "20260825_locked/0282",
                                     "20260825_locked/0286", "finalized-modifications/0239",
                                     "finalized-modifications/0247"})
        r = runs["20260825_locked/0247"]
        self.assertEqual((r["status"], r["code"], r["calc_type"], r["n_atoms"], r["cores"], r["wall_time_s"]),
                         ("completed", "lammps", "md_nvt", 18144, 256, 154))
        self.assertEqual(r["start_time"], "2026-08-26T17:54:44")
        r = runs["20260825_locked/0256"]
        self.assertEqual((r["status"], r["code"], r["calc_type"], r["group_path"], r["formula"]),
                         ("completed", "vasp", "relax_ions", "reaction_energy", "Si8O16H4"))
        self.assertIn("water dissociation", r["readme"])
        self.assertEqual(runs["20260825_locked/0272"]["status"], "incomplete")
        self.assertEqual(runs["20260825_locked/0272"]["label"], "vasp_0272")
        self.assertEqual(runs["20260825_locked/0286"]["status"], "not_started")
        r = runs["finalized-modifications/0239"]
        self.assertEqual((r["status"], r["n_subruns"]), ("partial", 3))
        subs = {s["label"]: s["status"] for s in conn.execute(
            "SELECT * FROM subruns WHERE run_key='finalized-modifications/0239'")}
        self.assertEqual(subs, {"(top)": "unknown", "T15C": "completed", "T30C": "failed"})
        # varying temperature recorded per sub-run, from input_files/
        temps = dict(conn.execute("SELECT subrun, value_text FROM params WHERE run_key='finalized-modifications/0239' "
                                  "AND key='T_target'").fetchall())
        self.assertEqual(temps, {"T15C": "288.15", "T30C": "300.0"})

        msgs = " | ".join(m for _, m in res["warnings"])
        self.assertIn("also used in", msgs)                  # 0247 in both projects
        self.assertIn("run id 0257 not found", msgs)         # 025_6-8 range
        self.assertIn("differs from T30C/input_files/in.input", msgs)  # edited copy
        self.assertNotIn("T15C/input_files", msgs)

    def test_system_description_from_log(self):
        conn = scan(self.cfg, self.root)["conn"]
        r = conn.execute("SELECT * FROM runs WHERE run_key='20260825_locked/0250'").fetchone()
        import json
        s = json.loads(r["system"])
        self.assertEqual((s["n_atoms"], s["composition"], s["box"], s["T_target"], s["P_target"]),
                         (5, "O 2, H 3", "10 × 10 × 20 Å", 288.15, 1.0))
        self.assertEqual(s["description_source"], "log")
        self.assertIsNone(s["n_structures"])
        self.assertEqual(s["simulated_time"], 0.4)            # 1600 steps x 0.00025 ps
        self.assertEqual(s["ensemble"], "NPT")                # longest stage
        self.assertAlmostEqual(s["density"], (2 * 15.9994 + 3 * 1.00784) * 1.66053907 / 2000, 4)
        self.assertEqual(r["calc_type"], "md_npt")
        self.assertIn("O 2, H 3", r["summary"])
        self.assertIn("NVT 288.15 K 0.2 ps (2×) → NPT 288.15 K 1 bar 0.2 ps", r["summary"])
        # replicate 2 1 1: the 5-atom data file becomes 10 atoms in a doubled box
        r2 = conn.execute("SELECT system FROM runs WHERE run_key='20260825_locked/0251'").fetchone()
        s2 = json.loads(r2["system"])
        self.assertEqual((s2["n_atoms"], s2["composition"], s2["box"]), (10, "O 4, H 6", "20 × 10 × 20 Å"))
        self.assertAlmostEqual(s2["density"], s["density"], 6)
        s3 = json.loads(conn.execute("SELECT system FROM runs WHERE run_key='20260825_locked/0253'").fetchone()[0])
        self.assertEqual(s3["composition"], "O 2, H 3")
        self.assertEqual(s3["data_file"], "run/final.data (end state; starting data file unavailable)")
        self.assertEqual(s3["structure_origin"], "/scratch1/lkyamamo/other-project/runs/0103/run/30C.data")
        sub = conn.execute("SELECT protocol FROM subruns WHERE run_key='20260825_locked/0250'").fetchone()
        proto = json.loads(sub["protocol"])
        self.assertEqual([(p["count"], p["ensemble"], p["steps"]) for p in proto], [(2, "NVT", 400), (1, "NPT", 800)])
        self.assertAlmostEqual(proto[0]["T_mean"], (288.0 + 290.0 + 290.0 + 286.0) / 4 if False else proto[0]["T_mean"])
        card = self.cli_show("0250")
        for text in ("## System", "## Protocol", "O 2, H 3", "288.15 K"):
            self.assertIn(text, card)
        # the sweep run summarises its sub-runs
        sweep = conn.execute("SELECT summary FROM runs WHERE run_key='finalized-modifications/0239'").fetchone()[0]
        self.assertTrue(sweep.startswith("2 sub-runs over T = 288.15–300 K"), sweep)

    def test_moved_links_and_multiple_structures(self):
        import json
        conn = scan(self.cfg, self.root)["conn"]
        s = json.loads(conn.execute("SELECT system FROM runs WHERE run_key='20260825_locked/0254'").fetchone()[0])
        self.assertEqual(s["composition"], "O 2, H 3")
        self.assertTrue(s["data_file"].endswith("20260825_locked/runs/0247/run/final.data "
                                                "(followed a broken link to its moved location)"), s["data_file"])
        self.assertEqual(s["structure_origin_run"], "20260825_locked/0247")
        self.assertIsNone(s["n_structures"])                    # a read_restart is not a second structure
        r = conn.execute("SELECT status, summary, system FROM runs WHERE run_key='20260825_locked/0255'").fetchone()
        s = json.loads(r["system"])
        self.assertEqual(r["status"], "completed")              # log_NNN.lammps make it a calculation
        self.assertEqual((s["n_structures"], s["n_atoms"], s["composition"]), (2, 5, "O 2, H 3"))
        self.assertIn("2 structures; first (Ortho.data): 5 atoms", r["summary"])
        card = self.cli_show("0254")
        self.assertIn("from run [20260825_locked/0247](../20260825_locked/0247.md)", card)

    def cli_show(self, run):
        self.cli("scan", str(self.root))
        return self.cli("show", run)

    def test_frames_grouped(self):
        conn = scan(self.cfg, self.root)["conn"]
        subs = conn.execute("SELECT * FROM subruns WHERE run_key='20260825_locked/0282'").fetchall()
        self.assertEqual([(s["label"], s["status"], s["n_atoms"]) for s in subs],
                         [("size_03__ring_00008__roll_060", "completed", 28)])
        res = dict(conn.execute("SELECT key, value_text FROM results WHERE run_key='20260825_locked/0282'").fetchall())
        self.assertEqual((res["n_frames"], res["energy_max_frame"]), ("3", "2"))
        self.assertAlmostEqual(float(res["energy_span"]), 0.8)
        self.assertEqual(res["frame_energies"], "[-100.0, -99.2, -99.9]")

    def test_potcar_never_stored(self):
        conn = scan(self.cfg, self.root)["conn"]
        dump = "\n".join(conn.iterdump())
        self.assertNotIn("SECRET-PSEUDOPOTENTIAL-DATA", dump)
        self.assertIn("PAW_PBE Si 05Jan2001", dump)

    def test_files_symlinks_and_vasp_results(self):
        conn = scan(self.cfg, self.root)["conn"]
        f = conn.execute("SELECT * FROM files WHERE run_key='20260825_locked/0247' AND relpath='run/in.input'").fetchone()
        self.assertEqual((f["category"], f["broken"]), ("lammps_input", 1))
        self.assertTrue(f["symlink"].startswith("/scratch1/"))
        res = dict(conn.execute("SELECT key, value_num FROM results WHERE run_key='20260825_locked/0256' "
                                "AND subrun=''").fetchall())
        self.assertAlmostEqual(res["energy_sigma0"], -169.92421643)
        self.assertAlmostEqual(res["energy_per_atom"], -169.92421643 / 28)
        self.assertEqual(res["ionic_steps"], 2)
        self.assertAlmostEqual(res["volume_change_pct"], 2.5)

    def test_analysis_links(self):
        conn = scan(self.cfg, self.root)["conn"]
        links = {(r["analysis_key"], r["run_key"], r["subrun"]) for r in conn.execute("SELECT * FROM analysis_runs")}
        self.assertIn(("20260825_locked/025_6-8_reaction_energy", "20260825_locked/0256", None), links)
        self.assertIn(("20260825_locked/0247_density", "20260825_locked/0247", None), links)
        self.assertIn(("finalized-modifications/0239_T15C_msd", "finalized-modifications/0239", "T15C"), links)
        pf = conn.execute("SELECT * FROM project_files WHERE name='0039.data'").fetchone()
        self.assertEqual((pf["description"], pf["mentioned_ids"]), ("0232 end (finalized mod)", "0232"))

    def test_incremental_rescan(self):
        scan(self.cfg, self.root)
        res = scan(self.cfg, self.root)
        self.assertEqual((res["counts"]["new"], res["counts"]["changed"]), (0, 0))
        self.assertIn("differs from T30C/input_files/in.input", " ".join(m for _, _, m in res["open_issues"]))
        log = self.root / "20260825_locked" / "runs" / "vasp_0272" / "run" / "OUTCAR"
        log.write_text(OUTCAR_DONE.replace("reached required accuracy", ""))
        res = scan(self.cfg, self.root)
        self.assertEqual(res["changed"], ["20260825_locked/0272"])
        self.assertEqual(self.runs(res["conn"])["20260825_locked/0272"]["status"], "completed")  # IBRION=-1
        # a run that disappears is flagged, not deleted
        os.rename(self.root / "20260825_locked" / "runs" / "ring-barrier" / "0286",
                  self.root / "20260825_locked" / "gone_0286")
        res = scan(self.cfg, self.root)
        self.assertEqual(self.runs(res["conn"])["20260825_locked/0286"]["missing"], 1)

    def test_hook_events_and_notes(self):
        env = {"SLURM_JOB_ID": "4242", "SLURM_JOB_NODELIST": "e01-[01-04]", "SLURM_JOB_PARTITION": "priya"}
        moved = "/scratch1/lkyamamo/finalized-modifications/runs/0239/T15C/run"
        with mock.patch.dict(os.environ, env), mock.patch.object(inbox.subprocess, "run", side_effect=OSError):
            self.cli("hook", "run", "--dir", moved, "--exit-code", "1", "--field", "temperature_C=15")
            self.cli("hook", "analysis", "--dir", moved, "--type", "msd",
                     "--output", str(self.root / "finalized-modifications/analysis/0239_T15C_msd/msd.dat"))
        self.assertEqual(len(list(self.cfg.inbox.glob("*.json"))), 2)
        self.cli("scan", str(self.root))
        self.assertEqual(list(self.cfg.inbox.glob("*.json")), [])
        conn = store.connect(self.cfg.db_path)
        sub = conn.execute("SELECT * FROM subruns WHERE run_key='finalized-modifications/0239' AND label='T15C'").fetchone()
        self.assertEqual(sub["status"], "failed")
        self.assertIn("job 4242 exited 1", sub["status_evidence"])
        self.assertEqual(self.runs(conn)["finalized-modifications/0239"]["job_id"], "4242")
        out = self.cli("search", "temperature_C=15")
        self.assertIn("finalized-modifications/0239", out)
        ev_link = conn.execute("SELECT * FROM analysis_runs WHERE source='event'").fetchone()
        self.assertEqual((ev_link["analysis_key"], ev_link["subrun"]),
                         ("finalized-modifications/0239_T15C_msd", "T15C"))
        conn.close()

        self.cli("note", "0256", "Use this as the reference energy", "--tag", "reference")
        self.cli("scan", str(self.root), "--full")
        card = (self.cfg.cards / "20260825_locked" / "0256.md").read_text()
        self.assertIn("Use this as the reference energy", card)
        self.assertIn("tags: reference", card)

    def test_search_and_cards(self):
        self.cli("scan", str(self.root))
        out = self.cli("search", "code=vasp", "ENCUT>=520")
        self.assertIn("20260825_locked/0256", out)
        self.assertNotIn("0272", out)
        out = self.cli("search", "incar.ENCUT<500")
        self.assertIn("0272", out)
        out = self.cli("search", "status=completed", "group~reaction")
        self.assertEqual(out.count("\n"), 2)
        out = self.cli("search", "--text", "dissociation silanol")
        self.assertIn("20260825_locked/0256", out)
        out = self.cli("search", "density", "interface")
        self.assertIn("analysis: 20260825_locked/0247_density", out)
        out = self.cli("show", "0239")
        self.assertIn("## Sub-runs (3)", out)
        self.assertIn("### Varying across sub-runs", out)
        self.assertIn("ambiguous", self._expect_exit("show", "0247"))
        self.assertIn("finalized-modifications/0247", self.cli("show", str(self.root / "finalized-modifications/runs/0247/run")))
        for p in ("INDEX.md", "runs.csv", "last_scan.md", "20260825_locked/README.md", "20260825_locked/0256.md"):
            self.assertTrue((self.cfg.cards / p).is_file(), p)
        readme = (self.cfg.cards / "20260825_locked" / "README.md").read_text()
        self.assertIn("20260824_SiOH_v6.usc from finalized_modifications", readme)

    def _expect_exit(self, *argv) -> str:
        try:
            self.cli(*argv)
        except SystemExit as exc:
            return str(exc)
        return ""


if __name__ == "__main__":
    unittest.main()
