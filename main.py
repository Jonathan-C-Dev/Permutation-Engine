"""
Command-line demo.

    python main.py verify [--deep]     recompute group orders / Lagrange indices
    python main.py solve3 [--n 5]      random-state 3x3x3 solves (Thistlethwaite)
    python main.py solve2 [--n 5]      random-state 2x2x2 optimal solves
    python main.py group  "(0 1 2)" "(0 1)"   order of a group from cycle notation
"""
from __future__ import annotations

import argparse
import random
import re
import time

from cosets import ColoringSpace
from permgroup import StabChain, from_cycles
from puzzles import Cube
from thistlethwaite import G_GENS, Thistlethwaite, fmt


def cmd_verify(args):
    t = time.time()
    Thistlethwaite().verify(deep=args.deep)
    c2 = Cube(2, "URF")
    o = StabChain(c2.N, c2.gens(["U", "R", "F"])).order()
    assert o == 3_674_160
    print(f"2x2x2 group order {o:,} OK")
    print(f"all checks passed in {time.time() - t:.1f}s")


def cmd_solve3(args):
    rng = random.Random(args.seed)
    T = Thistlethwaite()
    chain = StabChain(T.cube.N, T.cube.gens(G_GENS[0]))
    t0 = time.time()
    T.table3                                     # build/load phase-3 table
    print(f"phase-3 table ready ({time.time() - t0:.1f}s)\n")
    total = 0
    for i in range(args.n):
        state = chain.random_element(rng)        # uniformly random cube state
        t = time.time()
        phases = T.solve(state)
        dt = time.time() - t
        sol = sum(phases, [])
        assert T.cube.apply(state, sol) == T.cube.id
        total += len(sol)
        print(f"#{i + 1}: {len(sol)} moves in {dt * 1000:.0f} ms  phases {[len(p) for p in phases]}")
        print("   " + "  |  ".join(fmt(p) for p in phases))
    print(f"\naverage {total / args.n:.1f} moves")


def cmd_solve2(args):
    rng = random.Random(args.seed)
    c = Cube(2, "URF")
    names = [n for n in c.moves]
    space = ColoringSpace(c.label_identity(), [c.moves[n] for n in names])
    chain = StabChain(c.N, c.gens(["U", "R", "F"]))
    for i in range(args.n):
        state = chain.random_element(rng)
        t = time.time()
        path = space.solve(state)
        seq = [names[j] for j in path]
        assert c.apply(state, seq) == c.id
        print(f"#{i + 1}: {len(seq)} moves (optimal) in {(time.time() - t) * 1000:.0f} ms: {fmt(seq)}")


def cmd_group(args):
    cyc = [[int(x) for x in re.findall(r"\d+", s)] for s in args.cycles]
    n = max(max(c) for c in cyc) + 1
    ch = StabChain(n, [from_cycles(n, [c]) for c in cyc])
    print(f"degree {n}, order {ch.order():,}, orbit sizes along chain {ch.orbit_sizes()}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("verify"); a.add_argument("--deep", action="store_true"); a.set_defaults(f=cmd_verify)
    a = sub.add_parser("solve3"); a.add_argument("--n", type=int, default=5); a.add_argument("--seed", type=int, default=0); a.set_defaults(f=cmd_solve3)
    a = sub.add_parser("solve2"); a.add_argument("--n", type=int, default=5); a.add_argument("--seed", type=int, default=0); a.set_defaults(f=cmd_solve2)
    a = sub.add_parser("group"); a.add_argument("cycles", nargs="+"); a.set_defaults(f=cmd_group)
    args = p.parse_args()
    args.f(args)
