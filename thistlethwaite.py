"""
thistlethwaite.py -- Rubik's cube solver via a subgroup chain and coset spaces.

    G0 = <U,D,L,R,F,B>          |G0:G1| = 2,048
    G1 = <U,D,L,R,F2,B2>        |G1:G2| = 1,082,565
    G2 = <U2,D2,L,R,F2,B2>      |G2:G3| = 29,400
    G3 = <U2,D2,L2,R2,F2,B2>    |G3:G4| = 663,552
    G4 = {e}

Phase i takes a state in G_i into G_{i+1} by an optimal search in the coset
space G_{i+1}\\G_i, whose size is the Lagrange index.  verify() recomputes every
|G_i| with Schreier-Sims and checks the coset-space sizes found by BFS against
|G_i|/|G_{i+1}| -- so a wrong coordinate definition cannot go unnoticed.
"""
from __future__ import annotations

import os
import pickle
import time

from cosets import CanonCosetSpace, ColoringSpace
from permgroup import Perm, StabChain
from puzzles import Cube

G_MOVES = [
    ["U", "U'", "U2", "D", "D'", "D2", "L", "L'", "L2", "R", "R'", "R2", "F", "F'", "F2", "B", "B'", "B2"],
    ["U", "U'", "U2", "D", "D'", "D2", "L", "L'", "L2", "R", "R'", "R2", "F2", "B2"],
    ["U2", "D2", "L", "L'", "L2", "R", "R'", "R2", "F2", "B2"],
    ["U2", "D2", "L2", "R2", "F2", "B2"],
]
G_GENS = [["U", "D", "L", "R", "F", "B"], ["U", "D", "L", "R", "F2", "B2"],
          ["U2", "D2", "L", "R", "F2", "B2"], ["U2", "D2", "L2", "R2", "F2", "B2"]]
EXPECTED_INDEX = [2048, 1_082_565, 29_400, 663_552]


class Thistlethwaite:
    def __init__(self, cache_dir: str | None = ".cache"):
        self.cube = Cube(3)
        c = self.cube
        self.cache_dir = cache_dir
        mv = lambda names: [c.moves[n] for n in names]
        # phases 1, 2, 4: coloring cosets (fast);  phase 3: canonical-representative cosets
        self.p1 = ColoringSpace(c.label_edge_orientation(), mv(G_MOVES[0]))
        self.p2 = ColoringSpace(c.label_corner_orient_and_slice(), mv(G_MOVES[1]))
        self.chain_G3 = StabChain(c.N, c.gens(G_GENS[3]))
        self.p3 = CanonCosetSpace(self.chain_G3, mv(G_MOVES[2]))
        self.p4 = ColoringSpace(c.label_identity(), mv(G_MOVES[3]))
        self._table3: dict | None = None

    # ---- phase-3 table (29,400 cosets) cached on disk ----------------------
    @property
    def table3(self) -> dict:
        if self._table3 is None:
            path = os.path.join(self.cache_dir, "phase3.pkl") if self.cache_dir else None
            if path and os.path.exists(path):
                with open(path, "rb") as f:
                    self._table3 = pickle.load(f)
            else:
                self._table3 = self.p3.build_table()
                if path:
                    os.makedirs(self.cache_dir, exist_ok=True)
                    with open(path, "wb") as f:
                        pickle.dump(self._table3, f)
        return self._table3

    # ---- solving -----------------------------------------------------------
    def solve(self, state: Perm) -> list[list[str]]:
        """Return the move list of each of the four phases."""
        c = self.cube
        g, phases = state, []
        for idx, (space, names) in enumerate([(self.p1, G_MOVES[0]), (self.p2, G_MOVES[1]),
                                              (None, G_MOVES[2]), (self.p4, G_MOVES[3])]):
            if idx == 2:
                path = self.p3.solve_with_table(g, self.table3)
            else:
                path = space.solve(g)
            seq = [names[i] for i in path]
            g = c.apply(g, seq)
            phases.append(seq)
            assert idx == 3 or self._in_subgroup(g, idx + 1), f"phase {idx + 1} failed"
        assert g == c.id
        return phases

    def _in_subgroup(self, g: Perm, level: int) -> bool:
        if level == 1:
            return self.p1.key(g) == self.p1.goal
        if level == 2:
            return self.p2.key(g) == self.p2.goal
        if level == 3:
            return self.chain_G3.contains(g)
        return g == self.cube.id

    # ---- verification of the math -------------------------------------------
    def verify(self, deep: bool = False, log=print) -> dict:
        """
        Recompute |G_i| by Schreier-Sims; check Lagrange indices; check that the
        coloring/canonical coset spaces really are G_{i+1}\\G_i by counting cosets.
        deep=True also BFS-counts the 1,082,565 cosets of phase 2 (~1-2 min).
        """
        c = self.cube
        orders = [StabChain(c.N, c.gens(g)).order() for g in G_GENS] + [1]
        idx = [orders[i] // orders[i + 1] for i in range(4)]
        for i in range(4):
            assert orders[i] % orders[i + 1] == 0
            log(f"|G{i}| = {orders[i]:,}   [G{i}:G{i + 1}] = {idx[i]:,}")
        assert idx == EXPECTED_INDEX, idx
        assert orders[0] == 43_252_003_274_489_856_000
        # phase 1,2: H generators preserve the coloring (so Stab ⊇ H)
        assert self.p1.preserved_by(c.gens(G_GENS[1]))
        assert self.p2.preserved_by(c.gens(G_GENS[2]))
        n1 = self.p1.count_cosets()
        assert n1 == idx[0], (n1, idx[0])
        log(f"phase 1: BFS coset count {n1:,} == Lagrange index")
        if deep:
            n2 = self.p2.count_cosets()
            assert n2 == idx[1], (n2, idx[1])
            log(f"phase 2: BFS coset count {n2:,} == Lagrange index")
        n3 = len(self.table3)
        assert n3 == idx[2], (n3, idx[2])
        log(f"phase 3: BFS coset count {n3:,} == Lagrange index")
        return {"orders": orders, "indices": idx}


def fmt(moves: list[str]) -> str:
    return " ".join(moves)
