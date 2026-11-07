"""
cosets.py -- searching in coset spaces H\\G instead of in G.

Why cosets?  Solving "get g into the subgroup H" only depends on the right coset
Hg, because  Hg·w = H  <=>  g·w in H.  Right multiplication by a move is a
well-defined action on right cosets, so we can run BFS on the (much smaller)
coset graph of size [G:H] = |G|/|H|  (Lagrange) instead of on G.

Two ways to name a coset:

1. ColoringSpace -- H is "the states that look solved under a sticker coloring".
   Then Hg <-> the coloring label∘g, and moves act by permuting the coloring
   (no group elements needed; keys are 48-byte strings).  Fast.

2. CanonCosetSpace -- H is any subgroup given by generators.  Uses the
   Schreier-Sims chain of H to compute a canonical representative of g^-1·H.
   Slower but completely general (needed for Thistlethwaite phase 3, whose
   subgroup is NOT a coloring stabilizer).
"""
from __future__ import annotations

from collections import deque
from operator import itemgetter
from typing import Callable, Sequence

from permgroup import Perm, StabChain, compose, identity, inverse


# ---------------------------------------------------------------------------
# bidirectional BFS on an abstract graph with invertible edge labels
# ---------------------------------------------------------------------------
def bidirectional_search(start, goal, fwd: Callable, bwd: Callable, n_moves: int,
                         max_depth: int = 64) -> list[int]:
    """
    Shortest path (as a list of move indices) from `start` to `goal`.
    fwd(k, i): node reached from k by move i.
    bwd(k, i): node from which move i leads to k  (i.e. the inverse action).
    """
    if start == goal:
        return []
    F = {start: (None, -1, 0)}          # node -> (parent, move, depth)
    B = {goal: (None, -1, 0)}
    ff, bf = [start], [goal]
    depth = 0
    while ff and bf:
        depth += 1
        if depth > max_depth:
            break
        forward = len(ff) <= len(bf)
        frontier, mine, other, step = (ff, F, B, fwd) if forward else (bf, B, F, bwd)
        new, best, best_node = [], None, None
        for k in frontier:
            d = mine[k][2] + 1
            for i in range(n_moves):
                nk = step(k, i)
                if nk in mine:
                    continue
                mine[nk] = (k, i, d)
                new.append(nk)
                if nk in other:
                    total = d + other[nk][2]
                    if best is None or total < best:
                        best, best_node = total, nk
        if best is not None:
            return _stitch(F, B, best_node)
        if forward:
            ff = new
        else:
            bf = new
    raise RuntimeError("no path found (state not in the expected subgroup?)")


def _stitch(F, B, meet) -> list[int]:
    left = []
    k = meet
    while F[k][0] is not None:
        left.append(F[k][1])
        k = F[k][0]
    left.reverse()
    right = []
    k = meet
    while B[k][0] is not None:          # k = bwd(parent, i)  =>  move i leads k -> parent
        right.append(B[k][1])
        k = B[k][0]
    return left + right


# ---------------------------------------------------------------------------
# 1. colorings
# ---------------------------------------------------------------------------
class ColoringSpace:
    """Right cosets Hg where H = {states that look solved under `label`}."""

    def __init__(self, label: Sequence[int], moves: Sequence[Perm]):
        assert max(label) < 256
        self.label = bytes(label)
        self.moves = list(moves)
        self._fwd = [itemgetter(*m) for m in self.moves]
        self._bwd = [itemgetter(*inverse(m)) for m in self.moves]
        self.goal = self.label                      # label∘identity

    def key(self, g: Perm) -> bytes:
        lab = self.label
        return bytes(lab[x] for x in g)             # label∘g

    def fwd(self, k: bytes, i: int) -> bytes:       # coloring of g∘m
        return bytes(self._fwd[i](k))

    def bwd(self, k: bytes, i: int) -> bytes:       # coloring of g∘m^-1
        return bytes(self._bwd[i](k))

    def solve(self, g: Perm) -> list[int]:
        """Move indices taking g into H."""
        return bidirectional_search(self.key(g), self.goal, self.fwd, self.bwd, len(self.moves))

    def count_cosets(self, limit: int = 5_000_000) -> int:
        """Size of the orbit of the solved coloring = [G_i : Stab]. (Orbit-stabilizer.)"""
        seen = {self.goal}
        q = deque([self.goal])
        while q:
            k = q.popleft()
            for i in range(len(self.moves)):
                nk = self.fwd(k, i)
                if nk not in seen:
                    seen.add(nk)
                    q.append(nk)
                    if len(seen) > limit:
                        raise RuntimeError("coset space larger than limit")
        return len(seen)

    def preserved_by(self, gens: Sequence[Perm]) -> bool:
        """Do these generators of H keep the coloring looking solved?  label∘h == label."""
        lab = self.label
        return all(lab[h[i]] == lab[i] for h in gens for i in range(len(h)))


# ---------------------------------------------------------------------------
# 2. canonical-representative cosets (any subgroup)
# ---------------------------------------------------------------------------
class CanonCosetSpace:
    """
    Right cosets Hg for an arbitrary subgroup H (given as a StabChain).
    key(Hg) = canon_left(g^-1): a canonical element of the left coset g^-1·H.
    A node is stored as the pair (g, g^-1).
    """

    def __init__(self, chain_H: StabChain, moves: Sequence[Perm]):
        self.H = chain_H
        self.moves = list(moves)
        self.minv = [inverse(m) for m in self.moves]
        self.n = chain_H.n

    def key(self, ginv: Perm) -> Perm:
        return self.H.canon_left(ginv)

    def step(self, g: Perm, ginv: Perm, i: int):
        """(g, g^-1)  ->  (g∘m_i, m_i^-1∘g^-1)"""
        return compose(g, self.moves[i]), compose(self.minv[i], ginv)

    def build_table(self) -> dict[Perm, int]:
        """BFS from the goal coset H, applying inverse moves: dist[key] = #moves to reach H."""
        e = identity(self.n)
        dist = {self.key(e): 0}
        frontier = [(e, e)]
        d = 0
        while frontier:
            d += 1
            nxt = []
            for g, ginv in frontier:
                for i in range(len(self.moves)):
                    g2 = compose(g, self.minv[i])          # g∘m^-1
                    ginv2 = compose(self.moves[i], ginv)   # m∘g^-1
                    k = self.key(ginv2)
                    if k not in dist:
                        dist[k] = d
                        nxt.append((g2, ginv2))
            frontier = nxt
        return dist

    def solve_with_table(self, g: Perm, dist: dict[Perm, int]) -> list[int]:
        ginv = inverse(g)
        d = dist[self.key(ginv)]
        out = []
        while d > 0:
            for i in range(len(self.moves)):
                g2, ginv2 = self.step(g, ginv, i)
                if dist.get(self.key(ginv2)) == d - 1:
                    out.append(i)
                    g, ginv, d = g2, ginv2, d - 1
                    break
            else:
                raise RuntimeError("table inconsistent")
        return out
