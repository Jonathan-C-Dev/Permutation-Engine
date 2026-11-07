import math
import random

import pytest

from cosets import ColoringSpace
from permgroup import StabChain, compose, from_cycles, identity, inverse, perm_order
from puzzles import Cube
from thistlethwaite import G_GENS, Thistlethwaite

RNG = random.Random(0)


def test_symmetric_and_alternating_orders():
    for n in range(3, 9):
        s = StabChain(n, [from_cycles(n, [range(n)]), from_cycles(n, [(0, 1)])])
        assert s.order() == math.factorial(n)
    n = 7  # A_7 from 3-cycles
    a = StabChain(n, [from_cycles(n, [(0, 1, i)]) for i in range(2, n)])
    assert a.order() == math.factorial(n) // 2


def test_mathieu_m12():
    g1 = from_cycles(12, [(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)])
    g2 = from_cycles(12, [(2, 6, 10, 7), (3, 9, 4, 5)])
    g3 = from_cycles(12, [(0, 11), (1, 10), (2, 5), (3, 7), (4, 8), (6, 9)])
    assert StabChain(12, [g1, g2, g3]).order() == 95040


def test_membership_and_random_elements():
    c = Cube(3)
    ch = StabChain(c.N, c.gens("URFDLB"))
    assert ch.order() == 43252003274489856000
    for _ in range(20):
        assert ch.contains(ch.random_element(RNG))
    # a single swapped pair of stickers is not a legal cube state
    bad = list(c.id); bad[0], bad[1] = bad[1], bad[0]
    assert not ch.contains(tuple(bad))


def test_moves_have_right_order_and_group_laws():
    c = Cube(3)
    for f in "URFDLB":
        assert perm_order(c.moves[f]) == 4
        assert compose(c.moves[f], c.moves[f + "'"]) == c.id
    # commutator of R and U has order 6 on a 3x3
    comm = c.apply(c.id, "R U R' U'")
    assert perm_order(comm) == 6


def test_2x2_order():
    c = Cube(2, "URF")
    assert StabChain(c.N, c.gens(["U", "R", "F"])).order() == 3674160


def test_canon_left_is_coset_invariant():
    c = Cube(3)
    H = StabChain(c.N, c.gens(G_GENS[3]))
    G = StabChain(c.N, c.gens(G_GENS[2]))
    for _ in range(10):
        x = G.random_element(RNG)
        h = H.random_element(RNG)
        assert H.canon_left(compose(x, h)) == H.canon_left(x)       # same coset, same key
    x, y = G.random_element(RNG), G.random_element(RNG)
    same = H.contains(compose(inverse(x), y))
    assert (H.canon_left(x) == H.canon_left(y)) == same


def test_lagrange_indices_and_coloring_counts():
    T = Thistlethwaite(cache_dir=None)
    info = T.verify(deep=False, log=lambda *_: None)
    assert info["indices"] == [2048, 1082565, 29400, 663552]


def test_full_solver_roundtrip():
    T = Thistlethwaite(cache_dir=".cache")
    ch = StabChain(T.cube.N, T.cube.gens(G_GENS[0]))
    for _ in range(3):
        s = ch.random_element(RNG)
        sol = sum(T.solve(s), [])
        assert T.cube.apply(s, sol) == T.cube.id
    s, _ = T.cube.scramble(20, RNG)
    assert T.cube.apply(s, sum(T.solve(s), [])) == T.cube.id


def test_2x2_optimal_solver():
    c = Cube(2, "URF")
    names = list(c.moves)
    sp = ColoringSpace(c.label_identity(), [c.moves[n] for n in names])
    ch = StabChain(c.N, c.gens(["U", "R", "F"]))
    for _ in range(3):
        s = ch.random_element(RNG)
        seq = [names[i] for i in sp.solve(s)]
        assert c.apply(s, seq) == c.id and len(seq) <= 14
