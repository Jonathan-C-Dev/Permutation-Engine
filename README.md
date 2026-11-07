# groupsolver: permutation-group engine + coset-space cube solver

Pure Python, no dependencies (pytest for tests). Everything is derived from
permutation groups: puzzle moves are generated from sticker geometry, group
orders come from Schreier-Sims, and the solver searches **coset spaces** whose
sizes are checked against Lagrange's theorem.

## Layout
| file | what it does |
|---|---|
| `permgroup.py` | permutations, cycle decomposition, **Schreier-Sims** (BSGS), order, membership, uniform random elements, canonical coset representatives |
| `puzzles.py` | 2x2x2 / 3x3x3 move permutations generated from 3D rotations of stickers; sticker colorings that define subgroups |
| `cosets.py` | bidirectional BFS on coset graphs; `ColoringSpace` (fast) and `CanonCosetSpace` (any subgroup) |
| `thistlethwaite.py` | 4-phase solver along G0 > G1 > G2 > G3 > 1, plus `verify()` |
| `main.py` | CLI: `verify`, `solve3`, `solve2`, `group` |
| `test_all.py` | known group orders (S_n, A_7, Mathieu M12, Rubik), coset-key invariance, round-trip solves |

## Run
    python main.py verify --deep     # |G_i|, Lagrange indices, BFS coset counts (~1 min)
    python main.py solve3 --n 20     # random-state 3x3 solves
    python main.py solve2 --n 5      # optimal 2x2 solves
    python main.py group "(0 1 2 3 4)" "(0 1)"
    pytest -q

## The math
* **Lagrange / orbit-stabilizer**: |G| = prod of orbit sizes down a stabilizer chain, so Schreier-Sims
  returns exact orders (Rubik: 43,252,003,274,489,856,000 in ~0.25 s).
* **Cosets**: to move a state g into H, only the right coset Hg matters because Hg·w = H <=> g·w in H.
  Right multiplication by a move acts on right cosets, so each phase is a BFS on a graph of size [G_i : G_{i+1}].
* **Coloring cosets**: if H = {states that look solved under a sticker coloring}, then Hg <-> label∘g and moves
  just permute the coloring. Used for edge orientation, corner twist + slice edges, and (with the identity coloring) the last phase.
* **General cosets**: for subgroups that are not coloring stabilizers (phase 3), a coset is named by a canonical element
  of g^-1·H computed from H's stabilizer chain (min image at each level, multiply by that transversal element).
* **Self-checking**: `verify()` recomputes every |G_i| and confirms the BFS coset counts equal |G_i|/|G_{i+1}|
  (2,048 / 1,082,565 / 29,400), so a wrong coordinate definition fails loudly.

## Measured (CPython 3.12, one core)
* Schreier-Sims on the full cube group: 0.25 s (base length 18)
* Phase-3 table (29,400 cosets): built once (~10 s), cached to `.cache/`
* Random-state 3x3 solves: average 31.6 moves over 20 states, typically 10-60 ms per solve after the table is cached
* 2x2x2 (3,674,160 states): optimal solutions, 4-20 ms
