# Twisted primitive group association schemes

Computational files accompanying the paper
[**Twisted primitive group association schemes**](https://arxiv.org/abs/2608.16278)
by Akihiro Higashitani, Masanari Kamiya, and Hirotake Kurihara.

The paper studies whether the intersection numbers of a primitive group
association scheme determine it up to combinatorial isomorphism. Its main
theoretical construction gives, for a broad infinite family of groups
$\operatorname{PSL}(2,q)$, Cayley association schemes that are algebraically
isomorphic to the corresponding group association schemes but are not
combinatorially isomorphic to them. In particular, these group association
schemes are non-separable. The construction applies when $q$ is an odd prime
power with $q=11$ or $q\geq 17$, or when $q=2^f$ with $f\geq 3$.

This repository contains the GAP computations used for the finite cases in the
paper:

- explicit twists of the group association schemes of $\mathfrak A_6$ and
  $\mathfrak A_8$, together with witnesses for non-isomorphism and
  non-Schurity;
- rigidity checks for $\mathfrak A_4,\mathfrak A_5,\mathfrak A_7$, and
  $\mathfrak A_9$ under the paper's two-class twisting criterion; and
- an exhaustive search over all groups of order at most $200$.

The $\operatorname{PSL}(2,q)$ construction is proved theoretically in the
paper; this repository does not contain a separate script for that construction.

## Repository contents

| Script | Purpose | Direct package dependency | Output |
| --- | --- | --- | --- |
| [`A4.g`](A4.g) | Rigidity of the split 3-cycle classes of $\mathfrak A_4$ | GRAPE | [`out/result_A4.txt`](out/result_A4.txt) |
| [`A5.g`](A5.g) | Rigidity of the split 5-cycle classes of $\mathfrak A_5$ | GRAPE | [`out/result_A5.txt`](out/result_A5.txt) |
| [`A6.g`](A6.g) | Construction and verification of the $\mathfrak A_6$ twist, including the product identities, local graphs, triangle counts, and exhaustive check of the ten $3+3$ cut partitions | GRAPE | [`out/result_A6.txt`](out/result_A6.txt), [`out/A6_twisted_Schur_partition.txt`](out/A6_twisted_Schur_partition.txt) |
| [`checktwistA6.g`](checktwistA6.g) | Independent full relation-matrix check of algebraic isomorphism, combinatorial non-isomorphism, Schurity, and automorphism groups for $\mathfrak A_6$ | AssociationSchemes | [`out/result_checktwistA6.txt`](out/result_checktwistA6.txt) |
| [`A7.g`](A7.g) | Rigidity of the split 7-cycle classes of $\mathfrak A_7$, using a nullity computation over $\mathrm{GF}(1009)$ | GRAPE | [`out/result_A7.txt`](out/result_A7.txt) |
| [`A8.g`](A8.g) | Construction and verification of the $\mathfrak A_8$ twist without forming a full $20160\times20160$ relation matrix | GRAPE | [`out/result_A8.txt`](out/result_A8.txt), [`out/A8_twisted_Schur_partition.txt`](out/A8_twisted_Schur_partition.txt) |
| [`checktwistA8.g`](checktwistA8.g) | Optional dense relation-matrix check of algebraic isomorphism and Schurity for $\mathfrak A_8$ | AssociationSchemes | `out/result_checktwistA8.txt` (generated, not included) |
| [`A9.g`](A9.g) | Character-table, class-algebra, and centralizer-orbit quotient checks for the split 9-cycle and split $(5,3)$-classes of $\mathfrak A_9$ | cvec | [`out/result_A9.txt`](out/result_A9.txt) |
| [`small_groups_twist.g`](small_groups_twist.g) | Exhaustive search for two-class twists among all groups of order at most $200$ | SmallGrp, AssociationSchemes | [`out/result_small_groups_twist.txt`](out/result_small_groups_twist.txt), [`out/result_small_groups_twist.csv`](out/result_small_groups_twist.csv) |

The checked-in files under [`out/`](out/) are the reference outputs used in the
paper. The two `*_twisted_Schur_partition.txt` files are GAP-readable
definitions of the constructed partitions.

## Requirements

The computations in the paper were run with
[GAP 4.15.1](https://www.gap-system.org/). The scripts do not pin package
versions. The paper cites AssociationSchemes 3.1.0, SmallGrp 1.5.4, and cvec
2.8.4. GRAPE 4.9.3 was available in the verification environment used for this
README.

Install the packages needed by the scripts you intend to run:

- [GRAPE](https://gap-packages.github.io/grape/)
- [AssociationSchemes](https://www.jesselansdown.com/AssociationSchemes/)
- [SmallGrp](https://gap-packages.github.io/smallgrp/)
- [cvec](https://gap-packages.github.io/cvec/)

You can check that GAP can load all four packages with:

```bash
gap -q -c 'for p in ["grape", "AssociationSchemes", "SmallGrp", "cvec"] do Print(p, ": ", LoadPackage(p) <> fail, "\n"); od;'
```

The scripts call `mkdir -p` to create `out/`, so a POSIX-like environment is
assumed. On Windows, use WSL or create `out/` before running the scripts.

## Running the computations

Run all commands from the repository root because output paths are relative to
the current working directory.

### Quick start: the $\mathfrak A_6$ twist

```bash
git clone https://github.com/Kurihara190/Twisted-primitive-group-association-schemes.git
cd Twisted-primitive-group-association-schemes

gap -q -c 'Read("A6.g");'
gap -q -c 'Read("checktwistA6.g");'
```

The first command constructs the twisted partition and verifies the identities
and local-graph calculations used in the proof. The second command reads that
partition and independently confirms, using AssociationSchemes, that the two
schemes are algebraically isomorphic but not combinatorially isomorphic, and
that the twisted scheme is non-Schurian.

### Alternating-group checks

```bash
gap -q -c 'Read("A4.g");'
gap -q -c 'Read("A5.g");'
gap -q -c 'Read("A6.g");'
gap -q -c 'Read("checktwistA6.g");'
gap -q -c 'Read("A7.g");'
gap -q -c 'Read("A8.g");'
gap -q -c 'Read("A9.g");'
```

Run `A6.g` before `checktwistA6.g`. Similarly, `A8.g` must be run before the
optional `checktwistA8.g`, since each check script reads the partition generated
by the corresponding construction script.

### Small-group search

```bash
gap -q -c 'Read("small_groups_twist.g");'
```

This searches all $6{,}065$ isomorphism classes of groups of order at most
$200$. The CSV contains one row for each conjugacy-class pair satisfying
conditions (C1)--(C4), together with the linear-search statistics, numbers of
replacement partitions, and combinatorial-isomorphism flags. The text file
contains the final summary.

Running the scripts regenerates files under `out/`. Use `git diff -- out/` to
compare a new run with the checked-in reference results. GAP may wrap long
printed lists differently depending on the terminal width; these formatting
differences do not change the underlying GAP objects.

## Checked-in results at a glance

- For $\mathfrak A_6$ and $\mathfrak A_8$, the constructed Schur partitions
  are algebraically isomorphic, but not combinatorially isomorphic, to the
  respective conjugacy-class partitions. The resulting twisted association
  schemes are non-Schurian.
- For $\mathfrak A_4,\mathfrak A_5,\mathfrak A_7$, and $\mathfrak A_9$, the
  relevant one-dimensional eigenspace checks establish rigidity under the
  two-class twisting criterion used in the paper.
- The search through order $200$ finds $1{,}672$ groups with at least one
  class pair satisfying (C1)--(C4), $556$ non-original partitions satisfying
  (D1)--(D4), and $256$ partitions whose Cayley schemes are not
  combinatorially isomorphic to the original group association scheme. These
  new twists occur for $48$ group isomorphism types, with
  `SmallGroup(108,15)` as the smallest example.

## Resource notes

`A8.g` is the intended verification for the $\mathfrak A_8$ results: it works
directly with permutations and avoids constructing the full relation matrix.
In contrast, `checktwistA8.g` explicitly builds a $20160\times20160$ dense
matrix ($406{,}425{,}600$ entries), so it can require very large amounts of
memory and time. Its combinatorial-isomorphism call is intentionally disabled
in the current script, and no output from this optional check is committed.

The $\mathfrak A_9$ quotient calculations and the exhaustive small-group
search can also take substantially longer than the smaller alternating-group
checks.

## Citation

If you use this code, please cite the accompanying paper:

```bibtex
@misc{higashitani2026twisted,
  title         = {Twisted primitive group association schemes},
  author        = {Higashitani, Akihiro and Kamiya, Masanari and Kurihara, Hirotake},
  year          = {2026},
  eprint        = {2608.16278},
  archivePrefix = {arXiv},
  primaryClass  = {math.CO},
  url           = {https://arxiv.org/abs/2608.16278}
}
```
