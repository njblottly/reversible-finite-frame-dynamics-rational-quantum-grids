(* ::Package:: *)

(* Exact companion checks for Mathematica / Wolfram Language.
   Open a notebook in this folder and evaluate:
   Get["verify_gate_composition.wl"]
   Or: wolframscript -file verify_gate_composition.wl
   It needs no external packages. *)

ClearAll[assert, conditionals, admissibleQ, probabilityGrid, permuteValues,
  treePermutationQ, monomialMatrix, canonical, compose, nearestInteger,
  primitivePhases];

assert[test_, label_] := If[!TrueQ[test], Print["FAILED: ", label]; Abort[]];

conditionals[p_] := Module[{a = p[[1]] + p[[2]]},
  Join[{a}, If[a == 0, {}, {p[[1]]/a}],
    If[a == 1, {}, {p[[3]]/(1 - a)}]]];

admissibleQ[p_, L_, occupation_: False] :=
  Length[p] == 4 && Min[p] >= 0 && Total[p] == 1 &&
  AllTrue[conditionals[p], IntegerQ[L #] &] &&
  (!occupation || AllTrue[p, IntegerQ[L #] &]);

probabilityGrid[L_] := DeleteDuplicates[Flatten[
  Table[With[{a = m/L, b = n/L, c = r/L},
    {a b, a (1 - b), (1 - a) c, (1 - a) (1 - c)}],
    {m, 0, L}, {n, 0, L}, {r, 0, L}], 2]];

(* pi[[j]] is the destination of input basis state j (one-based). *)
permuteValues[p_, pi_] := p[[Ordering[pi]]];
treePermutationQ[pi_] :=
  Sort[Sort /@ Partition[pi, 2]] == {{1, 2}, {3, 4}};

allPermutations = Permutations[Range[4]];
expectedPermutations = Select[allPermutations, treePermutationQ];
assert[Length[expectedPermutations] == 8, "tree group has order eight"];

gridResults = Table[Module[{tg, kg, st, sk},
  tg = probabilityGrid[L];
  kg = Select[tg, admissibleQ[#, L, True] &];
  st = Select[allPermutations, Function[pi,
    AllTrue[tg, admissibleQ[permuteValues[#, pi], L] &]]];
  sk = Select[allPermutations, Function[pi,
    AllTrue[kg, admissibleQ[permuteValues[#, pi], L, True] &]]];
  assert[Sort[st] == Sort[expectedPermutations], "T permutation classification"];
  assert[Sort[sk] == Sort[expectedPermutations], "K permutation classification"];
  assert[Length[tg] == (L + 1) (L^2 + 1), "T probability count"];
  {L, Length[tg], Length[kg], Length[st], Length[sk]}],
  {L, {4, 8, 16}}];
Print["Rows: L, T probabilities, K probabilities, T permutations, K permutations"];
Print[Grid[Prepend[gridResults, {"L", "T", "K", "Perm T", "Perm K"}], Frame -> All]];
assert[gridResults[[All, 3]] == {27, 77, 233}, "K probability counts"];

witness = {1/4, 1/4, 1/2, 0};
swapped = permuteValues[witness, {1, 3, 2, 4}];
assert[conditionals[swapped] == {3/4, 1/3, 1}, "SWAP witness"];
assert[admissibleQ[witness, 8, True] && !admissibleQ[swapped, 8], "SWAP membership"];
Print["SWAP output conditionals: ", conditionals[swapped]];

(* Independent exact matrix check of the normal-form product at L=8. *)
canonical[pi_, d_, L_] := {pi, Mod[d - First[d], L]};
compose[g_, h_, L_] := Module[{pi = g[[1]], d = g[[2]],
    sigma = h[[1]], e = h[[2]]},
  canonical[pi[[sigma]], d + e[[Ordering[pi]]], L]];
monomialMatrix[g_, L_] := Normal[SparseArray[
  Table[{g[[1, j]], j} -> Exp[2 Pi I g[[2, g[[1, j]]]]/L], {j, 4}], {4, 4}]];

g = canonical[{3, 4, 1, 2}, {0, 1, 3, 6}, 8];
h = canonical[{1, 2, 4, 3}, {0, 2, 5, 7}, 8];
productMatrix = FullSimplify[monomialMatrix[g, 8] . monomialMatrix[h, 8]];
directMatrix = monomialMatrix[compose[g, h, 8], 8];
relativeMatrix = FullSimplify[productMatrix . ConjugateTranspose[directMatrix]];
assert[FullSimplify[relativeMatrix - relativeMatrix[[1, 1]] IdentityMatrix[4]] ==
  ConstantArray[0, {4, 4}], "composition equals sequential matrices up to global phase"];

(* Exact circuit, including its inadmissible-intermediate counterexample. *)
hadamard = {{1, 1}, {1, -1}}/Sqrt[2];
ha = KroneckerProduct[hadamard, IdentityMatrix[2]];
cnot = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 0}};
cz = DiagonalMatrix[{1, 1, 1, -1}];
states = FoldList[FullSimplify[#2 . #1] &, {1, 0, 0, 0}, {ha, cnot, cz, cnot, ha}];
probabilities = FullSimplify[Abs[#]^2] & /@ states;
assert[AllTrue[probabilities, admissibleQ[#, 8, True] &], "circuit admissible at every step"];
assert[Last[states] == {0, 0, 1, 0}, "circuit ends in |10>"];
assert[FullSimplify[2 Abs[states[[3, 1]] states[[3, 4]] - states[[3, 2]] states[[3, 3]]]] == 1,
  "intermediate Bell concurrence"];
badInput = {Sqrt[3]/2, 0, 1/2, 0};
badProbability = FullSimplify[Abs[(ha . badInput)[[1]]]^2];
assert[FullSimplify[badProbability - (2 + Sqrt[3])/4] == 0, "Hadamard irrational output"];
Print["Exact circuit probabilities: ", probabilities];

nearestInteger[x_] := Sign[x] Floor[Abs[x] + 1/2];
primitivePhases[a_, b_, c_, L_] := With[
  {qa = nearestInteger[a], qb = nearestInteger[b], qc = nearestInteger[c]},
  Mod[{0, -qb, -qa, -(qa + qb + qc)}, L]];
spectral = Mod[-(nearestInteger /@ {0, 2/5, 2/5, 4/5}), 8];
local = primitivePhases[2/5, 2/5, 0, 8];
assert[spectral == {0, 0, 0, 7} && local == {0, 0, 0, 0}, "spurious interaction witness"];
spuriousState = {1, 1, 1, Exp[-2 Pi I/8]}/2;
spuriousConcurrence = FullSimplify[2 Abs[
  spuriousState[[1]] spuriousState[[4]] - spuriousState[[2]] spuriousState[[3]]]];
assert[FullSimplify[spuriousConcurrence - Sin[Pi/8]] == 0, "spurious concurrence"];
Print["Whole-spectrum phases: ", spectral, "; local-control phases: ", local];
Print["All Mathematica checks passed."];

