(* ::Package:: *)

(* Nominal-gate verification, logging revision 0.2 (8 July 2026).
   Recommended: save this file, then evaluate the following in a notebook
   and choose this .wl file in the file picker:
   With[{f = SystemDialogInput["FileOpen"]}, If[StringQ[f], Get[f]]]
   Results and environment metadata are exported automatically to a new
   reproducibility/mathematica_... directory beside this file.
   No external packages; all definitions use a dedicated context. *)

Begin["RaQMNominalVerification`"];
Clear[allocationRows, localRows, coverageRows, residuals];
scriptVersion = "0.2";
inputFile = $InputFileName;
loadedFromFile = StringQ[inputFile] && StringLength[inputFile] > 0 &&
  FileExistsQ[inputFile];
scriptPath = If[loadedFromFile, ExpandFileName[inputFile], Null];
If[loadedFromFile,
  baseDirectory = DirectoryName[scriptPath];
  directorySource = "loaded script directory",
  notebookDirectory = If[$FrontEnd === Null, $Failed,
    Quiet[Check[NotebookDirectory[], $Failed]]];
  If[StringQ[notebookDirectory] && DirectoryQ[notebookDirectory],
    baseDirectory = notebookDirectory; directorySource = "saved notebook directory",
    baseDirectory = Directory[]; directorySource = "kernel working directory"]
];
utcString[] := DateString[{"Year", "-", "Month", "-", "Day", "T",
  "Hour", ":", "Minute", ":", "Second", "Z"}, TimeZone -> 0];
recordDirectory = FileNameJoin[{baseDirectory, "reproducibility",
  "mathematica_" <> DateString[{"Year", "Month", "Day", "T", "Hour", "Minute", "Second"},
    TimeZone -> 0] <> "_" <> StringTake[CreateUUID[], 8] <> "Z"}];
directoryCreated = Quiet[Check[
  CreateDirectory[recordDirectory, CreateIntermediateDirectories -> True], $Failed]];
recordReady = StringQ[directoryCreated] && DirectoryQ[recordDirectory];
checkRecords = {}; textLog = {};
sourceHash = If[loadedFromFile,
  Quiet[Check[FileHash[scriptPath, "SHA256", "HexString"], Null]], Null];
metadata = <|"started_utc" -> utcString[], "script" -> "verify_nominal_gates.wl",
  "script_version" -> scriptVersion, "loaded_source_path" -> scriptPath,
  "execution_mode" -> If[loadedFromFile, "loaded file", "interactive evaluation"],
  "script_sha256" -> sourceHash,
  "source_hash_status" -> If[StringQ[sourceHash], "captured from loaded file",
    "unavailable; no executed source file was identified or hashing failed"],
  "record_directory_source" -> directorySource,
  "record_directory" -> recordDirectory,
  "wolfram_version" -> $Version, "system_id" -> $SystemID,
  "operating_system" -> $OperatingSystem, "system_word_length" -> $SystemWordLength,
  "numeric_working_precision" -> 50, "numeric_tolerance" -> "1e-25",
  "scope" -> "Exact arithmetic checks and separately identified numerical matrix checks; not interval certification."|>;
log[x_] := (AppendTo[textLog, ToString[x, InputForm]]; Print[x]);
assert[test_, label_] := (AppendTo[checkRecords, <|"label" -> label, "passed" -> TrueQ[test]|>];
  If[!TrueQ[test], log["FAILED: " <> label]; Throw["failed", "verification"]]);

conditionals[p_] := Module[{a = p[[1]] + p[[2]]},
  Join[{a}, If[a == 0, {}, {p[[1]]/a}], If[a == 1, {}, {p[[3]]/(1 - a)}]]];
admissibleQ[p_, L_, occupation_: False] := Total[p] == 1 && Min[p] >= 0 &&
  AllTrue[conditionals[p], IntegerQ[L #] &] &&
  (!occupation || AllTrue[p, IntegerQ[L #] &]);
probabilityGrid[A_, B_] := DeleteDuplicates[Flatten[
  Table[With[{a = m/A, b = n/B, c = r/B},
    {a b, a (1 - b), (1 - a) c, (1 - a) (1 - c)}],
    {m, 0, A}, {n, 0, B}, {r, 0, B}], 2]];
aligned[phi_, chi_] := Module[{z = Conjugate[phi] . chi},
  If[Abs[z] == 0, chi, chi Conjugate[z]/Abs[z]]];
chordal[phi_, chi_] := Norm[phi - aligned[phi, chi]];
correction[phi_, target_] := Module[{chi = aligned[phi, target], c, s, e},
  If[Norm[phi - chi] < 10^-35, Return[IdentityMatrix[Length[phi]]]];
  c = Min[1, Max[0, Re[Conjugate[phi] . chi]]];
  e = chi - c phi; s = Norm[e]; e = e/s;
  IdentityMatrix[Length[phi]] +
    (c - 1) (Outer[Times, phi, Conjugate[phi]] + Outer[Times, e, Conjugate[e]]) +
    s (Outer[Times, e, Conjugate[phi]] - Outer[Times, phi, Conjugate[e]])];
gridQuantise[phi_, L_, kind_] := Module[{A, B, m, n, r, a, b, c, p, ref, phases, probs},
  If[kind == "T", A = L; B = L,
    A = 2^Floor[Log[2, L]/2]; B = L/A];
  p = Abs[phi]^2; p = p/Total[p]; a = p[[1]] + p[[2]];
  b = If[a == 0, 0, p[[1]]/a];
  c = If[p[[3]] + p[[4]] == 0, 0, p[[3]]/(p[[3]] + p[[4]])];
  m = Clip[Round[A a], {0, A}]; n = Clip[Round[B b], {0, B}]; r = Clip[Round[B c], {0, B}];
  {a, b, c} = {m/A, n/B, r/B};
  probs = {a b, a (1 - b), (1 - a) c, (1 - a) (1 - c)};
  ref = First[Select[phi, Abs[#] > 0 &]];
  phases = If[Abs[#] == 0, 0, Mod[Round[L Arg[#/ref]/(2 Pi)], L]] & /@ phi;
  {N[Sqrt[probs] Exp[2 Pi I phases/L], 50], probs, {A, B}}];

status = CheckAbort[Catch[
  allocationRows = Table[Module[{tg, kg, union},
    tg = probabilityGrid[L, L]; kg = Select[tg, admissibleQ[#, L, True] &];
    union = DeleteDuplicates[Join @@ Table[probabilityGrid[2^r, L/2^r], {r, 0, Log[2, L]}]];
    assert[Sort[union] == Sort[kg], "exact allocation union L=" <> ToString[L]];
    {L, Length[kg], Length[union]}], {L, {4, 8, 16}}];
  log["Allocation rows: L, K count, union count"]; log[allocationRows];

  localRows = Table[Module[{b = (L/2 + 1)/L, allowed},
    allowed = Select[Range[0, L]/L, Function[a,
      admissibleQ[{a b, a (1 - b), (1 - a) b, (1 - a) (1 - b)}, L, True]]];
    assert[allowed == {0, 1}, "exact local obstruction L=" <> ToString[L]];
    {L, b, allowed}], {L, {4, 8, 16, 64, 256}}];
  log["Local obstruction rows: L, fixed B probability, allowed A probabilities"]; log[localRows];
  h = {{1, 1}, {1, -1}}/Sqrt[2]; z = DiagonalMatrix[{1, -1}];
  assert[FullSimplify[(ConjugateTranspose[z - h] . (z - h)) -
    (2 - Sqrt[2]) IdentityMatrix[2]] == ConstantArray[0, {2, 2}], "exact sharp local norm gap"];

  SeedRandom[20260907]; tolerance = 10^-25; residuals = {};
  coverageRows = Flatten[Table[Module[{targets, largest = 0, chi, probs, allocation, bound, rot, residual},
    targets = Join[IdentityMatrix[4], Table[
      Normalize[N[RandomInteger[{-20, 20}, 4] + I RandomInteger[{-20, 20}, 4], 50]], {24}]];
    Do[
      {chi, probs, allocation} = gridQuantise[N[phi, 50], L, kind];
      assert[admissibleQ[probs, L, kind == "K"], "exact output probability certificate"];
      bound = N[Sqrt[Total[1/allocation]] + 2 Sin[Pi/(2 L)], 50];
      assert[chordal[phi, chi] <= bound + tolerance, "numerical covering bound"];
      rot = correction[N[phi, 50], chi];
      residual = Norm[ConjugateTranspose[rot] . rot - IdentityMatrix[4], "Frobenius"];
      AppendTo[residuals, residual];
      assert[residual < tolerance, "numerical correction unitarity"];
      assert[Norm[rot . phi - aligned[phi, chi]] < tolerance, "numerical correction action"];
      largest = Max[largest, chordal[phi, chi]], {phi, targets}];
    {L, kind, N[largest, 12], N[bound, 12]}], {L, {4, 8, 16, 64, 256}}, {kind, {"T", "K"}}], 1];
  log["Coverage rows: L, grid, largest observed distance, proved bound"]; log[coverageRows];

  (* Exact local recovery on an entangled, fixed-marginal example. *)
  f = h . DiagonalMatrix[{Sqrt[3]/2, 1/2}];
  w = DiagonalMatrix[{1, I}]; cmat = w . f;
  recovered = FullSimplify[cmat . Inverse[f]];
  assert[recovered == w, "exact fixed-marginal local recovery"];
  assert[FullSimplify[ConjugateTranspose[cmat] . cmat - ConjugateTranspose[f] . f] ==
    ConstantArray[0, {2, 2}], "exact fixed B marginal"];
  log["All nominal-gate Mathematica checks passed."];
  "passed", "verification"], "aborted"];

AssociateTo[metadata, {"finished_utc" -> utcString[], "status" -> status,
  "checks" -> checkRecords}];
If[ValueQ[allocationRows], AssociateTo[metadata, "allocation_rows" -> allocationRows]];
If[ValueQ[coverageRows], AssociateTo[metadata, "coverage_rows" -> coverageRows]];
If[ValueQ[residuals] && Length[residuals] > 0,
  AssociateTo[metadata, "maximum_unitarity_residual" -> N[Max[residuals], 16]]];
jsonPath = FileNameJoin[{recordDirectory, "mathematica_run_record.json"}];
logPath = FileNameJoin[{recordDirectory, "mathematica_assertion_log.txt"}];
jsonExport = $Failed; logExport = $Failed;
If[recordReady,
  jsonExport = Quiet[Check[Export[jsonPath, metadata, "RawJSON"], $Failed]];
  logExport = Quiet[Check[Export[logPath, StringRiffle[textLog, "\n"], "Text"], $Failed]]
];
exportSucceeded = StringQ[jsonExport] && StringQ[logExport] &&
  FileExistsQ[jsonPath] && FileExistsQ[logPath];
Print["Mathematical verification: ", status];
Print["Record export: ", If[exportSucceeded, "saved", "FAILED"]];
If[exportSucceeded,
  Print["Local record: ", recordDirectory],
  Print["The complete local record was not saved. Intended directory: ", recordDirectory]
];
If[!StringQ[sourceHash],
  Print["Source hash unavailable. To capture it, load the saved .wl file using Get."]
];
End[];
