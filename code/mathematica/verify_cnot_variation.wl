(* Standalone exact CNOT-matching counterexample verifier.
Load the saved file using Get for source hash and local JSON/text records.
With[{f=SystemDialogInput["FileOpen"]},If[StringQ[f],Get[f]]]
Exact integer, rational and algebraic arithmetic. Mathematica run pending.
The 184320 objects are retained labels, not distinct physical matrices.
*)
BeginPackage["RaQMCNOTVariationVerification`"];
Begin["`Private`"];
Clear[summary];
currentCase="initialisation";
scriptVersion = "0.1";
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
  "cnot_variation_mathematica_" <> DateString[{"Year", "Month", "Day", "T", "Hour", "Minute", "Second"},
    TimeZone -> 0] <> "_" <> StringTake[CreateUUID[], 8] <> "Z"}];
directoryCreated = Quiet[Check[
  CreateDirectory[recordDirectory, CreateIntermediateDirectories -> True], $Failed]];
recordReady = StringQ[directoryCreated] && DirectoryQ[recordDirectory];
checkRecords = {}; textLog = {};
sourceHash = If[loadedFromFile,
  Quiet[Check[FileHash[scriptPath, "SHA256", "HexString"], Null]], Null];
metadata = <|"started_utc" -> utcString[], "script" -> "verify_cnot_variation.wl",
  "script_version" -> scriptVersion, "loaded_source_path" -> scriptPath,
  "execution_mode" -> If[loadedFromFile, "loaded file", "interactive evaluation"],
  "script_sha256" -> sourceHash,
  "source_hash_status" -> If[StringQ[sourceHash], "captured from loaded file",
    "unavailable; no executed source file was identified or hashing failed"],
  "record_directory_source" -> directorySource,
  "record_directory" -> recordDirectory,
  "wolfram_version" -> $Version, "system_id" -> $SystemID,
  "operating_system" -> $OperatingSystem, "system_word_length" -> $SystemWordLength,
  "arithmetic" -> "exact integer, rational and symbolic",
  "scope" -> "Explicit phase-type CNOT rule and probability history; no exhaustive matching search or asymptotic onset claim."|>;
log[x_] := (AppendTo[textLog, ToString[x, InputForm]]; Print[x]);
assert[test_, label_] := (AppendTo[checkRecords, <|"label" -> label, "passed" -> TrueQ[test]|>];
  If[!TrueQ[test], log["FAILED: " <> label]; log["Case: " <> ToString[currentCase,InputForm]]; Throw["failed", "verification"]]);

scanHorizon=64; (* Seven suffices. No full catalogue enumeration is needed. *)
eq[x_,y_] := TrueQ[RootReduce[x-y] == 0] /; !ListQ[x];
eq[x_List,y_List] := Dimensions[x]===Dimensions[y] && And@@Flatten[MapThread[eq,{x,y},ArrayDepth[x]]];
cs[x_] := ToString[x,InputForm];
status=CheckAbort[Catch[Block[{$IterationLimit=Infinity},
 currentCase="gate and error certificates";
 assert[IntegerQ[scanHorizon]&&scanHorizon>=7,"scan horizon at least seven"];
 eye=IdentityMatrix[2];id4=IdentityMatrix[4];xx={{0,1},{1,0}};yy={{0,-I},{I,0}};zz=DiagonalMatrix[{1,-1}];
 hh={{1,1},{1,-1}}/Sqrt[2];tt=DiagonalMatrix[{1,(1+I)/Sqrt[2]}];ss=RootReduce[tt.tt];rr=RootReduce[zz.hh];
 ha=KroneckerProduct[hh,eye];ta=KroneckerProduct[tt,eye];sa=KroneckerProduct[ss,eye];rb=KroneckerProduct[eye,rr];
 cc={{1,0,0,0},{0,1,0,0},{0,0,0,1},{0,0,1,0}};
 paulis=Flatten[Table[KroneckerProduct[p,q],{p,{eye,xx,yy,zz}},{q,{eye,xx,yy,zz}}],1];
 names=Flatten[Table[p<>q,{p,{"I","X","Y","Z"}},{q,{"I","X","Y","Z"}}]];
 signedPauli[mat_] := Module[{coef,pos,j},coef=RootReduce[Table[Tr[p.mat]/4,{p,paulis}]];
 pos=Flatten[Position[coef,Except[0],{1},Heads->False]];
 assert[Length[pos]==1,"single nonzero Pauli coefficient"];j=First[pos];
 assert[MemberQ[{-1,1},coef[[j]]],"unit signed Pauli coefficient"];
 If[coef[[j]]==-1,"-",""]<>names[[j]]];
 multiplier[a_,b_] := If[a==1&&b==0,cc.rb,cc];
 assert[eq[ConjugateTranspose[rr].rr,eye]&&eq[rr,(eye+I yy)/Sqrt[2]],"R=ZH=(I+iY)/sqrt(2) is unitary"];
 assert[eq[rr.xx.ConjugateTranspose[rr],zz]&&eq[rr.zz.ConjugateTranspose[rr],-xx],"quarter-turn Clifford action"];
 assert[eq[ConjugateTranspose[eye-rr].(eye-rr),(2-Sqrt[2])eye]&&eq[Tr[rr],Sqrt[2]],"exact quarter-turn error certificate"];
 Do[ee=KroneckerProduct[If[a==0,eye,tt],If[b==0,eye,tt]];dd=multiplier[a,b];
 assert[eq[ConjugateTranspose[dd].dd,id4],"bijective unitary left multiplier"];
 Do[signedPauli[dd.p.ConjugateTranspose[dd]],{p,paulis}];
 vv=RootReduce[ee.dd.ConjugateTranspose[ee]];rel=RootReduce[ConjugateTranspose[cc].vv];
 If[a==1&&b==0,assert[eq[vv,cc.rb]&&eq[rel,rb],"modified CNOT type"],
 If[b==0,assert[eq[rel,id4],"exact calibrated CNOT type"],
 assert[eq[rel,DiagonalMatrix[{1,1,(1+I)/Sqrt[2],(1-I)/Sqrt[2]}]],"unchanged odd-target CNOT spectrum"]]],
 {a,0,1},{b,0,1}];
 assert[eq[multiplier[0,0],cc]&&eq[cc.cc,id4],"CNOT calibration I <-> C"];
 assert[eq[cc.cc.ha,ha],"CNOT calibration H_A <-> C H_A"];
 assert[eq[ha.ha,id4],"Hadamard calibration"];
 hrel=RootReduce[ConjugateTranspose[hh].tt.hh.ConjugateTranspose[tt]];
 assert[eq[ConjugateTranspose[eye-hrel].(eye-hrel),(1-1/Sqrt[2])eye],"unchanged Hadamard error bound"];
 currentCase="exact probability histories";
 expectedStabilisers={{"ZI","IZ"},{"XX","-IX"},{"-YI","-IX"},{"YX","-ZZ"},{"-ZX","YY"},{"YY","XZ"},{"-IY","-XY"},{"-ZY","-IY"}};
 expectedIdeal={0,1/2,1/2,1/4,1/4,3/8,5/8,1/16};
 gg=id4;idealMatrix=id4;physical=id4;phaseA=0;actual={};ideal={};joint={};stabiliserRows={};first=Null;
 nominal=cc.ha.ta;
 Do[frame=RootReduce[If[phaseA==0,id4,ta].gg];
 assert[eq[frame,physical],"frame carries equal direct physical composition"];
 q=RootReduce[Conjugate[frame[[{3,4},1]]].frame[[{3,4},1]]];
 p=RootReduce[Conjugate[idealMatrix[[{3,4},1]]].idealMatrix[[{3,4},1]]];
 assert[(IntegerQ[q]||Head[q]===Rational)&&(IntegerQ[p]||Head[p]===Rational),"rational probabilities"];
 AppendTo[actual,q];AppendTo[ideal,p];AppendTo[joint,RootReduce[Conjugate[frame[[All,1]]] frame[[All,1]]]];
 If[first===Null&&Abs[q-p]>1/3,first=n];
 If[n<=7,row={signedPauli[gg.paulis[[13]].ConjugateTranspose[gg]],signedPauli[gg.paulis[[4]].ConjugateTranspose[gg]]};
 assert[row===expectedStabilisers[[n+1]],"hand-derived stabilisers at depth "<>ToString[n]];
 assert[p==expectedIdeal[[n+1]],"ideal probability at depth "<>ToString[n]];AppendTo[stabiliserRows,row]];
 oldA=phaseA;phaseA=1-phaseA;carry=If[oldA==0,id4,sa];
 gg=RootReduce[multiplier[phaseA,0].ha.carry.gg];ee=If[phaseA==0,id4,ta];
 vh=RootReduce[ee.ha.ConjugateTranspose[ee]];vc=RootReduce[ee.multiplier[phaseA,0].ConjugateTranspose[ee]];
 physical=RootReduce[vc.vh.ta.physical];idealMatrix=RootReduce[nominal.idealMatrix],{n,0,scanHorizon}];
 assert[Take[actual,8]=={0,1/2,1/2,1/2,1/2,1/2,1/2,0},"counterexample history through seven"];
 assert[Max[Abs[Take[actual-ideal,8]]]==1/4,"whole-history maximum discrepancy is one quarter"];
 assert[Abs[actual[[8]]-ideal[[8]]]==1/16,"seven-block discrepancy is one sixteenth"];
 assert[joint[[8]]=={1/2,1/2,0,0},"joint probabilities at seven"];
 assert[first===Null||first>7,"no strict one-third crossing through seven"];
 summary=<|"counterexample_verified"->True,"modified_phase_type"->{1,0},
 "CNOT_Clifford_multiplier_at_modified_type"->"C (I tensor ZH)","other_phase_types"->"C",
 "H_error_bound"->"sqrt(1-1/sqrt(2))","C_error_bound"->"sqrt(2-sqrt(2))",
 "retained_labels"->184320,"additional_labels_relative_to_previous_stage"->0,
 "calibration_preserved"->True,"inverse_C_rule"->"left multiplication by D_ab^dagger; phase types unchanged",
 "actual_A1_n_0_to_7"->(cs /@ Take[actual,8]),"ideal_A1_n_0_to_7"->(cs /@ Take[ideal,8]),
 "maximum_discrepancy_through_seven"->"1/4","discrepancy_at_seven"->"1/16","joint_at_seven"->(cs /@ joint[[8]]),
 "signed_inner_stabilisers_n_0_to_7"->stabiliserRows,"scan_horizon"->scanHorizon,
 "first_gap_above_one_third_within_scan"->first,"actual_A1_full_scan"->(cs /@ actual),"ideal_A1_full_scan"->(cs /@ ideal),
 "all_fixed_phase_priors"->"analytic consequence of inherited right-phase equivariance and compensated preparation",
 "certification_scope"->"explicit full-table algebraic rule; no exhaustive optimisation or asymptotic onset claim"|>;
 "passed"],"verification"],"aborted"];
If[!StringQ[status],status="failed"];If[!AssociationQ[summary],summary=<||>];
record=Join[metadata,<|"finished_utc"->utcString[],"status"->status,"total_assertions"->Length[checkRecords],
 "passed_assertions"->Count[Lookup[checkRecords,"passed"],True],"results"->summary,"checks"->checkRecords|>];
If[recordReady,jp=FileNameJoin[{recordDirectory,"cnot_variation_run_record.json"}];tp=FileNameJoin[{recordDirectory,"cnot_variation_assertion_log.txt"}];
 jr=Quiet[Check[Export[jp,record,"RawJSON"],$Failed]];
 tr=Quiet[Check[Export[tp,StringRiffle[(If[TrueQ[#["passed"]],"PASS: ","FAIL: "]<>#["label"]& /@ checkRecords),"\n"],"Text"],$Failed]];
 exportStatus=If[StringQ[jr]&&StringQ[tr],"saved","failed"],exportStatus="failed"];
Print["Mathematical verification: ",status];Print["Record export: ",exportStatus];Print["Local record: ",recordDirectory];
If[!loadedFromFile,Print["Source hash unavailable. Load the saved .wl file using Get to capture it."]];
End[];EndPackage[];
