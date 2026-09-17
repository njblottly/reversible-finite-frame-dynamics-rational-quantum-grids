(* ::Package:: *)

(* Exact CNOT and local-frame benchmark v0.1. Mathematica 13.3.
   Standalone: enumerates 11520 Clifford frames using exact signed Pauli actions.
   With[{f=SystemDialogInput["FileOpen"]},If[StringQ[f],Get[f]]]
   Saves JSON and a text log under reproducibility/ beside the saved script.
   All verification arithmetic is exact. *)
BeginPackage["RaQMEntanglingControlVerification`"];
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
  "entangling_controls_mathematica_" <> DateString[{"Year", "Month", "Day", "T", "Hour", "Minute", "Second"},
    TimeZone -> 0] <> "_" <> StringTake[CreateUUID[], 8] <> "Z"}];
directoryCreated = Quiet[Check[
  CreateDirectory[recordDirectory, CreateIntermediateDirectories -> True], $Failed]];
recordReady = StringQ[directoryCreated] && DirectoryQ[recordDirectory];
checkRecords = {}; textLog = {};
sourceHash = If[loadedFromFile,
  Quiet[Check[FileHash[scriptPath, "SHA256", "HexString"], Null]], Null];
metadata = <|"started_utc" -> utcString[], "script" -> "verify_entangling_controls.wl",
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
  "scope" -> "Exact Clifford benchmark only; the refining H/T/CNOT construction is an analytic theorem, not a computed table."|>;
log[x_] := (AppendTo[textLog, ToString[x, InputForm]]; Print[x]);
assert[test_, label_] := (AppendTo[checkRecords, <|"label" -> label, "passed" -> TrueQ[test]|>];
  If[!TrueQ[test], log["FAILED: " <> label]; log["Case: " <> ToString[currentCase,InputForm]]; Throw["failed", "verification"]]);


(* Hermitian Pauli order II,IX,IY,IZ,XI,...,ZZ. Signed indices are one-based. *)
id=Range[16];
compose[p_,q_] := Sign[q] p[[Abs[q]]];
inverse[p_] := Module[{v=ConstantArray[0,16]},v[[Abs[p]]]=Sign[p] Range[16];v];
key[p_] := FromDigits[p+16,33];
act[p_,r_] := Module[{v=ConstantArray[0,16]},v[[Abs[p]]]=Sign[p] r;v];
pow[p_,n_Integer] := Nest[compose[p,#]&,id,n];
local[p_,side_] := Flatten[Table[If[side==1,
 Sign[p[[a]]] (4(Abs[p[[a]]]-1)+b),
 Sign[p[[b]]] (4(a-1)+Abs[p[[b]]])],{a,4},{b,4}]];
(* Fixed buffer prevents quadratic list growth during exhaustive enumeration. *)
enumerate[generators_] := Module[{out=ConstantArray[0,{12000,16}],seen=<||>,n=1,head=1,x,y,ky},
 out[[1]]=id;AssociateTo[seen,key[id]->1];
 While[head<=n,x=out[[head]];
  Do[y=compose[g,x];ky=key[y];If[!KeyExistsQ[seen,ky],n++;
    If[n>12000,Throw["group exceeded safety buffer","verification"]];
    out[[n]]=y;AssociateTo[seen,ky->n]],{g,generators}];head++];
 {Take[out,n],seen}];
state[rules_] := Normal[SparseArray[Join[{1->1},rules],16]];
mean[rs_] := Total[rs]/Length[rs];
status=CheckAbort[Catch[
 currentCase="exact matrix generator validation";
 eye=IdentityMatrix[2];sx={{0,1},{1,0}};sy={{0,-I},{I,0}};sz=DiagonalMatrix[{1,-1}];
 h={{1,1},{1,-1}}/Sqrt[2];s=DiagonalMatrix[{1,I}];
 cm={{1,0,0,0},{0,1,0,0},{0,0,0,1},{0,0,1,0}};
 paulis=Flatten[Table[KroneckerProduct[a,b],{a,{eye,sx,sy,sz}},{b,{eye,sx,sy,sz}}],1];
 names={"HA","HB","SA","SB","C"};
 mats={KroneckerProduct[h,eye],KroneckerProduct[eye,h],KroneckerProduct[s,eye],KroneckerProduct[eye,s],cm};
 actions=Table[Table[
    coeffs=FullSimplify[Table[Tr[paulis[[j]] . mats[[k]] . paulis[[i]] . ConjugateTranspose[mats[[k]]]]/4,{j,16}]];
    positions=Flatten[Position[coeffs,Except[0],{1},Heads->False]];
    assert[Length[positions]==1 && MemberQ[{-1,1},coeffs[[First[positions]]]],"matrix conjugation is a signed Pauli image"];
    First[positions] coeffs[[First[positions]]],{i,16}],{k,5}];
 gens=AssociationThread[names,actions];ha=gens["HA"];hb=gens["HB"];sa=gens["SA"];sb=gens["SB"];c=gens["C"];
 assert[ha==local[{1,4,-3,2},1] && hb==local[{1,4,-3,2},2],"independent local H action"];
 assert[sa==local[{1,3,-2,4},1] && sb==local[{1,3,-2,4},2],"independent local S action"];
 assert[compose[c,c]==id && c[[5]]==6 && c[[4]]==16,"CNOT involution and XI/IZ images"];
 currentCase="complete group and local cosets";log["Enumerating exact Clifford frames..."];
 {group,index}=enumerate[actions];{locals,localIndex}=enumerate[Take[actions,4]];
 assert[Length[group]==11520 && Length[locals]==576,"complete and local group orders"];
 reps={};coordinates=<||>;unique=True;
 Do[If[!KeyExistsQ[coordinates,key[g]],AppendTo[reps,g];j=Length[reps];
   Do[frame=compose[locals[[k]],g];ky=key[frame];
    If[KeyExistsQ[coordinates,ky],unique=False];AssociateTo[coordinates,ky->{k,j}],{k,Length[locals]}]],{g,group}];
 assert[unique && Length[reps]==20 && Length[coordinates]==11520,"unique 20-coset factorisation"];
 assert[First[reps]==id,"identity preparation chart"];
 tables=<||>;
 Do[p=gens[name];table=Lookup[index,key[compose[p,#]]]& /@ group;AssociateTo[tables,name->table];
  assert[Sort[table]==Range[11520],name<>" complete transition permutation"];
  assert[And@@Table[compose[inverse[p],group[[table[[i]]]]]==group[[i]],{i,11520}],name<>" exact inverse on all frames"];
  If[name!="C",
   assert[And@@Table[Last[Lookup[coordinates,key[group[[table[[i]]]]]]]==Last[Lookup[coordinates,key[group[[i]]]]],{i,11520}],name<>" retains nonlocal chart"];
   assert[And@@Table[compose[group[[table[[i]]]],inverse[group[[i]]]]==p,{i,11520}],name<>" strictly local physical action"]],{name,names}];
 ct=tables["C"];assert[ct[[ct]]==Range[11520],"CNOT echo restores all frames"];
 assert[And@@(compose[ha,compose[hb,#]]==compose[hb,compose[ha,#]]& /@ group),"opposite Hadamards commute on all frames"];
 assert[!KeyExistsQ[localIndex,key[c]] && !KeyExistsQ[localIndex,key[compose[c,compose[ha,c]]]],"CNOT and conjugated local H are nonlocal"];
 currentCase="compensated preparation and entangling word channels";
 zero=state[{4->1,13->1,16->1}];phi=state[{6->1,11->-1,16->1}];mixed=state[{}];
 plusplus=state[{2->1,5->1,6->1}];plusi0=state[{4->1,9->1,12->1}];
 assert[act[compose[c,ha],zero]==phi,"H_A then CNOT creates Bell state"];
 assert[phi[[{2,3,4,5,9,13}]]==ConstantArray[0,6],"Bell marginals exactly maximally mixed"];
 Do[inp=state[{13->(-1)^a,4->(-1)^b,16->(-1)^(a+b)}];bb=BitXor[a,b];
  out=state[{13->(-1)^a,4->(-1)^bb,16->(-1)^(a+bb)}];assert[act[c,inp]==out,"CNOT computational truth table"],{a,0,1},{b,0,1}];
 pa={id,local[{1,2,-3,-4},1],local[{1,-2,3,-4},1],local[{1,-2,-3,4},1]};
 bell=act[#,phi]& /@ pa;rotated=act[sa,#]& /@ bell;
 assert[mean[bell]==mixed && mean[rotated]==mixed,"equal-density Bell preparation recipes"];
 words={{"HA","C"},{"HA","C","SB","HB","C","SA"},{"HB","SA","C","HA","SB","C","HB"}};
 phases=Flatten[Table[compose[pow[sa,k],pow[sb,l]],{k,0,3},{l,0,3}],1];rows={};
 Do[u=Fold[compose[gens[#2],#1]&,id,word];
  assert[And@@(compose[inverse[u],compose[u,#]]==#& /@ group),"entangling inverse word restores every frame"];
  averages={};
  Do[frame=compose[u,f0];w=compose[frame,inverse[f0]];
   assert[w==u,"compensated common unitary independent of phase"];
   Do[seed=act[inverse[f0],r];direct=act[frame,seed];
    assert[direct==act[w,r],"direct compensation agrees with common channel"];
    assert[act[inverse[frame],direct]==seed,"seed recoverable after entangling word"],
    {r,Join[{zero,phi,plusplus,plusi0},bell,rotated]}];
   Do[AppendTo[averages,mean[act[frame,act[inverse[f0],#]]& /@ recipe]],{recipe,{bell,rotated}}],{f0,phases}];
  assert[And@@(#==mixed& /@ averages),"equal-density recipes agree in every branch"];
  rr=act[u,zero];probs=Flatten[Table[(1+(-1)^a rr[[13]]+(-1)^b rr[[4]]+(-1)^(a+b) rr[[16]])/4,{a,0,1},{b,0,1}]];
  assert[Total[probs]==1 && Min[probs]>=0,"normalised exact computational readout"];
  AppendTo[rows,<|"chronological_word"->word,"probabilities"->(ToString[#,InputForm]& /@ probs)|>],{word,words}];
 summary=<|"group_size"->11520,"local_subgroup_size"->576,"left_cosets"->20,"D"->2,"L"->4,
  "single_qubit_frames"->24,"initial_phase_pairs"->16,"complete_generator_transition_checks"->57600,
  "probability_rows"->rows,"scope"->"Exact Clifford benchmark only; no non-Clifford CNOT table or numerical quotient radius certified."|>;
 log[summary];log["All exact entangling-controller benchmark checks passed."];"passed","verification"],"aborted"];
AssociateTo[metadata,{"finished_utc"->utcString[],"status"->status,"checks"->checkRecords,
 "passed_assertions"->Count[Lookup[checkRecords,"passed"],True],"total_assertions"->Length[checkRecords],"last_case"->currentCase}];
If[ValueQ[summary],AssociateTo[metadata,"summary"->summary]];
jsonPath=FileNameJoin[{recordDirectory,"entangling_controls_run_record.json"}];logPath=FileNameJoin[{recordDirectory,"entangling_controls_assertion_log.txt"}];
jsonExport=$Failed;logExport=$Failed;
If[recordReady,jsonExport=Quiet[Check[Export[jsonPath,metadata,"RawJSON"],$Failed]];
 logExport=Quiet[Check[Export[logPath,StringRiffle[textLog,"\n"]<>"\n"<>
 StringRiffle[Map[If[TrueQ[#["passed"]],"PASS: ","FAIL: "]<>#["label"]&,checkRecords],"\n"],"Text"],$Failed]]];
exportSucceeded=StringQ[jsonExport] && StringQ[logExport] && FileExistsQ[jsonPath] && FileExistsQ[logPath];
Print["Mathematical verification: ",status];Print["Record export: ",If[exportSucceeded,"saved","FAILED"]];Print["Local record: ",recordDirectory];
If[!StringQ[sourceHash],Print["Source hash unavailable. Load the saved .wl file using Get to capture it."]];
End[];
EndPackage[];
