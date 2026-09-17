(* ::Package:: *)

(* RaQM context preparation and probabilities, v0.1. Mathematica 13.3.
   Save this file and load it using:
   With[{f=SystemDialogInput["FileOpen"]},If[StringQ[f],Get[f]]]
   Fresh JSON and text records are saved in reproducibility/ beside this file. *)
BeginPackage["RaQMPreparationVerification`"];
Begin["`Private`"];
Clear[summary,fibre,output];
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
  "preparation_mathematica_" <> DateString[{"Year", "Month", "Day", "T", "Hour", "Minute", "Second"},
    TimeZone -> 0] <> "_" <> StringTake[CreateUUID[], 8] <> "Z"}];
directoryCreated = Quiet[Check[
  CreateDirectory[recordDirectory, CreateIntermediateDirectories -> True], $Failed]];
recordReady = StringQ[directoryCreated] && DirectoryQ[recordDirectory];
checkRecords = {}; textLog = {};
sourceHash = If[loadedFromFile,
  Quiet[Check[FileHash[scriptPath, "SHA256", "HexString"], Null]], Null];
metadata = <|"started_utc" -> utcString[], "script" -> "verify_preparation_law.wl",
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
  "scope" -> "Exact finite checks of specified context preparation laws and terminal Born probabilities; not a unique physical prediction of RaQM."|>;
log[x_] := (AppendTo[textLog, ToString[x, InputForm]]; Print[x]);
assert[test_, label_] := (AppendTo[checkRecords, <|"label" -> label, "passed" -> TrueQ[test]|>];
  If[!TrueQ[test], log["FAILED: " <> label]; log["Case: " <> ToString[currentCase,InputForm]]; Throw["failed", "verification"]]);

tilde[a_] := reflection . Transpose[a] . reflection;
swap[a_] := Which[a==ident,phase,a==phase,ident,True,a];
obs[x_] := x[[1]] . tilde[x[[2]]];
action[x_,aa_,bb_] := {If[aa==1,swap[x[[1]]],x[[1]]],If[bb==1,swap[x[[2]]],x[[2]]]};
fibre[r_] := fibre[r]=Select[labels,obs[#]==r&];
signature[x_] := Table[obs[action[x,ab[[1]],ab[[2]]]],{ab,{{0,0},{1,0},{0,1},{1,1}}}];
overlap[r_,s_] := (1+Tr[Transpose[r] . s])/4;
safeLabels[r_] := Select[fibre[r],!MemberQ[{ident,phase},#[[1]]] && !MemberQ[{ident,phase},#[[2]]]&];
output[r_,aa_,bb_,safe_:False] := Mean[obs[action[#,aa,bb]]& /@ If[TrueQ[safe],safeLabels[r],fibre[r]]];
status=CheckAbort[Catch[
 currentCase="exact conditional preparation and quotient measure";
 ident=IdentityMatrix[3];reflection=DiagonalMatrix[{1,-1,1}];
 phase={{0,-1,0},{1,0,0},{0,0,1}};
 group=Sort[Select[Flatten[Table[DiagonalMatrix[sg] . ident[[pm]],
   {pm,Permutations[Range[3]]},{sg,Tuples[{-1,1},3]}],1],Det[#]==1&]];
 assert[Length[DeleteDuplicates[group]]==24,"24 distinct proper rotations"];
 labels=Tuples[group,2];signatures=signature /@ labels;
 classes=GatherBy[Range[Length[labels]],signatures[[#]]&];
 assert[Length[classes]==116,"116 retained context classes"];
 classSize[key_] := Count[signatures,key];
 Do[assert[signature[action[action[x,1,0],1,0]]==signature[x],"exact context echo"];
   Do[assert[classSize[signature[x]]==classSize[signature[action[x,ab[[1]],ab[[2]]]]],
     "quotient counting measure is invariant"],{ab,{{1,0},{0,1},{1,1}}}],{x,labels}];
 Do[assert[Length[fibre[r]]==24,"24 frame representatives per ray"];
   assert[Mean[obs /@ fibre[r]]==r,"specified pure ray is prepared"];
   one=output[r,1,0];
   expected=(22 r+phase . r+Transpose[phase] . r)/24;
   assert[one==expected,"one-sided law is a random-unitary channel"];
   assert[overlap[r,one]==23/24,"one-step survival is 23/24"];
   assert[Length[safeLabels[r]]>=20,"stationary-support preparation exists"];
   localClasses=Tally[signature /@ fibre[r]];
   Do[assert[output[r,ab[[1]],ab[[2]],True]==r,"stationary-support identity law"];
     pos=First[FirstPosition[{{0,0},{1,0},{0,1},{1,1}},ab]];
     weighted=Total[(#[[2]]/24 #[[1,pos]])& /@ localClasses];
     assert[weighted==output[r,ab[[1]],ab[[2]]],"induced quotient weights preserve statistics"],
     {ab,{{0,0},{1,0},{0,1},{1,1}}}],{r,group}];
 phiClasses=Tally[signature /@ fibre[ident]];
 assert[Sort[Last /@ phiClasses]=={1,1,1,21},"Phi context multiplicities"];
 wrong=Mean[#[[1,2]]& /@ phiClasses];
 assert[overlap[ident,wrong]==3/4,"uniform quotient weights change prediction to 3/4"];
 reset=Mean[output[obs[action[#,1,0]],1,0]& /@ fibre[ident]];
 assert[overlap[ident,reset]==265/288,"context resampling loses the exact echo"];
 currentCase="equal-density preparations with different two-sided outputs";
 pauliRotations={ident,DiagonalMatrix[{1,-1,-1}],DiagonalMatrix[{-1,1,-1}],DiagonalMatrix[{-1,-1,1}]};
 ensembles={pauliRotations,phase . #& /@ pauliRotations};
 Do[assert[Mean[en]==ConstantArray[0,{3,3}],"both ensembles prepare I/4"];
   assert[Mean[output[#,1,0]& /@ en]==ConstantArray[0,{3,3}],"one-sided mixture consistency"],{en,ensembles}];
 t0=Mean[output[#,1,1]& /@ ensembles[[1]]];
 t1=Mean[output[#,1,1]& /@ ensembles[[2]]];
 expected={{0,1/24,0},{-1/24,0,0},{0,0,0}};
 assert[t0==expected && t1==-expected,"two-sided exact correlation obstruction"];
 assert[{(1-t0[[1,2]])/2,(1-t1[[1,2]])/2}=={23/48,25/48},"binary outcome probabilities"];
 sx={{0,1},{1,0}};sy={{0,-I},{I,0}};sz={{1,0},{0,-1}};paulis={sx,sy,sz};
 rho[t_] := (IdentityMatrix[4]+Sum[t[[j,k]] KroneckerProduct[paulis[[j]],Transpose[paulis[[k]]]],{j,3},{k,3}])/4;
 effect=(IdentityMatrix[4]+KroneckerProduct[sx,sy])/2;
 assert[{Tr[effect . rho[t0]],Tr[effect . rho[t1]]}=={23/48,25/48},"independent exact Pauli-matrix probability check"];
 assert[Sort[Eigenvalues[rho[t0]]]=={11/48,1/4,1/4,13/48} && Sort[Eigenvalues[rho[t1]]]=={11/48,1/4,1/4,13/48},
   "both output density matrices are positive and normalised"];
 currentCase="finite fair-bit sampling";
 Do[n=2^bits;{q,rem}=QuotientRemainder[n,m];counts=Join[ConstantArray[q+1,rem],ConstantArray[q,m-rem]];
   tv=Total[Abs[counts/n-ConstantArray[1/m,m]]]/2;
   assert[Total[counts]==n && tv==rem (m-rem)/(m n),"balanced finite-bit variation formula"];
   assert[tv<=m/(2 n),"finite-bit error bound"];
   assert[(tv==0)==(Mod[n,m]==0),"exact uniformity requires divisibility"],{m,1,32},{bits,0,8}];
 summary=<|"label_count"->576,"retained_context_classes"->116,"phi_context_multiplicities"->{1,1,1,21},
   "one_step_survival"->"23/24","inverse_echo"->"1","two_steps_with_context_reset"->"265/288",
   "uniform_quotient_one_step"->"3/4","equal_density_mixture_probabilities"->{"23/48","25/48"},
   "mixture_probability_difference"->"1/24","scope"->"Conditional predictions, not a unique RaQM law"|>;
 log[summary];log["All preparation-law checks passed."];"passed","verification"],"aborted"];
AssociateTo[metadata,{"finished_utc"->utcString[],"status"->status,"checks"->checkRecords,
 "passed_assertions"->Count[Lookup[checkRecords,"passed"],True],"total_assertions"->Length[checkRecords],"last_case"->currentCase}];
If[ValueQ[summary],AssociateTo[metadata,"summary"->summary]];
jsonPath=FileNameJoin[{recordDirectory,"preparation_law_run_record.json"}];
logPath=FileNameJoin[{recordDirectory,"preparation_law_assertion_log.txt"}];
jsonExport=$Failed;logExport=$Failed;
If[recordReady,
 jsonExport=Quiet[Check[Export[jsonPath,metadata,"RawJSON"],$Failed]];
 logExport=Quiet[Check[Export[logPath,StringRiffle[textLog,"\n"]<>"\n"<>
 StringRiffle[Map[If[TrueQ[#["passed"]],"PASS: ","FAIL: "]<>#["label"]&,checkRecords],"\n"],"Text"],$Failed]]];
exportSucceeded=StringQ[jsonExport] && StringQ[logExport] && FileExistsQ[jsonPath] && FileExistsQ[logPath];
Print["Mathematical verification: ",status];Print["Record export: ",If[exportSucceeded,"saved","FAILED"]];
Print["Local record: ",recordDirectory];
If[!StringQ[sourceHash],Print["Source hash unavailable. Load the saved .wl file using Get to capture it."]];
End[];
EndPackage[];
