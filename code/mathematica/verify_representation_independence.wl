(* ::Package:: *)

(* RaQM representation independence v0.1. Mathematica 13.3 compatible.
   Select this SAVED file from Mathematica with:
   With[{f=SystemDialogInput["FileOpen"]},If[StringQ[f],Get[f]]]
   Produces JSON and a text assertion log under reproducibility/ beside it.
   Exact arithmetic only. *)
BeginPackage["RaQMRepresentationVerification`"];
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
  "representation_mathematica_" <> DateString[{"Year", "Month", "Day", "T", "Hour", "Minute", "Second"},
    TimeZone -> 0] <> "_" <> StringTake[CreateUUID[], 8] <> "Z"}];
directoryCreated = Quiet[Check[
  CreateDirectory[recordDirectory, CreateIntermediateDirectories -> True], $Failed]];
recordReady = StringQ[directoryCreated] && DirectoryQ[recordDirectory];
checkRecords = {}; textLog = {};
sourceHash = If[loadedFromFile,
  Quiet[Check[FileHash[scriptPath, "SHA256", "HexString"], Null]], Null];
metadata = <|"started_utc" -> utcString[], "script" -> "verify_representation_independence.wl",
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
  "scope" -> "Exact checks on the 24-frame Clifford Bell sector and rational moment identities; not exhaustive verification at all resolutions."|>;
log[x_] := (AppendTo[textLog, ToString[x, InputForm]]; Print[x]);
assert[test_, label_] := (AppendTo[checkRecords, <|"label" -> label, "passed" -> TrueQ[test]|>];
  If[!TrueQ[test], log["FAILED: " <> label]; log["Case: " <> ToString[currentCase,InputForm]]; Throw["failed", "verification"]]);

idx[m_] := First[FirstPosition[group,m]];
tilde[m_] := reflection . Transpose[m] . reflection;
colour[values_] := Module[{blocks=GatherBy[Range[Length[values]],values[[#]]&],out},
 out=ConstantArray[0,Length[values]];
 Do[out[[blocks[[k]]]]=ConstantArray[k,Length[blocks[[k]]]],{k,Length[blocks]}];out];
labelIndex[a_,b_] := 24(a-1)+b;
status=CheckAbort[Catch[
 currentCase="exact Clifford group and Bell labels";
 ident=IdentityMatrix[3];reflection=DiagonalMatrix[{1,-1,1}];
 phase={{0,-1,0},{1,0,0},{0,0,1}};
 group=Sort[Select[Flatten[Table[DiagonalMatrix[sg] . ident[[pm]],
   {pm,Permutations[Range[3]]},{sg,Tuples[{-1,1},3]}],1],Det[#]==1&]];
 assert[Length[DeleteDuplicates[group]]==24,"24 proper signed permutation rotations"];
 Do[assert[MemberQ[group,tilde[a]] && tilde[tilde[a]]==a,"transpose involution"];
   Do[assert[MemberQ[group,a . b],"exact group closure"];
     assert[tilde[a . b]==tilde[b] . tilde[a],"transpose reverses products"],{b,group}],{a,group}];
 labels=Tuples[Range[24],2];
 obs=Table[idx[group[[pair[[1]]]] . tilde[group[[pair[[2]]]]]],{pair,labels}];
 assert[Length[DeleteDuplicates[obs]]==24 && Union[Last /@ Tally[obs]]=={24},
   "24 initial rays with 24 labels each"];
 ii=idx[ident];si=idx[phase];
 swap=Table[Which[k==ii,si,k==si,ii,True,k],{k,24}];
 fa=Table[labelIndex[swap[[pair[[1]]]],pair[[2]]],{pair,labels}];
 fb=Table[labelIndex[pair[[1]],swap[[pair[[2]]]]],{pair,labels}];
 Do[assert[fa[[fa[[k]]]]==k && fb[[fb[[k]]]]==k,"local involutions"];
   assert[fa[[fb[[k]]]]==fb[[fa[[k]]]],"local commutation"],{k,576}];
 x=labelIndex[si,ii];y=labelIndex[ii,si];
 assert[obs[[x]]==obs[[y]] && obs[[fa[[x]]]]!=obs[[fa[[y]]]],"Bell aliases separate"];
 assert[Transpose[group[[obs[[fa[[x]]]]]]] . group[[obs[[fa[[y]]]]]]==phase . phase,
   "witness relative output rotation is phase squared"];
 Do[f=Table[idx[v . a],{a,group}];
   expected=Table[idx[v . group[[o]]],{o,obs}];
   actual=Table[idx[group[[f[[pair[[1]]]]]] . tilde[group[[pair[[2]]]]]],{pair,labels}];
   assert[actual==expected,"common left multiplication descends on all Bell labels"],{v,group}];
 currentCase="exact minimum future-observation quotient";
 current=colour[obs];blockCounts={Length[Union[current]]};stable=False;
 While[!stable,
   next=colour[Table[{current[[k]],current[[fa[[k]]]],current[[fb[[k]]]]},{k,576}]];
   stable=(next===current);
   If[!stable,current=next;AppendTo[blockCounts,Length[Union[current]]]]];
 oracle=colour[Table[obs[[{k,fa[[k]],fb[[k]],fa[[fb[[k]]]]}]],{k,576}]];
 assert[current==oracle,"refinement equals exhaustive four-word observation oracle"];
 quotientSize=Length[Union[current]];
 Do[qmap=ConstantArray[0,quotientSize];
   Do[If[qmap[[current[[k]]]]==0,qmap[[current[[k]]]]=current[[transition[[k]]]]];
     assert[qmap[[current[[k]]]]==current[[transition[[k]]]],"quotient well defined"],{k,576}];
   assert[Sort[qmap]==Range[quotientSize],"quotient is a permutation"],{transition,{fa,fb}}];
 fibreCounts=Table[Length[Union[Pick[current,obs,r]]],{r,24}];
 assert[blockCounts=={24,116} && Sort[fibreCounts]==Join[{3,4,4},ConstantArray[5,21]],
   "independent expected quotient and context counts"];
 currentCase="exact moment identities";
 Do[l=2^power;Do[d=2^dp;
   z2=2+l Sum[(2 n/d-1)^2,{n,1,d-1}];
   x2=l/2 Sum[1-(2 n/d-1)^2,{n,1,d-1}];
   assert[z2==2+l(d-1)(d-2)/(3 d),"axial moment"];
   assert[x2==l (d^2-1)/(3 d),"transverse moment"];
   assert[z2-x2==2-l+l/d,"moment difference"];
   assert[(z2==x2)==(d==2 && l==4),"unique sampled isotropic exception"],
   {dp,1,power}],{power,2,8}];
 currentCase="symbolic two-qubit witness and catalogue identification";
 bell={1,0,0,1}/Sqrt[2];sGate=DiagonalMatrix[{1,I}];
 inA=KroneckerProduct[sGate,IdentityMatrix[2]] . bell;
 inB=KroneckerProduct[IdentityMatrix[2],sGate] . bell;
 outA=bell;outB=KroneckerProduct[sGate,sGate] . bell;
 assert[inA==inB && Conjugate[outA] . outB==0,"exact equal input and orthogonal outputs at L=4"];
 remote[psi_] := Module[{m=Partition[psi,2]},Transpose[ConjugateTranspose[m] . m]];
 assert[remote[outA]==IdentityMatrix[2]/2 && remote[outB]==IdentityMatrix[2]/2,
   "both witness remote marginals are maximally mixed"];
 param[t_,b_,g_] := {{Sqrt[t],Sqrt[1-t]Exp[I b]},
   {Sqrt[1-t]Exp[I g],-Sqrt[t]Exp[I(b+g)]}};
 cat=FullSimplify[Join[Table[param[0,0,Pi k/2],{k,0,3}],
   Table[param[1,0,Pi k/2],{k,0,3}],
   Flatten[Table[param[1/2,Pi b/2,Pi g/2],{b,0,3},{g,0,3}],1]]];
 pauli={{{0,1},{1,0}},{{0,-I},{I,0}},{{1,0},{0,-1}}};
 rotation[u_] := FullSimplify[Table[Tr[pauli[[j]] . u . pauli[[k]] . ConjugateTranspose[u]]/2,{j,3},{k,3}]];
 assert[Sort[rotation /@ cat]==group,"C_2,4 is exactly the 24-rotation Clifford catalogue"];
 summary=<|"label_count"->576,"initial_ray_count"->24,"refinement_block_counts"->blockCounts,
   "quotient_class_count"->quotientSize,"context_counts_by_ray"->Sort[fibreCounts],
   "maximum_context_count"->Max[fibreCounts],"conditional_fixed_width_bits"->Ceiling[Log[2,Max[fibreCounts]]]|>;
 log[summary];log["All representation-independence checks passed."];
 "passed","verification"],"aborted"];
AssociateTo[metadata,{"finished_utc"->utcString[],"status"->status,"checks"->checkRecords,
 "passed_assertions"->Count[Lookup[checkRecords,"passed"],True],"total_assertions"->Length[checkRecords],
 "last_case"->currentCase}];
If[ValueQ[summary],AssociateTo[metadata,"summary"->summary]];
jsonPath=FileNameJoin[{recordDirectory,"representation_independence_run_record.json"}];
logPath=FileNameJoin[{recordDirectory,"representation_independence_assertion_log.txt"}];
jsonExport=$Failed;logExport=$Failed;
If[recordReady,
 jsonExport=Quiet[Check[Export[jsonPath,metadata,"RawJSON"],$Failed]];
 logExport=Quiet[Check[Export[logPath,StringRiffle[textLog,"\n"]<>"\n"<>
   StringRiffle[Map[If[TrueQ[#["passed"]],"PASS: ","FAIL: "]<>#["label"]&,checkRecords],"\n"],"Text"],$Failed]]];
exportSucceeded=StringQ[jsonExport] && StringQ[logExport] && FileExistsQ[jsonPath] && FileExistsQ[logPath];
Print["Mathematical verification: ",status];
Print["Record export: ",If[exportSucceeded,"saved","FAILED"]];
Print["Local record: ",recordDirectory];
If[!StringQ[sourceHash],Print["Source hash unavailable. Load the saved .wl file using Get to capture it."]];
End[];
EndPackage[];
