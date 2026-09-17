(* ::Package:: *)

(* Standalone exact probability-history verification.
Load the saved file using Get for source hashing and local JSON/text records.
With[{f=SystemDialogInput["FileOpen"]},If[StringQ[f],Get[f]]]
Certifies 192-block postponement and 193-block obstruction in the full
right-phase-equivariant local C_8,8 matching class.
*)
BeginPackage["RaQMProbabilityHistoryVerification`"];
Begin["`Private`"];
Clear[summary];
currentCase="initialisation";
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
  "probability_histories_mathematica_" <> DateString[{"Year", "Month", "Day", "T", "Hour", "Minute", "Second"},
    TimeZone -> 0] <> "_" <> StringTake[CreateUUID[], 8] <> "Z"}];
directoryCreated = Quiet[Check[
  CreateDirectory[recordDirectory, CreateIntermediateDirectories -> True], $Failed]];
recordReady = StringQ[directoryCreated] && DirectoryQ[recordDirectory];
checkRecords = {}; textLog = {};
sourceHash = If[loadedFromFile,
  Quiet[Check[FileHash[scriptPath, "SHA256", "HexString"], Null]], Null];
metadata = <|"started_utc" -> utcString[], "script" -> "verify_probability_histories.wl",
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
  "scope" -> "Exact finite-history witness and exhaustive obstruction for the phase-equivariant local catalogue; no full two-qubit or asymptotic onset claim."|>;
log[x_] := (AppendTo[textLog, ToString[x, InputForm]]; Print[x]);
assert[test_, label_] := (AppendTo[checkRecords, <|"label" -> label, "passed" -> TrueQ[test]|>];
  If[!TrueQ[test], log["FAILED: " <> label]; log["Case: " <> ToString[currentCase,InputForm]]; Throw["failed", "verification"]]);

witnessH={212, 268, 260, 252, 244, 236, 228, 220, 240, 248, 256, 264, 208, 216, 224, 232, 276, 397, 261, 125, 180, 108, 163, 348, 284, 341, 269, 133, 188, 116, 171, 356, 292, 349, 213, 141, 196, 124, 179, 364, 300, 357, 221, 85, 204, 132, 187, 372, 308, 365, 229, 93, 148, 140, 195, 380, 316, 373, 237, 101, 156, 84, 203, 388, 324, 381, 245, 109, 164, 92, 147, 396, 332, 389, 253, 117, 172, 100, 155, 340, 411, 461, 325, 189, 61, 43, 227, 283, 419, 405, 333, 197, 69, 51, 235, 291, 427, 413, 277, 205, 77, 59, 243, 299, 435, 421, 285, 149, 21, 67, 251, 307, 443, 429, 293, 157, 29, 75, 259, 315, 451, 437, 301, 165, 37, 19, 267, 323, 459, 445, 309, 173, 45, 27, 211, 331, 403, 453, 317, 181, 53, 35, 219, 275, 404, 390, 198, 70, 52, 107, 290, 347, 412, 398, 206, 78, 60, 115, 298, 355, 420, 342, 150, 22, 68, 123, 306, 363, 428, 350, 158, 30, 76, 131, 314, 371, 436, 358, 166, 38, 20, 139, 322, 379, 444, 366, 174, 46, 28, 83, 330, 387, 452, 374, 182, 54, 36, 91, 274, 395, 460, 382, 190, 62, 44, 99, 282, 339, 12, 454, 262, 134, 0, 98, 226, 354, 13, 462, 270, 142, 7, 106, 234, 362, 14, 406, 214, 86, 6, 114, 242, 370, 15, 414, 222, 94, 5, 122, 250, 378, 8, 422, 230, 102, 4, 130, 258, 386, 9, 430, 238, 110, 3, 138, 266, 394, 10, 438, 246, 118, 2, 82, 210, 338, 11, 446, 254, 126, 1, 90, 218, 346, 432, 447, 326, 79, 16, 34, 162, 418, 440, 455, 334, 23, 24, 42, 170, 426, 448, 463, 278, 31, 32, 50, 178, 434, 456, 407, 286, 39, 40, 58, 186, 442, 400, 415, 294, 47, 48, 66, 194, 450, 408, 423, 302, 55, 56, 74, 202, 458, 416, 431, 310, 63, 64, 18, 146, 402, 424, 439, 318, 71, 72, 26, 154, 410, 425, 383, 327, 207, 25, 89, 225, 361, 433, 391, 335, 151, 33, 97, 233, 369, 441, 399, 279, 159, 41, 105, 241, 377, 449, 343, 287, 167, 49, 113, 249, 385, 457, 351, 295, 175, 57, 121, 257, 393, 401, 359, 303, 183, 65, 129, 265, 337, 409, 367, 311, 191, 73, 137, 209, 345, 417, 375, 319, 199, 17, 81, 217, 353, 304, 376, 263, 143, 144, 88, 161, 297, 312, 384, 271, 87, 152, 96, 169, 305, 320, 392, 215, 95, 160, 104, 177, 313, 328, 336, 223, 103, 168, 112, 185, 321, 272, 344, 231, 111, 176, 120, 193, 329, 280, 352, 239, 119, 184, 128, 201, 273, 288, 360, 247, 127, 192, 136, 145, 281, 296, 368, 255, 135, 200, 80, 153, 289}+1;
sgn[x_] := Sign[RootReduce[x]];
phaseIndex[a_Integer] := 8 Quotient[a-1,8]+Mod[a,8]+1;
rightPhase[a_Integer] := Module[{i=a-1,n,b,g},Which[i<8,Mod[i-1,8]+1,
 i<16,9+Mod[i-8+1,8],True,n=Quotient[i-16,64];b=Quotient[Mod[i-16,64],8];g=Mod[i-16,8];17+64 n+8 Mod[b+1,8]+g]];
restrictGraph[graph_,rules_List] := Module[{aa=First /@ rules,bb=Last /@ rules},
 Table[If[MemberQ[aa,i],{i /. rules},Complement[graph[[i]],bb]],{i,Length[graph]}]];
(* Explicit success flag: Return inside Do would exit only the loop. *)
fullMatching[graph_] := Module[{nn=Length[graph],owner,augment,seen,p},
 owner=ConstantArray[0,nn];
 augment[i_] := Module[{j,k=1,found=False},
  While[k<=Length[graph[[i]]] && !found,
   j=graph[[i,k]];
   If[!seen[[j]],seen[[j]]=True;
    If[owner[[j]]==0 || TrueQ[augment[owner[[j]]]],owner[[j]]=i;found=True]];
   k++];found];
 Do[seen=ConstantArray[False,nn];augment[i],{i,nn}];p=ConstantArray[0,nn];
 Do[If[owner[[j]]>0,p[[owner[[j]]]]=j],{j,nn}];p];
verifyHall[graph_,p_] := Module[{nn=Length[graph],owner,left,nb={},queue,i,j},
 owner=ConstantArray[0,nn];Do[If[p[[i]]>0,owner[[p[[i]]]]=i],{i,nn}];
 left=Flatten[Position[p,0]];queue=left;
 While[Length[queue]>0,i=First[queue];queue=Rest[queue];
 Do[If[!MemberQ[nb,j],AppendTo[nb,j];If[owner[[j]]>0 && !MemberQ[left,owner[[j]]],
 AppendTo[left,owner[[j]]];AppendTo[queue,owner[[j]]]]],{j,graph[[i]]}]];
 Length[nb]<Length[left] && Sort[nb]==Union[Flatten[graph[[left]]]]];
cycle[p_,start_] := Module[{out={},x=start},While[!MemberQ[out,x],AppendTo[out,x];x=p[[x]]];
 assert[x==start,"cycle closes at its start"];out];
(* One shared assignment and fixed probability bands over the whole history.
   Each recursive call has its own Catch tag; successful search cannot be
   swallowed by Do, and verification failures retain their separate tag. *)
historySearch[graph_,bands_,start_,anchors_,horizon_] := Module[
 {nodes=0,newAssignments=0,bandReject=0,injectionReject=0,hallReject=0,longest=0,explore,survives},
 explore[path_,assigned_] := Module[{n=Length[path],x=Last[path],choices,y,trial,mat,successTag},Catch[
  nodes++;longest=Max[longest,n-1];If[n>horizon,Throw[True,successTag]];
  choices=If[MemberQ[First /@ assigned,x],{x/.assigned},graph[[x]]];
  Do[If[!bands[[n+1,y]],bandReject++;Continue[]];
   If[MemberQ[First /@ assigned,x],If[TrueQ[explore[Append[path,y],assigned]],Throw[True,successTag]],
    If[MemberQ[Last /@ assigned,y],injectionReject++;Continue[]];
    trial=Append[assigned,x->y];mat=fullMatching[restrictGraph[graph,trial]];
    If[MemberQ[mat,0],assert[verifyHall[restrictGraph[graph,trial],mat],"independent Hall obstruction for search prune"];hallReject++;Continue[]];
    newAssignments++;If[TrueQ[explore[Append[path,y],trial]],Throw[True,successTag]]],{y,choices}];
  False,successTag]];
 survives=If[TrueQ[bands[[1,start]]],explore[{start},anchors],False];
 <|"survives"->survives,"longest"->longest,"counts"->{nodes,newAssignments,bandReject,injectionReject,hallReject}|>];
status=CheckAbort[Catch[Block[{$RecursionLimit=10000,$IterationLimit=Infinity},
 currentCase="small matching and shared-history regression checks";
 assert[fullMatching[{{1,2},{1}}]=={2,1},"augmenting path reassigns an earlier owner"];
 Do[tinyGraph=Table[Select[Range[3],BitAnd[mask,2^(3(i-1)+#-1)]!=0&],{i,3}];
  tinyMatch=fullMatching[tinyGraph];
  brute=AnyTrue[Permutations[Range[3]],Function[perm,And@@Table[MemberQ[tinyGraph[[i]],perm[[i]]],{i,3}]]];
  assert[(Sort[tinyMatch]==Range[3])===brute,"small matching agrees with exhaustive permutations"];
  If[!brute,assert[verifyHall[tinyGraph,tinyMatch],"small failed matching has a Hall witness"]],{mask,0,511}];
 tinyResult=historySearch[{{1,2},{1,2}},{{True,False},{False,True},{True,False}},1,{},2];
 assert[TrueQ[tinyResult["survives"]],"positive shared-history branch propagates success"];
 tinyResult=historySearch[{{1,2},{1,2}},{{True,False},{False,True},{False,True}},1,{},2];
 assert[tinyResult["survives"]===False,"shared bijection rejects a noninjective history"];
 Print["Small matching and history regression checks passed."];
 currentCase="exact catalogue, phase action and quotient graph";
 phases={1,(1+I)/Sqrt[2],I,(-1+I)/Sqrt[2],-1,(-1-I)/Sqrt[2],-I,(1-I)/Sqrt[2]};
 param[n_,b_,g_] := {{Sqrt[n/8],Sqrt[1-n/8] phases[[1+Mod[b,8]]]},
 {Sqrt[1-n/8] phases[[1+Mod[g,8]]],-Sqrt[n/8] phases[[1+Mod[b+g,8]]]}};
 cat=Join[Table[param[0,0,k],{k,0,7}],Table[param[8,0,k],{k,0,7}],
 Flatten[Table[param[n,b,g],{n,1,7},{b,0,7},{g,0,7}],2]];
 nn=Length[cat];hGate=cat[[209]];tGate=DiagonalMatrix[{1,phases[[2]]}];
 pt=phaseIndex /@ Range[nn];reps=Join[{1,9},Flatten[Table[17+64 n+g,{n,0,6},{g,0,7}]]];
 orbits=NestList[rightPhase,#,7]& /@ reps;decode=ConstantArray[{0,0},nn];
 Do[Do[decode[[orbits[[i,k]]]]={i,k-1},{k,8}],{i,58}];
 assert[Sort[Flatten[orbits]]==Range[nn],"free right-phase orbit partition"];
 Do[assert[RootReduce[ConjugateTranspose[cat[[a]]] . cat[[a]]]==IdentityMatrix[2],"unitary catalogue"];
 assert[RootReduce[tGate . cat[[a]]-cat[[pt[[a]]]]]==ConstantArray[0,{2,2}],"exact left T"];
 zz=Tr[ConjugateTranspose[cat[[rightPhase[a]]]] . cat[[a]] . tGate];
 assert[RootReduce[Conjugate[zz] zz]==4,"projective right-phase action"];
 assert[pt[[rightPhase[a]]]==rightPhase[pt[[a]]],"commuting left and right phase actions"],{a,nn}];
 conj=Conjugate[Flatten[#]]& /@ cat;threshold=(2-(13/50)^2)^2;qedges={};
 Do[target=Flatten[hGate . cat[[pt[[a]]]]];row={};
 Do[zz=conj[[j]] . target;ss=sgn[Conjugate[zz] zz-threshold];
 assert[MemberQ[{-1,0,1},ss],"resolved exact quotient edge comparison"];
 If[ss>=0,AppendTo[row,decode[[j]]]],{j,nn}];AppendTo[qedges,Sort[row]];
 If[Mod[Length[qedges],10]==0,Print["Exact quotient rows: ",Length[qedges]," / 58"]],{a,reps}];
 forced={2->27,34->2};forcedShift={2->5,34->4};
 Do[assert[MemberQ[qedges[[First[rule]]],{Last[rule],First[rule]/.forcedShift}],"allowed calibration edge and phase shift"],{rule,forced}];
 assert[decode[[12]]=={2,3} && decode[[209]]=={27,0} && decode[[216]]=={34,0} && decode[[13]]=={2,4},"quotient calibration derived from H anchors"];
 qgraph=restrictGraph[(Union[First /@ #]& /@ qedges),forced];
 assert[Total[Length /@ qgraph]==174,"complete constrained quotient graph has 174 edges"];
 initial=fullMatching[qgraph];assert[Sort[initial]==Range[58],"quotient graph has a full matching"];
 currentCase="full equivariant witness";block=witnessH[[pt]];
 assert[Sort[witnessH]==Range[nn] && witnessH[[13]]==209 && witnessH[[209]]==13,"full calibrated H permutation"];
 Do[assert[witnessH[[rightPhase[a]]]==rightPhase[witnessH[[a]]],"H phase equivariance"];
 zz=Tr[ConjugateTranspose[cat[[witnessH[[a]]]]] . hGate . cat[[a]]];
 assert[sgn[Conjugate[zz] zz-threshold]>=0,"every H edge satisfies the error certificate"],{a,nn}];
 qp=Table[decode[[block[[a]],1]],{a,reps}];shifts=Table[decode[[block[[a]],2]],{a,reps}];
 assert[Sort[qp]==Range[58] && And@@Table[MemberQ[qgraph[[a]],qp[[a]]] && MemberQ[qedges[[a]],{qp[[a]],shifts[[a]]}],{a,58}],"valid quotient witness"];
 cyc=cycle[qp,2];assert[cyc=={1,26,56,35,40,52,17,53,18,48,43,41,50,42,49,51,33}+1,"certified quotient seventeen-cycle"];
 assert[Length[cycle[block,13]]==17,"full frame period seventeen"];
 readq=Table[RootReduce[1-Conjugate[cat[[a,1,1]]] cat[[a,1,1]]],{a,reps}];
 Do[assert[RootReduce[1-Conjugate[cat[[a,1,1]]] cat[[a,1,1]]]==readq[[i]],"readout constant on phase orbits"],{i,58},{a,orbits[[i]]}];
 cc=-(2+Sqrt[2])/4;cosines={1,cc};Do[AppendTo[cosines,Expand[2 cc cosines[[-1]]-cosines[[-2]]]],{192}];
 ideal=Expand[(6-Sqrt[2])(1-#)/17]& /@ cosines;
 allowed=Table[sgn[(readq[[i]]-ideal[[n+1]])^2-1/9]<=0,{n,0,193},{i,58}];
 assert[And@@Table[allowed[[n+1,cyc[[1+Mod[n,17]]]]],{n,0,192}],"survival through block 192"];
 assert[!allowed[[194,cyc[[1+Mod[193,17]]]]],"first strict threshold failure at block 193"];
 assert[And@@Table[sgn[(readq[[cyc[[1+Mod[n,17]]]]]-ideal[[n+1]])^2-(13/40)^2]<0,{n,0,192}],"margin below 13/40 through 192"];
 assert[sgn[(readq[[cyc[[1+Mod[193,17]]]]]-ideal[[194]])^2-(17/50)^2]>0,"gap above 17/50 at 193"];
 states=Range[9,16];base=13;
 Do[Do[branch=cat[[states[[k]]]] . ConjugateTranspose[cat[[8+k]]];zz=Tr[ConjugateTranspose[cat[[base]]] . branch];
 assert[RootReduce[Conjugate[zz] zz]==4,"all compensated phase contexts agree projectively"],{k,8}];
 states=block[[states]];base=block[[base]],{n,0,193}];
 currentCase="exhaustive exact fixed-table history obstruction";
 Print["Beginning exhaustive exact 193-block search"];
 searchResult=historySearch[qgraph,allowed,2,forced,193];
 survives=searchResult["survives"];longest=searchResult["longest"];
 {nodes,newAssignments,bandReject,injectionReject,hallReject}=searchResult["counts"];
 assert[!survives && longest==192,"sharp maximal horizon: 192 feasible, 193 infeasible"];
 assert[{nodes,newAssignments,bandReject,injectionReject,hallReject}=={36975,12157,8038,11104,3897},"independent search reproduces exact traversal counts"];
 summary=<|"scope"->"Full phase-equivariant local C_8,8 matching class; all fixed phase priors; not the full two-qubit class.",
 "largest_feasible_horizon"->192,"first_unavoidable_curve_endpoint"->193,"D_curve"->194,
 "search_nodes"->nodes,"new_assignments"->newAssignments,"band_rejections"->bandReject,
 "injection_rejections"->injectionReject,"Hall_rejections"->hallReject,
 "witness_H_permutation"->witnessH-1,"witness_quotient_cycle"->cyc-1,
 "probabilities_n_0_to_193"->Table[ToString[readq[[cyc[[1+Mod[n,17]]]]],InputForm],{n,0,193}]|>;
 "passed"],"verification"],"aborted"];
If[!StringQ[status],status="failed"];If[!AssociationQ[summary],summary=<||>];
record=Join[metadata,<|"finished_utc"->utcString[],"status"->status,"last_case"->currentCase,"total_assertions"->Length[checkRecords],
 "passed_assertions"->Count[Lookup[checkRecords,"passed"],True],"results"->summary,"checks"->checkRecords|>];
If[recordReady,
 jsonPath=FileNameJoin[{recordDirectory,"probability_histories_run_record.json"}];
 textPath=FileNameJoin[{recordDirectory,"probability_histories_assertion_log.txt"}];
 jr=Quiet[Check[Export[jsonPath,record,"RawJSON"],$Failed]];
 tr=Quiet[Check[Export[textPath,StringRiffle[(If[TrueQ[#["passed"]],"PASS: ","FAIL: "]<>#["label"]& /@ checkRecords),"\n"],"Text"],$Failed]];
 exportStatus=If[StringQ[jr]&&StringQ[tr],"saved","failed"],exportStatus="failed"];
Print["Mathematical verification: ",status];Print["Record export: ",exportStatus];Print["Local record: ",recordDirectory];
If[!loadedFromFile,Print["Source hash unavailable. Load the saved .wl file using Get to capture it."]];
End[];EndPackage[];
