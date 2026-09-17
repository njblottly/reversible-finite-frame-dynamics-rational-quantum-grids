(* ::Package:: *)

(* Reachable-cycle certificates, standalone Wolfram Language verifier.
Load the saved file using Get to capture its hash and save records beside it.
Example: With[{f=SystemDialogInput["FileOpen"]},If[StringQ[f],Get[f]]]
Uses exact algebraic arithmetic.
*)
BeginPackage["RaQMReachableCycleVerification`"];
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
  "reachable_cycles_mathematica_" <> DateString[{"Year", "Month", "Day", "T", "Hour", "Minute", "Second"},
    TimeZone -> 0] <> "_" <> StringTake[CreateUUID[], 8] <> "Z"}];
directoryCreated = Quiet[Check[
  CreateDirectory[recordDirectory, CreateIntermediateDirectories -> True], $Failed]];
recordReady = StringQ[directoryCreated] && DirectoryQ[recordDirectory];
checkRecords = {}; textLog = {};
sourceHash = If[loadedFromFile,
  Quiet[Check[FileHash[scriptPath, "SHA256", "HexString"], Null]], Null];
metadata = <|"started_utc" -> utcString[], "script" -> "verify_reachable_cycles.wl",
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
  "scope" -> "Exact reachable-cycle restrictions, local matching graphs and onset certificates; no asymptotic onset claimed."|>;
log[x_] := (AppendTo[textLog, ToString[x, InputForm]]; Print[x]);
assert[test_, label_] := (AppendTo[checkRecords, <|"label" -> label, "passed" -> TrueQ[test]|>];
  If[!TrueQ[test], log["FAILED: " <> label]; log["Case: " <> ToString[currentCase,InputForm]]; Throw["failed", "verification"]]);

oldH={212,268,260,252,244,236,228,220,240,248,256,264,208,216,224,232,276,396,197,124,180,43,163,283,284,277,333,132,188,51,299,291,292,348,149,140,196,59,307,427,300,356,285,84,204,67,187,435,308,364,165,92,148,75,323,315,316,372,301,100,156,19,203,388,324,380,309,108,164,27,275,331,332,325,317,116,172,35,155,340,461,397,261,189,61,171,227,347,405,341,269,133,69,179,235,355,413,349,213,205,77,123,243,363,421,293,221,85,21,195,251,371,429,365,229,157,29,139,259,379,437,373,237,101,37,147,267,387,445,381,245,173,45,91,211,395,453,389,253,181,53,99,219,339,404,390,198,125,52,107,162,411,412,398,206,78,60,115,170,419,420,342,150,141,68,50,178,370,428,357,158,30,76,131,186,378,436,358,166,93,20,66,194,443,444,366,174,46,28,83,202,451,452,374,182,109,36,18,146,459,460,382,190,117,44,26,154,403,12,454,262,134,0,98,226,418,13,462,270,142,7,106,234,426,14,406,214,86,6,114,242,434,15,414,222,94,5,122,250,442,8,422,230,102,4,130,258,450,9,430,238,110,3,138,266,458,10,438,246,118,2,82,210,402,11,446,254,126,1,90,218,410,432,383,326,70,16,34,290,354,440,391,334,23,24,42,298,362,448,399,278,22,32,41,306,377,456,350,286,39,40,58,314,385,400,351,294,38,48,57,322,386,408,359,302,55,56,74,330,394,416,367,310,54,64,73,274,338,424,375,318,62,72,17,282,346,425,447,327,143,79,89,161,361,433,455,207,87,33,97,297,369,463,335,151,95,31,105,305,313,407,279,159,103,49,113,185,321,415,287,295,111,47,121,193,393,401,423,175,119,65,129,329,337,409,303,183,127,63,137,273,281,417,311,319,135,71,81,153,289,304,376,263,136,144,25,225,360,312,384,271,80,152,96,233,368,320,392,215,88,160,169,241,441,328,343,223,167,168,177,249,449,272,344,231,104,176,120,257,457,280,352,239,112,184,128,265,336,288,431,247,191,192,201,209,345,296,439,255,199,200,145,217,353}+1;
newH={212,268,260,252,244,236,228,220,240,248,256,264,208,216,224,232,276,396,261,124,180,108,227,283,284,340,269,132,188,51,235,356,292,348,149,140,196,123,243,299,300,421,285,84,204,67,187,371,308,364,229,29,148,139,259,315,316,372,237,100,156,83,267,323,324,380,181,45,164,92,211,395,332,388,189,116,172,99,155,339,461,397,325,125,43,107,291,411,405,277,333,197,69,115,171,419,413,349,213,205,77,59,179,427,435,293,221,85,21,131,251,307,429,301,165,157,75,203,195,443,437,309,173,37,19,147,331,451,445,317,245,109,27,91,275,459,453,389,253,117,35,163,219,403,404,390,326,61,52,34,290,347,412,341,334,133,60,42,298,355,420,342,278,141,68,50,306,363,428,357,286,30,76,58,314,378,436,365,294,93,20,66,322,379,444,373,302,101,28,74,330,387,452,381,310,54,36,18,274,338,460,382,318,53,44,26,282,346,12,454,262,134,0,98,226,418,13,398,270,142,7,106,234,426,14,406,214,86,6,114,242,434,15,414,222,94,5,122,250,442,8,422,230,102,4,130,258,450,9,430,238,110,3,138,266,458,10,438,246,118,2,82,210,402,11,446,254,126,1,90,218,410,432,383,198,70,16,89,162,354,440,462,206,78,24,97,170,362,448,399,150,22,32,105,178,370,456,350,158,103,40,113,186,385,400,358,166,38,48,57,194,386,408,366,174,46,56,129,202,394,416,374,182,63,64,137,146,345,424,375,190,62,72,81,154,417,425,447,199,207,79,25,161,361,455,391,335,151,23,33,297,433,441,463,215,95,41,169,241,377,407,343,223,167,49,177,249,449,415,351,295,175,47,121,257,457,423,359,239,119,65,193,265,337,409,431,183,127,73,201,209,281,439,311,319,135,17,145,217,353,304,376,263,143,144,88,225,360,312,327,271,87,152,96,233,369,320,392,279,31,160,104,305,313,328,336,159,39,168,112,185,384,272,287,231,111,176,120,321,393,280,352,303,55,184,128,329,401,288,367,247,191,192,136,273,344,296,368,255,71,200,80,153,289}+1;
unsupported={0,148}+1;
hallLeft={18,22,25,27,29,31,34,38,41,43,45,47,50,52,54,57,59,61,63,66,70,73,75,77,79,82,86,88,89,91,92,93,95,98,102,104,105,107,108,109,111,114,118,120,121,123,124,125,127,130,134,136,137,139,140,141,143,152,153,155,156,157,159,168,169,171,172,173,175,180,184,185,187,188,189,191,200,201,203,204,205,207,216,217,219,220,221,223,232,233,235,236,237,239,244,248,249,251,252,253,255,264,265,267,268,269,271,280,281,283,284,285,287,296,297,299,300,301,303,312,313,315,316,317,319,328,329,331,332,333,335,338,342,344,345,347,348,349,351,354,358,360,361,363,364,365,367,370,374,376,377,379,380,381,383,386,390,392,393,395,396,397,399,402,406,409,411,413,415,418,422,425,427,429,431,434,438,441,443,445,447,450,454,457,459,461,463}+1;
hallRight={1,3,4,5,7,9,11,13,15,17,19,20,21,23,24,26,28,30,33,35,37,39,40,42,44,46,49,51,53,55,56,58,60,62,65,67,69,71,72,74,76,78,80,81,83,84,85,87,90,94,96,97,99,100,101,103,106,110,112,113,115,116,117,119,122,126,128,129,131,132,133,135,138,142,145,147,149,151,161,163,165,167,177,179,181,183,193,195,197,199,209,211,213,215,225,227,229,231,241,243,245,247,257,259,261,263,273,275,277,279,289,291,293,295,305,307,309,311,321,323,325,327,336,337,339,340,341,343,346,350,352,353,355,356,357,359,362,366,368,369,371,372,373,375,378,382,384,385,387,388,389,391,394,398,401,403,405,407,408,410,412,414,417,419,421,423,424,426,428,430,433,435,437,439,440,442,444,446,449,451,453,455,456,458,460,462}+1;
sgn[x_] := Sign[RootReduce[x]];
simp[x_] := RootReduce[x];
phaseIndex[a_Integer] := 8 Quotient[a-1,8]+Mod[a,8]+1;
rightPhase[a_Integer] := Module[{i=a-1,n,b,g},
 Which[i<8,Mod[i-1,8]+1,i<16,9+Mod[i-8+1,8],True,
 n=Quotient[i-16,64]; b=Quotient[Mod[i-16,64],8];g=Mod[i-16,8];
 17+64 n+8 Mod[b+1,8]+g]];
restrictGraph[graph_,rules_List] := Module[{aa=First /@ rules,bb=Last /@ rules},
 Table[If[MemberQ[aa,i],{i /. rules},Complement[graph[[i]],bb]],{i,Length[graph]}]];
sccLabels[graph_] := Module[{nn=Length[graph],rev,seen,order={},visit,labels,label=0,stack,x,u},
 rev=Table[{}, {nn}];
 Do[Scan[Function[j,AppendTo[rev[[j]],i]],graph[[i]]],{i,nn}];
 seen=ConstantArray[False,nn];
 visit[v_] := (seen[[v]]=True;Scan[Function[y,If[!seen[[y]],visit[y]]],graph[[v]]];AppendTo[order,v]);
 Do[If[!seen[[i]],visit[i]],{i,nn}];labels=ConstantArray[0,nn];
 Do[If[labels[[x]]==0,label++;stack={x};labels[[x]]=label;
 While[Length[stack]>0,u=Last[stack];stack=Most[stack];
 Scan[Function[y,If[labels[[y]]==0,labels[[y]]=label;AppendTo[stack,y]]],rev[[u]]]]],{x,Reverse[order]}];labels];
cycle[p_,start_] := Module[{out={},x=start},While[!MemberQ[out,x],AppendTo[out,x];x=p[[x]]];
 assert[x==start,"cycle returns to its start"];out];
cosSequence[s_,count_] := Module[{c=-(2+s Sqrt[2])/4,out={1,-(2+s Sqrt[2])/4}},
 Do[AppendTo[out,Expand[2 c out[[-1]]-out[[-2]]]],{count-1}];out];
status=CheckAbort[Catch[Block[{$RecursionLimit=10000},
 currentCase="exact catalogue and phase symmetries";
 phases={1,(1+I)/Sqrt[2],I,(-1+I)/Sqrt[2],-1,(-1-I)/Sqrt[2],-I,(1-I)/Sqrt[2]};
 param[n_,b_,g_] := {{Sqrt[n/8],Sqrt[1-n/8] phases[[1+Mod[b,8]]]},
 {Sqrt[1-n/8] phases[[1+Mod[g,8]]],-Sqrt[n/8] phases[[1+Mod[b+g,8]]]}};
 cat=Join[Table[param[0,0,k],{k,0,7}],Table[param[8,0,k],{k,0,7}],
 Flatten[Table[param[n,b,g],{n,1,7},{b,0,7},{g,0,7}],2]];
 nn=Length[cat];hGate=cat[[209]];tGate=DiagonalMatrix[{1,phases[[2]]}];
 assert[nn==464 && cat[[13]]==IdentityMatrix[2],"catalogue and identity label"];
 pt=phaseIndex /@ Range[nn];
 Do[assert[simp[ConjugateTranspose[cat[[a]]] . cat[[a]]]==IdentityMatrix[2],"exact unitarity"];
 assert[simp[tGate . cat[[a]]-cat[[pt[[a]]]]]==ConstantArray[0,{2,2}],"exact left phase action"];
 zz=Tr[ConjugateTranspose[cat[[rightPhase[a]]]] . cat[[a]] . tGate];
 assert[simp[Conjugate[zz] zz]==4,"right phase symmetry up to global phase"],{a,nn}];
 currentCase="complete exact graph via right phase symmetry";
 hc=Flatten[hGate . #]& /@ cat;conj=Conjugate[Flatten[#]]& /@ cat;
 reps=Join[{1,9},Flatten[Table[17+64 n+g,{n,0,6},{g,0,7}]]];
 graph=Table[{}, {nn}];assigned={};threshold=(2-(13/50)^2)^2;
 Do[row={};Do[zz=conj[[j]] . hc[[rep]];edgeSign=sgn[Conjugate[zz] zz-threshold];
 assert[MemberQ[{-1,0,1},edgeSign],"resolved exact representative edge comparison"];
 If[edgeSign>=0,AppendTo[row,j]],{j,nn}];a=rep;
 Do[assert[!MemberQ[assigned,a],"right orbits partition catalogue"];AppendTo[assigned,a];
 graph[[a]]=Sort[row];a=rightPhase[a];row=rightPhase /@ row,{8}];
 If[Mod[Length[assigned],80]==0,Print["Exact graph labels: ",Length[assigned]," / ",nn]],{rep,reps}];
 graph=restrictGraph[graph,{13->209,209->13}];
 assert[Sort[oldH]==Range[nn] && And@@Table[MemberQ[graph[[a]],oldH[[a]]],{a,nn}],"existing full matching"];
 inv=Ordering[oldH];labels=sccLabels[Map[inv[[#]]&,graph,{2}]];
 supported=Table[Select[graph[[a]],labels[[a]]==labels[[inv[[#]]]]&],{a,nn}];
 assert[Total[Length /@ graph]==1420,"1420 calibrated geometric edges"];
 assert[Total[Length /@ supported]==1264,"1264 supported edges"];
 assert[MemberQ[graph[[unsupported[[1]]]],unsupported[[2]]] &&
 !MemberQ[supported[[unsupported[[1]]]],unsupported[[2]]],"unsupported geometric edge"];
 bad=restrictGraph[graph,{unsupported[[1]]->unsupported[[2]]}];
 assert[Union[Flatten[bad[[hallLeft]]]]==hallRight && Length[hallRight]<Length[hallLeft],"explicit Hall obstruction"];
 bg=graph[[pt]];sbg=supported[[pt]];
 assert[Length[Union[sccLabels[sbg]]]==1,"supported block graph has one SCC"];
 distances=ConstantArray[-1,nn];distances[[13]]=0;queue={13};short=0;
 While[Length[queue]>0 && short==0,a=First[queue];queue=Rest[queue];
 Do[If[j==13,short=distances[[a]]+1;Break[]];
 If[distances[[j]]<0,distances[[j]]=distances[[a]]+1;AppendTo[queue,j]],{j,sbg[[a]]}]];
 assert[short==5,"no shorter supported cycle through calibrated label"];
 assert[Sort[newH]==Range[nn] && newH[[13]]==209 && newH[[209]]==13 &&
 And@@Table[MemberQ[graph[[a]],newH[[a]]],{a,nn}],"full calibrated five-cycle matching"];
 oldB=oldH[[pt]];newB=newH[[pt]];path={12,216,398,353,215}+1;
 assert[cycle[newB,13]==path && Length[cycle[oldB,13]]==312,"exact cycles five and 312"];
 currentCase="exact spectral screens";cp=cosSequence[1,512];cm=cosSequence[-1,512];screens={};
 Do[den=item[[1]];single={};strong={};
 Do[e2=(q/den)^2;sp=(1-cp[[q+1]])/2;sm=(1-cm[[q+1]])/2;
 mu=2 (2+cp[[q+1]]+cm[[q+1]]);
 If[EvenQ[q],mu+=8 (-1)^(q/2) cp[[q/2+1]] cm[[q/2+1]]];
 assert[sgn[mu]>=0,"nonnegative exact trace squared"];
 pass=sgn[e2-sm]>=0;If[pass,AppendTo[single,q]];
 If[pass && sgn[e2-sp]>=0 && (e2>=2 || sgn[mu-(4-2 e2)^2]>=0),AppendTo[strong,q]],{q,1,256}];
 assert[{Length[single],First[single],Length[strong],First[strong]}==Rest[item],"exact spectral-screen counts"];
 AppendTo[screens,<|"error"->ToString[1/den,InputForm],"single_count"->Length[single],
 "single_first"->First[single],"strong_count"->Length[strong],"strong_first"->First[strong],"strong_survivors"->strong|>],
 {item,{{16,245,7,241,15},{64,216,11,192,22},{256,95,11,20,128},{1024,22,95,1,172}}}];
 currentCase="exact calibrated probability onset";onsetRows={};
 Do[b=item[[2]];a=13;first={-1,-1};probs={};
 Do[actual=simp[1-Conjugate[cat[[a,1,1]]] cat[[a,1,1]]];ideal=(6-Sqrt[2])(1-cp[[n+1]])/17;
 diff=actual-ideal;
 Do[If[first[[k]]<0 && sgn[diff^2-{1/20,1/3}[[k]]^2]>0,first[[k]]=n],{k,2}];
 AppendTo[probs,actual];a=b[[a]],{n,0,64}];
 assert[first==item[[3]],"first strict discrepancy thresholds: "<>item[[1]]];
 AppendTo[onsetRows,<|"table"->item[[1]],"first_above_1_over_20"->first[[1]],
 "first_above_1_over_3"->first[[2]],"probabilities_n_0_to_64"->(ToString[#,InputForm]& /@ probs)|>],
 {item,{{"frozen",oldB,{4,6}},{"five_cycle",newB,{2,19}}}}];
 summary=<|"N"->nn,"geometric_edges"->1420,"supported_edges"->1264,"block_scc_sizes"->{464},
 "spectral_screens"->screens,"onset_rows"->onsetRows,"five_cycle_H_permutation"->newH-1,
 "scope"->"Local (T,H) benchmark, calibrated preparation; no asymptotic or full two-qubit onset claim."|>;
 "passed"],"verification"],"aborted"];
If[!StringQ[status],status="failed"];
If[!AssociationQ[summary],summary=<||>];
record=Join[metadata,<|"finished_utc"->utcString[],"status"->status,
 "total_assertions"->Length[checkRecords],"passed_assertions"->Count[Lookup[checkRecords,"passed"],True],
 "results"->summary,"checks"->checkRecords|>];
If[recordReady,
 jsonPath=FileNameJoin[{recordDirectory,"reachable_cycles_run_record.json"}];
 textPath=FileNameJoin[{recordDirectory,"reachable_cycles_assertion_log.txt"}];
 jsonResult=Quiet[Check[Export[jsonPath,record,"RawJSON"],$Failed]];
 textResult=Quiet[Check[Export[textPath,StringRiffle[
 (If[TrueQ[#["passed"]],"PASS: ","FAIL: "]<>#["label"]& /@ checkRecords),"\n"],"Text"],$Failed]];
 exportStatus=If[StringQ[jsonResult]&&StringQ[textResult],"saved","failed"],exportStatus="failed"];
Print["Mathematical verification: ",status];Print["Record export: ",exportStatus];Print["Local record: ",recordDirectory];
If[!loadedFromFile,Print["Source hash unavailable. Load the saved .wl file using Get to capture it."]];
End[];EndPackage[];
