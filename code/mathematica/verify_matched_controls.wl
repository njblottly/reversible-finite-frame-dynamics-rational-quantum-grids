(* ::Package:: *)

(* Matched H/T controls and compensated preparation v0.1. Mathematica 13.3.
   Standalone: the exact frozen permutation is embedded below.
   With[{f=SystemDialogInput["FileOpen"]},If[StringQ[f],Get[f]]]
   Saves JSON and a text log under reproducibility/ beside the saved script.
   All verification arithmetic is exact. *)
BeginPackage["RaQMMatchedControlVerification`"];
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
  "matched_controls_mathematica_" <> DateString[{"Year", "Month", "Day", "T", "Hour", "Minute", "Second"},
    TimeZone -> 0] <> "_" <> StringTake[CreateUUID[], 8] <> "Z"}];
directoryCreated = Quiet[Check[
  CreateDirectory[recordDirectory, CreateIntermediateDirectories -> True], $Failed]];
recordReady = StringQ[directoryCreated] && DirectoryQ[recordDirectory];
checkRecords = {}; textLog = {};
sourceHash = If[loadedFromFile,
  Quiet[Check[FileHash[scriptPath, "SHA256", "HexString"], Null]], Null];
metadata = <|"started_utc" -> utcString[], "script" -> "verify_matched_controls.wl",
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
  "scope" -> "Exact frozen-table, calibrated preparation and phase-compensated channel checks; no exact optimality or unique RaQM prediction claimed."|>;
log[x_] := (AppendTo[textLog, ToString[x, InputForm]]; Print[x]);
assert[test_, label_] := (AppendTo[checkRecords, <|"label" -> label, "passed" -> TrueQ[test]|>];
  If[!TrueQ[test], log["FAILED: " <> label]; log["Case: " <> ToString[currentCase,InputForm]]; Throw["failed", "verification"]]);


p0={212,268,260,252,244,236,228,220,240,248,256,264,208,216,224,232,276,396,197,124,180,43,163,283,284,277,333,132,188,51,299,291,292,348,149,140,196,59,307,427,300,356,285,84,204,67,187,435,308,364,165,92,148,75,323,315,316,372,301,100,156,19,203,388,324,380,309,108,164,27,275,331,332,325,317,116,172,35,155,340,461,397,261,189,61,171,227,347,405,341,269,133,69,179,235,355,413,349,213,205,77,123,243,363,421,293,221,85,21,195,251,371,429,365,229,157,29,139,259,379,437,373,237,101,37,147,267,387,445,381,245,173,45,91,211,395,453,389,253,181,53,99,219,339,404,390,198,125,52,107,162,411,412,398,206,78,60,115,170,419,420,342,150,141,68,50,178,370,428,357,158,30,76,131,186,378,436,358,166,93,20,66,194,443,444,366,174,46,28,83,202,451,452,374,182,109,36,18,146,459,460,382,190,117,44,26,154,403,12,454,262,134,0,98,226,418,13,462,270,142,7,106,234,426,14,406,214,86,6,114,242,434,15,414,222,94,5,122,250,442,8,422,230,102,4,130,258,450,9,430,238,110,3,138,266,458,10,438,246,118,2,82,210,402,11,446,254,126,1,90,218,410,432,383,326,70,16,34,290,354,440,391,334,23,24,42,298,362,448,399,278,22,32,41,306,377,456,350,286,39,40,58,314,385,400,351,294,38,48,57,322,386,408,359,302,55,56,74,330,394,416,367,310,54,64,73,274,338,424,375,318,62,72,17,282,346,425,447,327,143,79,89,161,361,433,455,207,87,33,97,297,369,463,335,151,95,31,105,305,313,407,279,159,103,49,113,185,321,415,287,295,111,47,121,193,393,401,423,175,119,65,129,329,337,409,303,183,127,63,137,273,281,417,311,319,135,71,81,153,289,304,376,263,136,144,25,225,360,312,384,271,80,152,96,233,368,320,392,215,88,160,169,241,441,328,343,223,167,168,177,249,449,272,344,231,104,176,120,257,457,280,352,239,112,184,128,265,336,288,431,247,191,192,201,209,345,296,439,255,199,200,145,217,353};
p=p0+1;
simp[x_] := FullSimplify[Expand[x]];
param[n_,b_,g_] := {{Sqrt[n/8],Sqrt[1-n/8] Exp[I Pi b/4]},
 {Sqrt[1-n/8] Exp[I Pi g/4],-Sqrt[n/8] Exp[I Pi (b+g)/4]}};
end[a_,word_] := Fold[Function[{current,gate},Switch[gate,"H",p[[current]],"h",ih[[current]],"T",pt[[current]],"t",it[[current]]]],a,Characters[word]];
outer[v_] := Outer[Times,v,Conjugate[v]];
status=CheckAbort[Catch[
 currentCase="exact frozen matching and error certificate";
 cat=simp[Join[Table[param[0,0,k],{k,0,7}],Table[param[8,0,k],{k,0,7}],
 Flatten[Table[param[n,b,g],{n,1,7},{b,0,7},{g,0,7}],2]]];
 nn=Length[cat];eye=IdentityMatrix[2];hGate={{1,1},{1,-1}}/Sqrt[2];tGate=DiagonalMatrix[{1,Exp[I Pi/4]}];
 pt=Table[8 Quotient[a-1,8]+Mod[a,8]+1,{a,nn}];ih=Ordering[p];it=Ordering[pt];
 assert[nn==464 && Sort[p]==Range[nn],"frozen H table is a permutation"];
 assert[Sort[pt]==Range[nn],"exact T table is a permutation"];
 assert[cat[[13]]==eye && cat[[209]]==hGate && p[[13]]==209 && p[[209]]==13,"calibrated H anchors"];
 threshold=(2-(13/50)^2)^2;
 Do[assert[simp[ConjugateTranspose[cat[[a]]] . cat[[a]]]==eye,"exact catalogue unitarity"];
   assert[simp[tGate . cat[[a]]-cat[[pt[[a]]]]]==ConstantArray[0,{2,2}],"exact T left multiplication"];
   zz=Tr[ConjugateTranspose[cat[[p[[a]]]]] . hGate . cat[[a]]];qq=simp[Conjugate[zz] zz];
   assert[FullSimplify[qq>=threshold],"exact algebraic certificate: H error at most 13/50"];
   assert[ih[[p[[a]]]]==a && it[[pt[[a]]]]==a,"exact inverse commands"],{a,nn}];
 words={"HTH","HTHTH","THHTHT","HTtHhT","HTHTHTHT"};
 Do[inv=StringJoin[Reverse[Characters[StringReplace[word,{"H"->"h","h"->"H","T"->"t","t"->"T"}]]]];
   Do[assert[end[end[a,word],inv]==a,"every frame restored by inverse word"],{a,nn}],{word,words}];
 currentCase="exact non-Clifford word probabilities";
 starts=Table[9+Mod[k+4,8],{k,0,7}];rows={};
 Do[word=item[[1]];expected=item[[2]];probs={};ends={};
   Do[a=starts[[k]];b=end[a,word];AppendTo[ends,b-1];vv=cat[[b]] . ConjugateTranspose[cat[[a]]];
     qq=simp[Conjugate[vv[[1,1]]] vv[[1,1]]];
     assert[qq==expected[[k]],"phase-conditioned output probability"];AppendTo[probs,qq],{k,8}];
   AppendTo[rows,<|"word"->word,"zero_based_final_indices"->ends,
     "phase_probabilities"->(ToString[#,InputForm]& /@ probs),"calibrated"->ToString[First[probs],InputForm],
     "uniform_phase_average"->ToString[Mean[probs],InputForm]|>],
   {item,{{"HTH",ConstantArray[7/8,8]},{"HTHTH",{3/4,3/4,3/4,3/4,7/8,7/8,7/8,3/4}}}}];
 ideal=simp[hGate . tGate . hGate];assert[simp[Abs[ideal[[1,1]]]^2]==(2+Sqrt[2])/4,"ideal HTH probability"];
 ideal=simp[hGate . tGate . hGate . tGate . hGate];assert[simp[Abs[ideal[[1,1]]]^2]==3/4,"ideal HTHTH probability"];
 currentCase="two-sided compensated preparation and equal-density recipes";
 sx={{0,1},{1,0}};sy={{0,-I},{I,0}};sz=DiagonalMatrix[{1,-1}];sGate=DiagonalMatrix[{1,I}];
 rs={eye,sx,sy,sz};ensembles={rs,sGate . #& /@ rs};target=IdentityMatrix[4]/4;
 Do[wa=pair[[1]];wb=pair[[2]];averages={ConstantArray[0,{4,4}],ConstantArray[0,{4,4}]};
   Do[ap=cat[[end[a,wa]]];bp=cat[[end[b,wb]]];va=ap . ConjugateTranspose[cat[[a]]];vb=bp . ConjugateTranspose[cat[[b]]];
     Do[density=ConstantArray[0,{4,4}];
       Do[psi=r/Sqrt[2];seed=ConjugateTranspose[cat[[a]]] . psi . Transpose[ConjugateTranspose[cat[[b]]]];
         actual=simp[ap . seed . Transpose[bp]];common=simp[va . psi . Transpose[vb]];
         assert[simp[actual-common]==ConstantArray[0,{2,2}],"compensated and common-channel endpoints agree"];
         density+=outer[Flatten[actual]]/4,{r,ensembles[[ei]]}];
       density=simp[density];assert[density==target,"each calibration preserves both I/4 recipes"];
       averages[[ei]]+=density/64,{ei,2}],{a,starts},{b,starts}];
   assert[simp[averages]=={target,target},"uniform phase law is two-sided mixture consistent"],
   {pair,{{"HTH","HTHTH"},{"HTHTH","THHTHT"}}}];
 summary=<|"D"->8,"L"->8,"N"->464,"certified_H_error_bound"->"13/50",
   "phase_control"->"T=diag(1,exp(i*pi/4)); exact","probability_rows"->rows,
   "arithmetic"->"exact algebraic arithmetic; no floating-point verification",
   "scope"->"Frozen candidate verified; all-word preparation consistency follows from the theorem"|>;
 log[summary];log["All matched-control preparation checks passed."];"passed","verification"],"aborted"];
AssociateTo[metadata,{"finished_utc"->utcString[],"status"->status,"checks"->checkRecords,
 "passed_assertions"->Count[Lookup[checkRecords,"passed"],True],"total_assertions"->Length[checkRecords],"last_case"->currentCase}];
If[ValueQ[summary],AssociateTo[metadata,"summary"->summary]];
jsonPath=FileNameJoin[{recordDirectory,"matched_controls_run_record.json"}];logPath=FileNameJoin[{recordDirectory,"matched_controls_assertion_log.txt"}];
jsonExport=$Failed;logExport=$Failed;
If[recordReady,jsonExport=Quiet[Check[Export[jsonPath,metadata,"RawJSON"],$Failed]];
 logExport=Quiet[Check[Export[logPath,StringRiffle[textLog,"\n"]<>"\n"<>
 StringRiffle[Map[If[TrueQ[#["passed"]],"PASS: ","FAIL: "]<>#["label"]&,checkRecords],"\n"],"Text"],$Failed]]];
exportSucceeded=StringQ[jsonExport] && StringQ[logExport] && FileExistsQ[jsonPath] && FileExistsQ[logPath];
Print["Mathematical verification: ",status];Print["Record export: ",If[exportSucceeded,"saved","FAILED"]];Print["Local record: ",recordDirectory];
If[!StringQ[sourceHash],Print["Source hash unavailable. Load the saved .wl file using Get to capture it."]];
End[];
EndPackage[];
