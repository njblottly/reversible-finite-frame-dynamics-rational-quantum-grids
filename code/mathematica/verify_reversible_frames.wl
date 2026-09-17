(* ::Package:: *)

(* Reversible finite-frame verification v0.2. Mathematica 13.3 compatible.
   Save and load using:
   With[{f=SystemDialogInput["FileOpen"]},If[StringQ[f],Get[f]]]
   Exact finite arithmetic and separate 50-digit numerical matching checks. *)
BeginPackage["RaQMReversibleVerification`"];
Begin["`Private`"];
Clear[matchingRows, historyRows, residuals];
currentCase="initialisation";matchingRegressionGraphs=0;
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
  "reversible_mathematica_" <> DateString[{"Year", "Month", "Day", "T", "Hour", "Minute", "Second"},
    TimeZone -> 0] <> "_" <> StringTake[CreateUUID[], 8] <> "Z"}];
directoryCreated = Quiet[Check[
  CreateDirectory[recordDirectory, CreateIntermediateDirectories -> True], $Failed]];
recordReady = StringQ[directoryCreated] && DirectoryQ[recordDirectory];
checkRecords = {}; textLog = {};
sourceHash = If[loadedFromFile,
  Quiet[Check[FileHash[scriptPath, "SHA256", "HexString"], Null]], Null];
metadata = <|"started_utc" -> utcString[], "script" -> "verify_reversible_frames.wl",
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
  If[!TrueQ[test], log["FAILED: " <> label]; log["Case: " <> ToString[currentCase,InputForm]]; Throw["failed", "verification"]]);


unitaryParam[t_,b_,g_] := {{Sqrt[t],Sqrt[1-t] Exp[I b]},
  {Sqrt[1-t] Exp[I g],-Sqrt[t] Exp[I(b+g)]}};
randomUnitary[] := N[Exp[I RandomInteger[{-1000,1000}]/1000]
  unitaryParam[RandomInteger[{1,999}]/1000,
  Pi RandomInteger[{-1000,1000}]/1000,Pi RandomInteger[{-1000,1000}]/1000],50];
unitaryQuantise[u_,d_,L_] := Module[{t,lam,b,g,n,bi,gi},
  t=Clip[Abs[u[[1,1]]]^2,{0,1}];
  Which[Abs[u[[1,1]]]<10^-35,
    t=0;lam=Arg[u[[1,2]]];b=0;g=Arg[u[[2,1]]]-lam,
    Abs[u[[1,2]]]<10^-35,
    t=1;lam=Arg[u[[1,1]]];b=0;g=Arg[-u[[2,2]]]-lam,
    True,lam=Arg[u[[1,1]]];b=Arg[u[[1,2]]]-lam;g=Arg[u[[2,1]]]-lam];
  n=Clip[Round[d t],{0,d}];bi=Mod[Round[L b/(2 Pi)],L];gi=Mod[Round[L g/(2 Pi)],L];
  {N[Exp[I lam] unitaryParam[n/d,2 Pi bi/L,2 Pi gi/L],50],n/d,{bi,gi}}];
unitaryBound[d_,L_] := N[1/Sqrt[d]+4 Sin[Pi/(2 L)],50];
stateDistance[a_,b_] := Module[{z=Conjugate[Flatten[a]] . Flatten[b]},
  Norm[a-If[Abs[z]==0,b,b Conjugate[z]/Abs[z]],"Frobenius"]];



catalogue[d_,l_] := Join[
 Table[unitaryParam[0,0,2 Pi k/l],{k,0,l-1}],
 Table[unitaryParam[1,0,2 Pi k/l],{k,0,l-1}],
 Flatten[Table[unitaryParam[n/d,2 Pi b/l,2 Pi g/l],
   {n,1,d-1},{b,0,l-1},{g,0,l-1}],2]];
pdistance[v_,w_] := Module[{z=Tr[ConjugateTranspose[v] . w],delta,gram},
 If[Abs[z]<10^-35,Return[N[Sqrt[2],50]]];
 delta=v-w Conjugate[z]/Abs[z];gram=ConjugateTranspose[delta] . delta;
 Sqrt[Max[0,(Re[gram[[1,1]]]+Re[gram[[2,2]]]+
   Sqrt[(Re[gram[[1,1]]]-Re[gram[[2,2]]])^2+4 Abs[gram[[1,2]]]^2])/2]]];
(* Explicit success flags avoid Return being caught by an enclosing loop. *)
validPermutation[p_,n_] := ListQ[p] && Length[p]==n && Sort[p]===Range[n];
perfectMatching[cost_,threshold_] := Module[
 {n=Length[cost],right,edges,seen,aug,p,i=1,ok=True},
 right=ConstantArray[0,n];
 edges=Table[Select[Range[n],TrueQ[cost[[row,#]]<=threshold]&],{row,n}];
 aug[vertex_] := Module[{j,pos=1,found=False},
   While[pos<=Length[edges[[vertex]]] && !found,
     j=edges[[vertex,pos]];
     If[!seen[[j]],seen[[j]]=True;
       If[right[[j]]==0 || TrueQ[aug[right[[j]]]],
         right[[j]]=vertex;found=True]];
     pos++];found];
 While[i<=n && ok,
   seen=ConstantArray[False,n];ok=TrueQ[aug[i]];i++];
 If[ok,
   p=ConstantArray[0,n];Do[p[[right[[j]]]]=j,{j,n}];p,
   $Failed]];
bottleneck[cost_] := Module[{values=Sort[DeleteDuplicates[Flatten[cost]]],lo=0,hi,mid},
 hi=Length[values];While[hi-lo>1,mid=Floor[(hi+lo)/2];
   If[perfectMatching[cost,values[[mid]]]===$Failed,lo=mid,hi=mid]];
 {perfectMatching[cost,values[[hi]]],values[[hi]],If[lo==0,Null,values[[lo]]]}];
cleanExtension[q_] := Module[{n=Length[q],groups,m,p,used,unusedIn,unusedOut},
 groups=Table[Flatten[Position[q,y]],{y,1,n}];m=Max[Length /@ groups];
 p=ConstantArray[0,n m];
 Do[Do[p[[(groups[[y,k]]-1)m+1]]=(y-1)m+k,{k,Length[groups[[y]]]}],{y,n}];
 unusedIn=Flatten[Position[p,0]];used=DeleteCases[p,0];
 unusedOut=Complement[Range[n m],used];
 Do[p[[unusedIn[[j]]]]=unusedOut[[j]],{j,Length[unusedIn]}];{m,p}];
status=Block[{$RecursionLimit=10000},CheckAbort[Catch[
 currentCase="exact matching regression";
 assert[perfectMatching[{{0,1},{1,0}},0]==={1,2},"exact regression: identity matching"];
 assert[perfectMatching[{{0,1},{0,1}},0]===$Failed,"exact regression: deficient graph rejected"];
 assert[perfectMatching[{{0,0},{0,1}},0]==={2,1},"exact regression: augmenting-path reassignment"];
 Do[Module[{mat=Partition[bits,3],candidate,exists},
   candidate=perfectMatching[mat,0];
   exists=Or@@Table[And@@Table[mat[[i,perm[[i]]]]==0,{i,3}],
     {perm,Permutations[Range[3]]}];
   assert[validPermutation[candidate,3]===exists,
     "exact regression: matching agrees with exhaustive permutation oracle"];
   If[exists,assert[And@@Table[mat[[i,candidate[[i]]]]==0,{i,3}],
     "exact regression: all selected edges belong to graph"]];
   matchingRegressionGraphs++],{bits,Tuples[{0,1},9]}];
 log["Exact matching regression passed on all 512 three-by-three graphs."];
 currentCase="finite-map and Haar-cell checks";
 mapCount=0;
 Do[Do[Module[{m,p}, {m,p}=cleanExtension[q];mapCount++;
   assert[Sort[p]==Range[n m],"exact clean-ancilla extension is a permutation"];
   assert[And@@Table[Quotient[p[[(x-1)m+1]]-1,m]+1==q[[x]],{x,n}],
     "exact clean-input rounding reproduced"];
   assert[m==Max[Values[Counts[q]]],"exact largest-fibre auxiliary alphabet"]],
   {q,Tuples[Range[n],n]}],{n,1,4}];
 Do[Module[{nn=2 l+(d-1)l^2,hh=d-1+2/l,cc,aa},cc=1/hh;aa=cc/l;
   assert[2 aa+(d-1)cc==1,"exact Haar slab volumes sum to one"];
   assert[cc/l^2==1/nn && aa/l==1/nn,"exact equal Haar cell volumes"];
   assert[And@@Table[aa+(n-1)cc<=n/d<=aa+n cc,{n,1,d-1}],
     "exact interior frames in assigned slabs"]],{l,{4,8,16,32,64,256}},{d,{2,l}}];
 gates={{{1,1},{1,-1}}/Sqrt[2],
   {{Cos[Pi/7],-Sin[Pi/7]},{Sin[Pi/7],Cos[Pi/7]}},unitaryParam[1/3,Pi/5,Pi/7]};
 names={"Hadamard","real rotation","complex gate"};
 matchingRows={};historyRows={};residuals={};tolerance=10^-25;
 Do[Module[{d=pair[[1]],l=pair[[2]],cat,nn,bound,tables={},cost,p,b,prev,pinv,u,
    ai,bi,start,actual,ideal,seed=N[{{1,1},{1,I}}/2,50],word={},side,gi,next,v,
    remote,previous,budget=0,finalError,res},
   cat=N[catalogue[d,l],50];nn=Length[cat];
   bound=N[2(Sqrt[2/(d-1+2/l)]+4 Sin[Pi/(2 l)]),50];
   assert[nn==2 l+(d-1)l^2,"exact catalogue size"];
   Do[currentCase=<|"D"->d,"L"->l,"frames"->nn,"gate"->names[[g]]|>;
     log["Matching: D="<>ToString[d]<>" L="<>ToString[l]<>" gate="<>names[[g]]];
     u=N[gates[[g]],50];cost=Table[pdistance[cat[[j]],u . cat[[i]]],{i,nn},{j,nn}];
     assert[MatrixQ[cost,NumericQ] && Dimensions[cost]=={nn,nn},"numerical cost matrix well formed"];
     {p,b,prev}=bottleneck[cost];
     If[!validPermutation[p,nn],log["Returned matching: "<>ToString[p,InputForm]]];
     assert[validPermutation[p,nn],"exact matched label permutation"];
     pinv=Ordering[p];
     assert[And@@Table[pinv[[p[[i]]]]==i,{i,nn}],"exact inverse permutation"];
     assert[Max[Table[cost[[i,p[[i]]]],{i,nn}]]<=b+tolerance,"numerical matched edge bound"];
     assert[prev===Null || perfectMatching[cost,prev]===$Failed,
       "numerical preceding threshold has no perfect matching"];
     assert[b<=bound+tolerance,"numerical matching obeys analytic bound"];
     assert[Max[Table[pdistance[cat[[pinv[[j]]]],ConjugateTranspose[u] . cat[[j]]],{j,nn}]]<=b+tolerance,
       "numerical inverse command bound"];
     AppendTo[tables,{u,p,pinv,b}];AppendTo[matchingRows,{d,l,nn,names[[g]],N[b,12],N[bound,12]}],{g,3}];
   ai=l+l/2+1;bi=ai;start={ai,bi};actual=seed;ideal=seed;
   Do[gi=Mod[step-1,3]+1;side=Mod[step-1,2];{u,p,pinv,b}=tables[[gi]];
     previous=If[side==0,ConjugateTranspose[actual] . actual,actual . ConjugateTranspose[actual]];
     If[side==0,next=p[[ai]];v=cat[[next]] . ConjugateTranspose[cat[[ai]]];ai=next;
       actual=v . actual;ideal=u . ideal,
       next=p[[bi]];v=cat[[next]] . ConjugateTranspose[cat[[bi]]];bi=next;
       actual=actual . Transpose[v];ideal=ideal . Transpose[u]];
     AppendTo[word,{side,gi}];budget+=b;
     res=Norm[ConjugateTranspose[v] . v-IdentityMatrix[2],"Frobenius"];AppendTo[residuals,res];
     assert[res<tolerance,"numerical selected local matrix unitary"];
     remote=If[side==0,ConjugateTranspose[actual] . actual,actual . ConjugateTranspose[actual]];
     assert[Norm[remote-previous,"Frobenius"]<tolerance,"numerical remote marginal preserved"];
     assert[Norm[actual-cat[[ai]] . seed . Transpose[cat[[bi]]],"Frobenius"]<tolerance,
       "numerical labelled endpoint"];
     assert[stateDistance[actual,ideal]<=Min[Sqrt[2],budget]+tolerance,"numerical word error bound"],
     {step,24}];finalError=stateDistance[actual,ideal];
   Do[{side,gi}=item;pinv=tables[[gi,3]];
     If[side==0,next=pinv[[ai]];v=cat[[next]] . ConjugateTranspose[cat[[ai]]];ai=next;actual=v . actual,
       next=pinv[[bi]];v=cat[[next]] . ConjugateTranspose[cat[[bi]]];bi=next;actual=actual . Transpose[v]],
     {item,Reverse[word]}];
   assert[{ai,bi}==start,"exact full word reversal restores labels"];
   assert[stateDistance[actual,seed]<tolerance,"numerical full word reversal restores state"];
   AppendTo[historyRows,{d,l,24,N[finalError,12],N[stateDistance[actual,seed],12]}]],
   {pair,{{2,4},{4,4},{4,8}}}];
 log["Exact finite maps exhausted: "<>ToString[mapCount]];
 log["Matching rows: D, L, N, gate, numerical bottleneck, proved bound"];log[matchingRows];
 log["History rows: D, L, forward steps, forward error, return residual"];log[historyRows];
 log["All reversible-frame Mathematica checks passed."];
 "passed","verification"],"aborted"]];
AssociateTo[metadata,{"finished_utc"->utcString[],"status"->status,"checks"->checkRecords,
 "matching_regression_graphs"->matchingRegressionGraphs,"last_case"->currentCase}];
If[ValueQ[matchingRows],AssociateTo[metadata,"matching_rows"->matchingRows]];
If[ValueQ[mapCount],AssociateTo[metadata,"finite_maps_exhausted"->mapCount]];

If[ValueQ[historyRows],AssociateTo[metadata,"history_rows"->historyRows]];
If[ValueQ[residuals] && Length[residuals]>0,
  AssociateTo[metadata,"maximum_unitarity_residual"->N[Max[residuals],16]]];
jsonPath=FileNameJoin[{recordDirectory,"reversible_frames_run_record.json"}];
logPath=FileNameJoin[{recordDirectory,"reversible_frames_assertion_log.txt"}];
jsonExport=$Failed;logExport=$Failed;
If[recordReady,
  jsonExport=Quiet[Check[Export[jsonPath,metadata,"RawJSON"],$Failed]];
  logExport=Quiet[Check[Export[logPath,StringRiffle[textLog,"\n"],"Text"],$Failed]]];
exportSucceeded=StringQ[jsonExport] && StringQ[logExport] &&
  FileExistsQ[jsonPath] && FileExistsQ[logPath];
Print["Mathematical verification: ",status];
Print["Record export: ",If[exportSucceeded,"saved","FAILED"]];
If[exportSucceeded,Print["Local record: ",recordDirectory],
  Print["The complete local record was not saved. Intended directory: ",recordDirectory]];
If[!StringQ[sourceHash],Print["Source hash unavailable. Load the saved .wl file using Get to capture it."]];
End[];
EndPackage[];
