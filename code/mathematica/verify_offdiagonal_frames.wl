(* ::Package:: *)

(* NEW off-diagonal and finite-frame verification, version 0.1.
   Not yet executed in Mathematica by the authoring environment.
   Save the file, then select it using this notebook cell:
   With[{f = SystemDialogInput["FileOpen"]}, If[StringQ[f], Get[f]]]
   The exact enumeration uses RootReduce and may take longer than previous
   checks depending on your hardware. 
   Fresh JSON and text records are written under reproducibility. *)

BeginPackage["RaQMOffdiagonalVerification`"];
Begin["`Private`"];
Clear[fibreRows, historyRows, residuals];
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
  "offdiagonal_mathematica_" <> DateString[{"Year", "Month", "Day", "T", "Hour", "Minute", "Second"},
    TimeZone -> 0] <> "_" <> StringTake[CreateUUID[], 8] <> "Z"}];
directoryCreated = Quiet[Check[
  CreateDirectory[recordDirectory, CreateIntermediateDirectories -> True], $Failed]];
recordReady = StringQ[directoryCreated] && DirectoryQ[recordDirectory];
checkRecords = {}; textLog = {};
sourceHash = If[loadedFromFile,
  Quiet[Check[FileHash[scriptPath, "SHA256", "HexString"], Null]], Null];
metadata = <|"started_utc" -> utcString[], "script" -> "verify_offdiagonal_frames.wl",
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


exactFibre[L_,r_,k_,occupation_] := Module[{count=0,a,b,c,p,coh,zz=Exp[2 Pi I/L]},
  Do[a=m/L;b=n/L;c=(1/2-a b)/(1-a);
    If[0<=c<=1 && IntegerQ[L c],
      p={a b,a(1-b),(1-a)c,(1-a)(1-c)};
      If[!occupation || AllTrue[L p,IntegerQ],
        Do[coh=RootReduce[a Sqrt[b(1-b)] zz^(-u)+(1-a)Sqrt[c(1-c)] zz^(-v)-r zz^k/2];
          If[coh===0,count+=L],{u,0,L-1},{v,0,L-1}]]],
    {m,1,L-1},{n,0,L}];count];

status=CheckAbort[Catch[
  fibreRows={};
  Do[Module[{count=exactFibre[L,r,k,kind=="K"],expected},
    expected=If[Mod[L,If[kind=="T",2,4] Denominator[r]]==0,2 L,0];
    assert[count==expected,"exact off-diagonal row/phase enumeration L="<>ToString[L]<>
      " r="<>ToString[r]<>" k="<>ToString[k]<>" "<>kind];
    AppendTo[fibreRows,{L,ToString[r,InputForm],k,kind,count}]],
    {L,{4,8}},{r,{1/2,1/4}},{k,{0,1}},{kind,{"T","K"}}];
  log["Exact fibre rows: L, r, k, grid, ray count"];log[fibreRows];
  Do[Module[{p=DeleteDuplicates[PowerMod[5,#,2^m]& /@ Range[0,2^(m-2)-1]]},
    assert[Length[p]==2^(m-2),"exact order of 5 modulo 2^M"];
    assert[Sort[Union[p,Mod[-p,2^m]]]==Range[1,2^m-1,2],"exact dyadic unit-group generators"]],
    {m,3,10}];
  Do[assert[Mod[n(2^s-n),4]==3,"exact nonhalf dyadic radicand congruence"],
    {s,2,8},{n,1,2^s-1,2}];

  h={{1,1},{1,-1}}/Sqrt[2];z=DiagonalMatrix[{1,-1}];
  nu={{1,1},{1,I}}/2;rho=nu . ConjugateTranspose[nu];
  assert[rho[[1,2]]==(1-I)/4,"exact irrational-Schmidt example coherence"];
  assert[FullSimplify[Det[rho]]==1/8,"exact irrational-Schmidt example determinant"];
  assert[FullSimplify[CharacteristicPolynomial[rho,x]-(x^2-x+1/8)]==0,
    "exact irrational-Schmidt characteristic polynomial"];
  assert[FullSimplify[(ConjugateTranspose[z-h] . (z-h))-(2-Sqrt[2])IdentityMatrix[2]]==
    ConstantArray[0,{2,2}],"exact sharp off-diagonal operator gap attained"];
  sq=DiagonalMatrix[{Sqrt[3]/2,1/2}];omega=sq . Transpose[h];target=h . omega;
  assert[FullSimplify[Abs[Conjugate[Flatten[z . omega]] . Flatten[target]]^2]==1/2,
    "exact sharp off-diagonal state gap attained"];
  assert[FullSimplify[Abs[Flatten[target]]^2-
    {(2+Sqrt[3])/8,(2-Sqrt[3])/8,(2-Sqrt[3])/8,(2+Sqrt[3])/8}]=={0,0,0,0},
    "exact finite-atlas endpoint outside original grid"];

  SeedRandom[20260909];tolerance=10^-25;historyRows={};residuals={};
  Do[Module[{aa=N[IdentityMatrix[2],50],bb=N[IdentityMatrix[2],50],
    seed=N[If[seedName=="irrational Schmidt",nu,omega],50],actual,ideal,
    bound=unitaryBound[L,L],u,v,next,t,phase,side,previous,remote,res,largest=0},
    actual=seed;ideal=seed;
    Do[u=randomUnitary[];side=If[OddQ[j],"A","B"];
      previous=If[side=="A",ConjugateTranspose[actual] . actual,actual . ConjugateTranspose[actual]];
      If[side=="A",
        {next,t,phase}=unitaryQuantise[u . aa,L,L];v=next . ConjugateTranspose[aa];aa=next;
        actual=v . actual;ideal=u . ideal,
        {next,t,phase}=unitaryQuantise[u . bb,L,L];v=next . ConjugateTranspose[bb];bb=next;
        actual=actual . Transpose[v];ideal=ideal . Transpose[u]];
      assert[Norm[v-u,2]<=bound+tolerance,"numerical finite-frame local step bound"];
      assert[IntegerQ[L t],"exact frame mixing certificate"];
      res=Norm[ConjugateTranspose[v] . v-IdentityMatrix[2],"Frobenius"];
      AppendTo[residuals,res];assert[res<tolerance,"numerical selected-gate unitarity"];
      assert[Norm[actual-aa . seed . Transpose[bb],"Frobenius"]<tolerance,
        "numerical finite-atlas endpoint reconstruction"];
      remote=If[side=="A",ConjugateTranspose[actual] . actual,actual . ConjugateTranspose[actual]];
      res=Norm[remote-previous,"Frobenius"];largest=Max[largest,res];
      assert[res<tolerance,"numerical unchanged remote marginal"];
      assert[stateDistance[actual,ideal]<=Min[Sqrt[2],j bound]+tolerance,
        "numerical finite-frame word bound"],{j,1,20}];
    AppendTo[historyRows,{L,seedName,20,N[largest,12],N[stateDistance[actual,ideal],12]}]],
    {L,{16,64,256}},{seedName,{"irrational Schmidt","off-diagonal seed"}}];
  log["Frame history rows: L, seed, steps, largest remote-marginal residual, final error"];log[historyRows];
  log["All off-diagonal/frame Mathematica checks passed."];
  "passed","verification"],"aborted"];
AssociateTo[metadata,{"finished_utc"->utcString[],"status"->status,"checks"->checkRecords}];
If[ValueQ[fibreRows],AssociateTo[metadata,"fibre_rows"->fibreRows]];

If[ValueQ[historyRows],AssociateTo[metadata,"history_rows"->historyRows]];
If[ValueQ[residuals] && Length[residuals]>0,
  AssociateTo[metadata,"maximum_unitarity_residual"->N[Max[residuals],16]]];
jsonPath=FileNameJoin[{recordDirectory,"offdiagonal_frames_run_record.json"}];
logPath=FileNameJoin[{recordDirectory,"offdiagonal_frames_assertion_log.txt"}];
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
