(* ::Package:: *)

(* Entangled local-control checks, version 0.1.
   Save the file and load it using this notebook cell, selecting this .wl:
   With[{f = SystemDialogInput["FileOpen"]}, If[StringQ[f], Get[f]]]
   Fresh JSON and text records are created under reproducibility beside
   the script. Interactive evaluation falls back to the notebook directory
   and records the missing executed-source hash honestly. *)

Begin["RaQMEntangledLocalVerification`"];
Clear[fibreRows, coverageRows, historyRows, residuals];
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
  "entangled_mathematica_" <> DateString[{"Year", "Month", "Day", "T", "Hour", "Minute", "Second"},
    TimeZone -> 0] <> "_" <> StringTake[CreateUUID[], 8] <> "Z"}];
directoryCreated = Quiet[Check[
  CreateDirectory[recordDirectory, CreateIntermediateDirectories -> True], $Failed]];
recordReady = StringQ[directoryCreated] && DirectoryQ[recordDirectory];
checkRecords = {}; textLog = {};
sourceHash = If[loadedFromFile,
  Quiet[Check[FileHash[scriptPath, "SHA256", "HexString"], Null]], Null];
metadata = <|"started_utc" -> utcString[], "script" -> "verify_entangled_local.wl",
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


probabilityGrid[L_] := DeleteDuplicates[Flatten[Table[
  With[{a=m/L,b=n/L,c=r/L}, {a b,a(1-b),(1-a)c,(1-a)(1-c)}],
  {m,0,L},{n,0,L},{r,0,L}],2]];
fibreProbabilities[grid_,q_,side_] := If[side=="A",
  Select[grid, #[[1]]+#[[3]]==q && #[[1]] #[[2]]==#[[3]] #[[4]] &],
  Select[grid, #[[1]]+#[[2]]==q && #[[1]] #[[3]]==#[[2]] #[[4]] &]];
expectedMixing[L_,q_,side_,kind_] := Module[{d},
  If[side=="A" && q!=1/2,
    Return[If[kind=="T" || Mod[L,2 Denominator[q]]==0,{0,1/2,1},{0,1}]]];
  d=If[kind=="T",L,If[side=="A",L/2,L/Denominator[q]]];
  Range[0,d]/d];
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

status=CheckAbort[Catch[
  fibreRows={};
  Do[Module[{tg=probabilityGrid[L],grid,ps,ts,expected},
    Do[grid=If[kind=="T",tg,Select[tg,AllTrue[L #,IntegerQ]&]];
      Do[ps=fibreProbabilities[grid,q,side];
        ts=Sort[DeleteDuplicates[#[[1]]/q & /@ ps]];
        expected=expectedMixing[L,q,side,kind];
        assert[ts==expected,"exact full probability-grid fibre L="<>ToString[L]<>
          " q="<>ToString[q]<>" "<>side<>" "<>kind];
        If[MemberQ[{1/2,3/4},q],AppendTo[fibreRows,
          {L,ToString[q,InputForm],side,kind,Length[ts],2 L+(Length[ts]-2)L^2}]],
        {q,Range[1,L-1]/L},{side,{"A","B"}}],{kind,{"T","K"}}]],{L,{4,8,16}}];
  log["Fibre rows: L, q, active qubit, grid, mixing count, ray count"];log[fibreRows];

  theta=Pi/8;gap=2 Sin[Pi/16];
  ustar={{Cos[theta],-Sin[theta]},{Sin[theta],Cos[theta]}};
  assert[FullSimplify[ConjugateTranspose[IdentityMatrix[2]-ustar] .
    (IdentityMatrix[2]-ustar)-gap^2 IdentityMatrix[2]]==ConstantArray[0,{2,2}],
    "exact sharp operator gap attained by identity"];
  Do[assert[FullSimplify[Cos[theta] Sqrt[t]+Sin[theta] Sqrt[1-t]<=Cos[theta]],
    "exact column overlap bound"],{t,{0,1/2,1}}];

  SeedRandom[20260908];tolerance=10^-25;residuals={};coverageRows={};
  Do[Module[{bound=unitaryBound[d,L],targets,v,t,phases,err,largest=0,res},
    targets=Join[{N[IdentityMatrix[2],50],N[unitaryParam[0,0,1/3],50],
      N[unitaryParam[1,0,1/3],50]},Table[randomUnitary[],{20}]];
    Do[{v,t,phases}=unitaryQuantise[u,d,L];
      res=Norm[ConjugateTranspose[v] . v-IdentityMatrix[2],"Frobenius"];
      AppendTo[residuals,res];
      assert[res<tolerance,"numerical quantised unitary"];
      err=Norm[v-u,2];largest=Max[largest,err];
      assert[err<=bound+tolerance,"numerical unitary covering bound"];
      assert[IntegerQ[d t],"exact rounded mixing certificate"],{u,targets}];
    AppendTo[coverageRows,{L,d,N[largest,12],N[bound,12]}]],
    {L,{4,16,64,256}},{d,DeleteDuplicates[{1,L/4,L/2,L}]}];
  log["Unitary covering rows: L, D, largest sampled error, proved bound"];log[coverageRows];

  historyRows={};
  Do[Module[{d=If[kind=="T",L,L/2],bound,v=N[IdentityMatrix[2],50],
    actual=N[IdentityMatrix[2]/Sqrt[2],50],ideal=N[IdentityMatrix[2]/Sqrt[2],50],
    u,vp,t,phases,local,target,side,maximum=0,res},
    bound=unitaryBound[d,L];
    Do[u=randomUnitary[];side=If[OddQ[j],"A","B"];
      target=If[side=="A",u . v,v . Transpose[u]];
      {vp,t,phases}=unitaryQuantise[target,d,L];
      local=If[side=="A",vp . ConjugateTranspose[v],Transpose[ConjugateTranspose[v] . vp]];
      assert[Norm[local-u,2]<=bound+tolerance,"numerical Bell local step bound"];
      actual=If[side=="A",local . actual,actual . Transpose[local]];
      ideal=If[side=="A",u . ideal,ideal . Transpose[u]];v=vp;
      assert[Norm[actual-v/Sqrt[2],"Frobenius"]<tolerance,"numerical Bell action"];
      res=Max[Norm[actual . ConjugateTranspose[actual]-IdentityMatrix[2]/2,"Frobenius"],
        Norm[ConjugateTranspose[actual] . actual-IdentityMatrix[2]/2,"Frobenius"]];
      maximum=Max[maximum,res];assert[res<tolerance,"numerical Bell marginals"];
      assert[MemberQ[expectedMixing[L,1/2,"A",kind],t],"exact Bell endpoint certificate"];
      assert[stateDistance[actual,ideal]<=Min[Sqrt[2],j bound]+tolerance,
        "numerical Bell word bound"],{j,1,20}];
    AppendTo[historyRows,{L,kind,"Bell alternating",20,N[maximum,12]}]],
    {L,{16,64}},{kind,{"T","K"}}];

  Do[Module[{q=3/4,d=If[kind=="T",L,L/4],bound,v=N[IdentityMatrix[2],50],
    sq=N[DiagonalMatrix[{Sqrt[3]/2,1/2}],50],actual,u,vp,t,phases,local,maximum=0,res},
    actual=sq;bound=unitaryBound[d,L];
    Do[u=randomUnitary[];{vp,t,phases}=unitaryQuantise[u . v,d,L];
      local=vp . ConjugateTranspose[v];actual=actual . Transpose[local];v=vp;
      assert[Norm[local-u,2]<=bound+tolerance,"numerical B-only local step bound"];
      assert[Norm[actual-sq . Transpose[v],"Frobenius"]<tolerance,"numerical B-only action"];
      res=Norm[actual . ConjugateTranspose[actual]-DiagonalMatrix[{q,1-q}],"Frobenius"];
      maximum=Max[maximum,res];assert[res<tolerance,"numerical unchanged A marginal"];
      assert[MemberQ[expectedMixing[L,q,"B",kind],t],"exact B-only endpoint certificate"],{20}];
    AppendTo[historyRows,{L,kind,"q=3/4 B-only",20,N[maximum,12]}]],
    {L,{16,64}},{kind,{"T","K"}}];
  log["History rows: L, grid, sector, steps, largest marginal residual"];log[historyRows];
  log["All entangled-local Mathematica checks passed."];
  "passed","verification"],"aborted"];

AssociateTo[metadata,{"finished_utc"->utcString[],"status"->status,"checks"->checkRecords}];
If[ValueQ[fibreRows],AssociateTo[metadata,"fibre_rows"->fibreRows]];
If[ValueQ[coverageRows],AssociateTo[metadata,"coverage_rows"->coverageRows]];
If[ValueQ[historyRows],AssociateTo[metadata,"history_rows"->historyRows]];
If[ValueQ[residuals] && Length[residuals]>0,
  AssociateTo[metadata,"maximum_unitarity_residual"->N[Max[residuals],16]]];
jsonPath=FileNameJoin[{recordDirectory,"entangled_local_run_record.json"}];
logPath=FileNameJoin[{recordDirectory,"entangled_local_assertion_log.txt"}];
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
