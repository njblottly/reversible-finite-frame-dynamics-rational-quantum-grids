(* ::Package:: *)

(* Certified non-Clifford two-qubit controller v0.1. Mathematica 13.3.
   Standalone: enumerates 11520 Clifford frames using exact signed Pauli actions.
   With[{f=SystemDialogInput["FileOpen"]},If[StringQ[f],Get[f]]]
   Saves JSON and a text log under reproducibility/ beside the saved script.
   All verification arithmetic is exact. *)
BeginPackage["RaQMTwoQubitModelVerification`"];
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
  "two_qubit_mathematica_" <> DateString[{"Year", "Month", "Day", "T", "Hour", "Minute", "Second"},
    TimeZone -> 0] <> "_" <> StringTake[CreateUUID[], 8] <> "Z"}];
directoryCreated = Quiet[Check[
  CreateDirectory[recordDirectory, CreateIntermediateDirectories -> True], $Failed]];
recordReady = StringQ[directoryCreated] && DirectoryQ[recordDirectory];
checkRecords = {}; textLog = {};
sourceHash = If[loadedFromFile,
  Quiet[Check[FileHash[scriptPath, "SHA256", "HexString"], Null]], Null];
metadata = <|"started_utc" -> utcString[], "script" -> "verify_two_qubit_model.wl",
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
  "scope" -> "Coarse 46080-frame non-Clifford table; exact algebraic certification, not a bottleneck optimum or a refining family."|>;
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

simp[x_] := Expand[x];
phase[r_,a_,b_] := Module[{v=r,co={1,1/Sqrt[2],0,-1/Sqrt[2],-1,-1/Sqrt[2],0,1/Sqrt[2]},si={0,1/Sqrt[2],1,1/Sqrt[2],0,-1/Sqrt[2],-1,-1/Sqrt[2]},cx,sx0,xx,yy,rx,ry},
 Do[cx=co[[Mod[If[side==1,a,b],8]+1]];sx0=si[[Mod[If[side==1,a,b],8]+1]];
  Do[{xx,yy}=If[side==1,{5+j,9+j},{4j+2,4j+3}];{rx,ry}=v[[{xx,yy}]];
   v[[xx]]=simp[cx rx-sx0 ry];v[[yy]]=simp[sx0 rx+cx ry],{j,0,3}],{side,2}];v];
compProbs[r_] := Flatten[Table[simp[(1+(-1)^a r[[13]]+(-1)^b r[[4]]+(-1)^(a+b) r[[16]])/4],{a,0,1},{b,0,1}]];
expectedHashes=<|"HA"->"3728a0f3677152b9705b4f2b880c4401b2187f44aae7c092dc6ab1717e1b3a7a","HB"->"367acd37d052590e206708557ffd5b3ccbe716cb9526056911c5df3b2323dffe","TA"->"f82e962c784534125eb116a9bbc8d5da242607cd2f7c6555e9f3eb84381852f2","TB"->"fea55f9c9b1dc459b4e4b0dedb9110226c1f9c6a7d7cd739c740e4018421cc24","C"->"0e008ad66cc235bf1e0a74037fd4109c40d3968c9446427cb2e5e70618fba19c"|>;
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

 currentCase="complete non-Clifford table construction";log["Constructing the 46080-frame table..."];
 {group,index}=enumerate[actions];nn=Length[group];assert[nn==11520,"exact Clifford base size"];
 base=AssociationThread[names,Table[Lookup[index,key[compose[p,#]]]& /@ group,{p,actions}]];
 enc[a_,b_,g_] := (2a+b)nn+g;
 dec[x_] := {Quotient[Quotient[x-1,nn],2],Mod[Quotient[x-1,nn],2],Mod[x-1,nn]+1};
 controls={"HA","HB","TA","TB","C"};tables=<||>;
 Do[tab=Flatten[Table[Switch[name,
  "HA",enc[a,b,base["HA"][[g]]],"HB",enc[a,b,base["HB"][[g]]],"C",enc[a,b,base["C"][[g]]],
  "TA",enc[1-a,b,If[a==0,g,base["SA"][[g]]]],"TB",enc[a,1-b,If[b==0,g,base["SB"][[g]]]]],
  {a,0,1},{b,0,1},{g,nn}]];AssociateTo[tables,name->tab],{name,controls}];
 hashes=<||>;
 Do[p=tables[name];assert[Sort[p]==Range[4nn],name<>" full frame permutation"];
  ip=Ordering[p];assert[ip[[p]]==Range[4nn],name<>" inverse restores every label"];
  If[MemberQ[{"HA","HB","C"},name],assert[p[[p]]==Range[4nn],name<>" exact involution"],
   assert[Nest[p[[#]]&,Range[4nn],8]==Range[4nn],name<>" eighth power is identity"]];
  If[recordReady,
   hashPath=FileNameJoin[{recordDirectory,name<>"_canonical.csv"}];stream=OpenWrite[hashPath,BinaryFormat->True];
   BinaryWrite[stream,ToCharacterCode[StringRiffle[ToString /@ (p-1),","],"ASCII"],"Byte"];Close[stream];
   digest=FileHash[hashPath,"SHA256","HexString"];AssociateTo[hashes,name->digest];
   assert[digest==expectedHashes[name],name<>" frozen table hash matches Python"]],{name,controls}];
 assert[tables["HA"][[tables["HB"]]]==tables["HB"][[tables["HA"]]],"opposite Hadamards commute on all labels"];
 assert[tables["TA"][[tables["TB"]]]==tables["TB"][[tables["TA"]]],"opposite phase gates commute on all labels"];
 currentCase="sharp exact control certificates";
 t=DiagonalMatrix[{1,(1+I)/Sqrt[2]}];assert[simp[t . t]==s,"T carry equals exact S"];
 actual=simp[t . h . ConjugateTranspose[t]];rel=simp[ConjugateTranspose[h] . actual];delta=eye-rel;
 assert[simp[ConjugateTranspose[delta] . delta]==(1-1/Sqrt[2])eye,"sharp H squared-error certificate"];
 Do[e=KroneckerProduct[If[a==0,eye,t],If[b==0,eye,t]];
  actual=simp[e . cm . ConjugateTranspose[e]];rel=simp[ConjugateTranspose[cm] . actual];delta=IdentityMatrix[4]-rel;
  expected=DiagonalMatrix[{0,0,If[b==1,2-Sqrt[2],0],If[b==1,2-Sqrt[2],0]}];
  assert[simp[ConjugateTranspose[delta] . delta]==expected,"sharp CNOT squared-error spectrum"];
  assert[simp[actual . actual]==IdentityMatrix[4],"selected CNOT unitary is involutive"],{a,0,1},{b,0,1}];
 Do[If[{a,b}!={aa,bb},v=ConstantArray[0,16];v[[If[a!=aa,5,2]]]=1;
  assert[Count[phase[v,aa-a,bb-b],Except[0]]==2,"phase translates of Clifford group are distinct"]],
  {a,0,1},{b,0,1},{aa,0,1},{bb,0,1}];
 initial[k_,l_] := enc[Mod[k,2],Mod[l,2],Lookup[index,key[compose[pow[sa,Quotient[k,2]],pow[sb,Quotient[l,2]]]]]];
 finish[x_,word_] := Fold[tables[#2][[#1]]&,x,word];
 frameAct[x_,r_] := Module[{z=dec[x]},phase[act[group[[z[[3]]]],r],z[[1]],z[[2]]]];
 direct[word_,k_,l_,r_] := frameAct[finish[initial[k,l],word],phase[r,-k,-l]];
 physical[word_,k_,l_,input_] := Module[{x=initial[k,l],v=input,aa,bb,gg},
  Do[{aa,bb,gg}=dec[x];v=Switch[gate,
   "TA",phase[v,1,0],"TB",phase[v,0,1],
   "HA",phase[act[ha,phase[v,-aa,0]],aa,0],"HB",phase[act[hb,phase[v,0,-bb]],0,bb],
   "C",phase[act[c,phase[v,-aa,-bb]],aa,bb]];x=tables[gate][[x]],{gate,word}];v];
 zero=state[{4->1,13->1,16->1}];phi=state[{6->1,11->-1,16->1}];mixed=state[{}];
 pp=state[{2->1,5->1,6->1}];pi0=state[{4->1,9->1,12->1}];
 pa={id,local[{1,2,-3,-4},1],local[{1,-2,3,-4},1],local[{1,-2,-3,4},1]};
 bell=act[#,phi]& /@ pa;rotated=act[sa,#]& /@ bell;
 assert[direct[{"HA","C"},0,0,zero]==phi,"calibrated prefix creates a Bell state"];
 currentCase="exact preparation and channel verification";
 cases={{"HA","C","TB","C","HA"},{"TA","HB","C","HA","TB","C","TA","HB"},{"HB","TA","C","TB","HA","C"}};
 Do[Do[Do[got=direct[word,k,l,r];
   assert[got==physical[word,k,l,r],"retained-seed and physical-gate propagation agree"];
   assert[got==direct[word,Mod[k,2],Mod[l,2],r],"phase parities determine the channel"],{r,{zero,phi,pp,pi0}}];
  Do[avg=simp[mean[direct[word,k,l,#]& /@ recipe]];
   assert[avg==mixed,"equal-density recipes have equal outputs"],{recipe,{bell,rotated}}],{k,0,7},{l,0,7}],{word,cases}];
 currentCase="complete entangling phase scan";
 expectedCal={1,1,1/2,1/2,0,0,1/2,1/2,1};expectedUni={1,3/4,1/2,1/4,0,1/4,1/2,3/4,1};rows={};
 Do[word=Join[{"HA","C"},ConstantArray["TB",n],{"C","HA"}];probs={};
  Do[p=compProbs[direct[word,k,l,zero]];AppendTo[probs,p];m=Quotient[Mod[l,2]+n,2];ep=(1+{1,0,-1,0}[[Mod[m,4]+1]])/2;
   assert[p=={ep,0,1-ep,0},"exact phase-scan probability for every context"],{k,0,7},{l,0,7}];
  avg=simp[mean[probs]];assert[First[probs][[1]]==expectedCal[[n+1]] && avg[[1]]==expectedUni[[n+1]],"calibrated and uniform fringe row"];
  ideal=Fold[If[#2=="TB",phase[#1,0,1],act[gens[#2],#1]]&,zero,word];idealP=compProbs[ideal][[1]];
  assert[simp[idealP-(1+Cos[n Pi/4])/2]==0,"independent ideal interference formula"];
  If[OddQ[n],assert[simp[(idealP-avg[[1]])^2]==(3-2Sqrt[2])/16,"exact odd-phase gap squared"],
   assert[idealP==avg[[1]]==First[probs][[1]],"even phase values agree exactly"]];
  AppendTo[rows,<|"n"->n,"calibrated_p00"->ToString[expectedCal[[n+1]],InputForm],"uniform_p00"->ToString[expectedUni[[n+1]],InputForm],"ideal_p00"->ToString[idealP,InputForm]|>],{n,0,8}];
 currentCase="prior-independent contrast and finite context";
 Do[p1=compProbs[direct[{"HA","C","TB","C","HA"},k,l,zero]][[1]];
  p3=compProbs[direct[{"HA","C","TB","TB","TB","C","HA"},k,l,zero]][[1]];
  assert[p1-p3==1/2,"two-setting contrast independent of every initial context"],{k,0,7},{l,0,7}];
 contexts=Flatten[Table[direct[{"HA","HB"},a,b,zero],{a,0,1},{b,0,1}],1];
 assert[Length[DeleteDuplicates[contexts]]==4,"four phase contexts are future-distinguishable"];
 Do[weight=j/16;e1=(1-1/Sqrt[2]-weight)/2;e3=(1/Sqrt[2]-weight)/2;
  extreme=If[weight<=1/2,e3,-e1];
  assert[simp[extreme-((Sqrt[2]-1)/4+Abs[weight-1/2]/2)]==0,"sharp fixed-prior worst-error formula"],{j,0,16}];
 summary=<|"frames"->4nn,"clifford_base"->nn,"local_frames"->48,"nonlocal_charts"->20,"total_forward_table_edges"->5*4nn,
  "table_hashes"->hashes,"phase_scan"->rows,"effective_preparation_branches"->4,
  "prior_independent_p1_minus_p3"->"1/2","sharp_minimax_prior_error"->"(sqrt(2)-1)/4 + abs(w-1/2)/2",
  "scope"->"Coarse structured non-Clifford matching; no optimality, refinement or unique RaQM prediction claimed."|>;
 If[recordReady,tableExport=Quiet[Check[Export[FileNameJoin[{recordDirectory,"two_qubit_matching_table.json"}],
  <|"schema"->"raqm-two-qubit-structured-matching-v1","label_encoding"->"zero-based index=(2*a+b)*11520+g_index",
    "clifford_signed_pauli_actions"->group,"permutations"->AssociationThread[controls,(tables[#]-1)& /@ controls],
    "permutation_sha256_comma_joined_decimal"->hashes|>,"RawJSON"],$Failed]];
  log["Matching table export: "<>If[StringQ[tableExport],"saved","FAILED"]]];
 log[summary];log["All non-Clifford two-qubit checks passed."];"passed","verification"],"aborted"];
AssociateTo[metadata,{"finished_utc"->utcString[],"status"->status,"checks"->checkRecords,
 "passed_assertions"->Count[Lookup[checkRecords,"passed"],True],"total_assertions"->Length[checkRecords],"last_case"->currentCase}];
If[ValueQ[summary],AssociateTo[metadata,"summary"->summary]];
jsonPath=FileNameJoin[{recordDirectory,"two_qubit_run_record.json"}];logPath=FileNameJoin[{recordDirectory,"two_qubit_assertion_log.txt"}];
jsonExport=$Failed;logExport=$Failed;
If[recordReady,jsonExport=Quiet[Check[Export[jsonPath,metadata,"RawJSON"],$Failed]];
 logExport=Quiet[Check[Export[logPath,StringRiffle[textLog,"\n"]<>"\n"<>
 StringRiffle[Map[If[TrueQ[#["passed"]],"PASS: ","FAIL: "]<>#["label"]&,checkRecords],"\n"],"Text"],$Failed]]];
exportSucceeded=StringQ[jsonExport] && StringQ[logExport] && FileExistsQ[jsonPath] && FileExistsQ[logPath];
Print["Mathematical verification: ",status];Print["Record export: ",If[exportSucceeded,"saved","FAILED"]];Print["Local record: ",recordDirectory];
If[!StringQ[sourceHash],Print["Source hash unavailable. Load the saved .wl file using Get to capture it."]];
End[];
EndPackage[];
