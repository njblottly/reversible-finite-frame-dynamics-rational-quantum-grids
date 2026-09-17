(* ::Package:: *)

(* Improved matchings and contrast robustness v0.1. Mathematica 13.3.
   Standalone: enumerates 11520 Clifford frames using exact signed Pauli actions.
   With[{f=SystemDialogInput["FileOpen"]},If[StringQ[f],Get[f]]]
   Saves JSON and a text log under reproducibility/ beside the saved script.
   All verification arithmetic is exact. *)
BeginPackage["RaQMContrastRefinementVerification`"];
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
  "contrast_refinement_mathematica_" <> DateString[{"Year", "Month", "Day", "T", "Hour", "Minute", "Second"},
    TimeZone -> 0] <> "_" <> StringTake[CreateUUID[], 8] <> "Z"}];
directoryCreated = Quiet[Check[
  CreateDirectory[recordDirectory, CreateIntermediateDirectories -> True], $Failed]];
recordReady = StringQ[directoryCreated] && DirectoryQ[recordDirectory];
checkRecords = {}; textLog = {};
sourceHash = If[loadedFromFile,
  Quiet[Check[FileHash[scriptPath, "SHA256", "HexString"], Null]], Null];
metadata = <|"started_utc" -> utcString[], "script" -> "verify_contrast_refinement.wl",
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
expectedHashes=<|"improved"-><|"HA"->"d28dd9591b8cbf82d8b7782e51d2300257824ec788d868dafcb63a55cdf70332","HB"->"b1a43923c6011d37c255968f0f67e87121116ac037d4f9e6d1aaf4be3f5e9298","TA"->"f82e962c784534125eb116a9bbc8d5da242607cd2f7c6555e9f3eb84381852f2","TB"->"fea55f9c9b1dc459b4e4b0dedb9110226c1f9c6a7d7cd739c740e4018421cc24","C"->"411bb85521b38374bf50992baf39a641f562b32ae5623937db51999fc44bd0f9"|>,"anchor_preserving"-><|"HA"->"3728a0f3677152b9705b4f2b880c4401b2187f44aae7c092dc6ab1717e1b3a7a","HB"->"367acd37d052590e206708557ffd5b3ccbe716cb9526056911c5df3b2323dffe","TA"->"f82e962c784534125eb116a9bbc8d5da242607cd2f7c6555e9f3eb84381852f2","TB"->"fea55f9c9b1dc459b4e4b0dedb9110226c1f9c6a7d7cd739c740e4018421cc24","C"->"d972aaba99bfef46dbcbc72d28f75c483154337d838219a9d182c26640d2ae67"|>|>;
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


 currentCase="same-catalogue matching variants";log["Constructing the matching variants..."];
 {group,index}=enumerate[actions];nn=Length[group];assert[nn==11520,"exact Clifford base size"];
 base=AssociationThread[names,Table[Lookup[index,key[compose[p,#]]]& /@ group,{p,actions}]];
 twisted=Lookup[index,key[compose[compose[inverse[sb],compose[c,sb]],#]]]& /@ group;
 enc[a_,b_,g_] := (2a+b)nn+g;
 dec[x_] := {Quotient[Quotient[x-1,nn],2],Mod[Quotient[x-1,nn],2],Mod[x-1,nn]+1};
 controls={"HA","HB","TA","TB","C"};original=<||>;improved=<||>;
 Do[tab=Flatten[Table[Switch[name,"HA",enc[a,b,base["HA"][[g]]],"HB",enc[a,b,base["HB"][[g]]],"C",enc[a,b,base["C"][[g]]],
   "TA",enc[1-a,b,If[a==0,g,base["SA"][[g]]]],"TB",enc[a,1-b,If[b==0,g,base["SB"][[g]]]]],{a,0,1},{b,0,1},{g,nn}]];
  AssociateTo[original,name->tab],{name,controls}];
 improved=original;
 AssociateTo[improved,{"HA"->Flatten[Table[enc[1-a,b,base["HA"][[g]]],{a,0,1},{b,0,1},{g,nn}]],
 "HB"->Flatten[Table[enc[a,1-b,base["HB"][[g]]],{a,0,1},{b,0,1},{g,nn}]],
 "C"->Flatten[Table[enc[a,1-b,base["C"][[g]]],{a,0,1},{b,0,1},{g,nn}]]}];
 anchored=original;AssociateTo[anchored,"C"->Flatten[Table[enc[a,b,If[b==0,base["C"][[g]],twisted[[g]]]],{a,0,1},{b,0,1},{g,nn}]]];
 variants=<|"improved"->improved,"anchor_preserving"->anchored|>;hashes=<||>;
 Do[tables=variants[model];digests=<||>;
  Do[p=tables[name];assert[Sort[p]==Range[4nn],model<>" "<>name<>" full permutation"];
   If[MemberQ[{"HA","HB","C"},name],assert[p[[p]]==Range[4nn],model<>" "<>name<>" involution"],
    assert[Nest[p[[#]]&,Range[4nn],8]==Range[4nn],model<>" "<>name<>" eighth power identity"]];
   If[recordReady,hp=FileNameJoin[{recordDirectory,model<>"_"<>name<>"_canonical.csv"}];stream=OpenWrite[hp,BinaryFormat->True];
    BinaryWrite[stream,ToCharacterCode[StringRiffle[ToString /@ (p-1),","],"ASCII"],"Byte"];Close[stream];digest=FileHash[hp,"SHA256","HexString"];
    AssociateTo[digests,name->digest];assert[digest==expectedHashes[model][name],model<>" "<>name<>" frozen hash"]],{name,controls}];
  Do[assert[tables[pair[[1]]][[tables[pair[[2]]]]]==tables[pair[[2]]][[tables[pair[[1]]]]],model<>" opposite local commands commute"],
    {pair,{{"HA","HB"},{"HA","TB"},{"TA","HB"},{"TA","TB"}}}];AssociateTo[hashes,model->digests],{model,Keys[variants]}];
 Do[x=enc[0,0,Lookup[index,key[g]]];assert[anchored["C"][[x]]==enc[0,0,Lookup[index,key[compose[c,g]]]],"exact C calibration retained"],{g,{id,c,ha,compose[c,ha]}}];
 currentCase="exact spectral error certificates";
 t=DiagonalMatrix[{1,(1+I)/Sqrt[2]}];zeta=(1+I)/Sqrt[2];eye4=IdentityMatrix[4];
 Do[actual=If[a==0,t . h,h . ConjugateTranspose[t]];rel=simp[ConjugateTranspose[h] . actual];z=If[a==0,zeta,Conjugate[zeta]];
  assert[simp[(rel-eye) . (rel-z eye)]==ConstantArray[0,{2,2}] && simp[Tr[rel]-1-z]==0,"improved H endpoint spectrum"],{a,0,1}];
 Do[before=KroneckerProduct[If[a==0,eye,t],If[b==0,eye,t]];after=KroneckerProduct[If[a==0,eye,t],If[b==1,eye,t]];
  actual=simp[after . cm . ConjugateTranspose[before]];rel=simp[ConjugateTranspose[cm] . actual];z=If[b==0,zeta,Conjugate[zeta]];
  assert[simp[(rel-eye4) . (rel-z eye4)]==ConstantArray[0,{4,4}] && simp[Tr[rel]-2-2z]==0,"improved C endpoint spectrum"];
  mid=If[b==0,cm,KroneckerProduct[eye,ConjugateTranspose[s]] . cm . KroneckerProduct[eye,s]];
  actual=simp[before . mid . ConjugateTranspose[before]];rel=simp[ConjugateTranspose[cm] . actual];delta=eye4-rel;
  assert[simp[ConjugateTranspose[delta] . delta]==DiagonalMatrix[{0,0,If[b==1,2-Sqrt[2],0],If[b==1,2-Sqrt[2],0]}],"anchored C error unchanged"],{a,0,1},{b,0,1}];
 assert[FullSimplify[2Sin[Pi/16]<Sqrt[1-1/Sqrt[2]]<Sqrt[2-Sqrt[2]]],"strict worst-error improvement"];
 initial[k_,l_] := enc[Mod[k,2],Mod[l,2],Lookup[index,key[compose[pow[sa,Quotient[k,2]],pow[sb,Quotient[l,2]]]]]];
 finish[x_,word_,ts_] := Fold[ts[#2][[#1]]&,x,word];
 direct[word_,k_,l_,r_,ts_] := Module[{z=dec[finish[initial[k,l],word,ts]]},phase[act[group[[z[[3]]]],phase[r,-k,-l]],z[[1]],z[[2]]]];
 zero=state[{4->1,13->1,16->1}];phi=state[{6->1,11->-1,16->1}];pp=state[{2->1,5->1,6->1}];
 currentCase="exact interference and preparation checks";rows={};
 Do[word=Join[{"HA","C"},ConstantArray["TB",n],{"C","HA"}];values={};
  Do[p=compProbs[direct[word,k,l,zero,improved]][[1]];pa=compProbs[direct[word,k,l,zero,anchored]][[1]];
   oldFlipped=compProbs[direct[word,k,BitXor[l,1],zero,original]][[1]];
   assert[p==pa==oldFlipped,"new fringes exchange original target-parity branches"];AppendTo[values,p],{k,0,7},{l,0,7}];
  AppendTo[rows,<|"n"->n,"calibrated_p00"->ToString[First[values],InputForm],"uniform_p00"->ToString[Mean[values],InputForm]|>],{n,0,8}];
 Do[tables=variants[model];
  Do[p1=compProbs[direct[{"HA","C","TB","C","HA"},k,l,zero,tables]][[1]];
   p3=compProbs[direct[{"HA","C","TB","TB","TB","C","HA"},k,l,zero,tables]][[1]];
   assert[p1-p3==1/2,model<>" prior-independent contrast"],{k,0,7},{l,0,7}];
  Do[Do[Do[assert[direct[word,k,l,r,tables]==direct[word,Mod[k,2],Mod[l,2],r,tables],model<>" phase-parity preparation law"],{r,{zero,phi,pp}}],
   {k,0,7},{l,0,7}],{word,{{"TA","HA","C","TB","HB","C"},{"HB","C","TA","HA","TB","C","HB"}}}],{model,Keys[variants]}];
 assert[1/2>(5/8)^2,"refinement threshold arithmetic; not an instantiated refined catalogue"];
 summary=<|"frames"->4nn,"improved_forward_edges"->5*4nn,"hashes"->hashes,"fringe_rows"->rows,
  "improved_sharp_H_and_C_error"->"2*sin(pi/16)","contrast_all_tested_models"->"1/2",
  "refinement_theorem"->"abs(contrast-1/sqrt(2)) <= 4*(epsH+epsC), uniformly over priors",
  "scope"->"Same-catalogue variants checked exactly. Genuine nested refinement is analytic, not enumerated here."|>;
 If[recordReady,Export[FileNameJoin[{recordDirectory,"contrast_matching_tables.json"}],
  <|"schema"->"raqm-contrast-matching-variants-v1","clifford_signed_pauli_actions"->group,
    "improved_permutations"->AssociationThread[controls,(improved[#]-1)& /@ controls],
    "anchor_preserving_C_permutation"->(anchored["C"]-1),"canonical_permutation_hashes"->hashes|>,"RawJSON"]];
 log[summary];log["All finite contrast-robustness checks passed."];"passed","verification"],"aborted"];
AssociateTo[metadata,{"finished_utc"->utcString[],"status"->status,"checks"->checkRecords,
 "passed_assertions"->Count[Lookup[checkRecords,"passed"],True],"total_assertions"->Length[checkRecords],"last_case"->currentCase}];
If[ValueQ[summary],AssociateTo[metadata,"summary"->summary]];
jsonPath=FileNameJoin[{recordDirectory,"contrast_refinement_run_record.json"}];logPath=FileNameJoin[{recordDirectory,"contrast_refinement_assertion_log.txt"}];
jsonExport=$Failed;logExport=$Failed;
If[recordReady,jsonExport=Quiet[Check[Export[jsonPath,metadata,"RawJSON"],$Failed]];
 logExport=Quiet[Check[Export[logPath,StringRiffle[textLog,"\n"]<>"\n"<>
 StringRiffle[Map[If[TrueQ[#["passed"]],"PASS: ","FAIL: "]<>#["label"]&,checkRecords],"\n"],"Text"],$Failed]]];
exportSucceeded=StringQ[jsonExport] && StringQ[logExport] && FileExistsQ[jsonPath] && FileExistsQ[logPath];
Print["Mathematical verification: ",status];Print["Record export: ",If[exportSucceeded,"saved","FAILED"]];Print["Local record: ",recordDirectory];
If[!StringQ[sourceHash],Print["Source hash unavailable. Load the saved .wl file using Get to capture it."]];
End[];
EndPackage[];
