(* ::Package:: *)

(* Repeated mixed-control recurrence v0.1. Mathematica 13.3.
   Standalone: enumerates 11520 Clifford frames using exact signed Pauli actions.
   With[{f=SystemDialogInput["FileOpen"]},If[StringQ[f],Get[f]]]
   Saves JSON and a text log under reproducibility/ beside the saved script.
   All verification arithmetic is exact. *)
BeginPackage["RaQMMixedDepthVerification`"];
Begin["`Private`"];
Clear[summary, resultTable];
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
  "mixed_depth_mathematica_" <> DateString[{"Year", "Month", "Day", "T", "Hour", "Minute", "Second"},
    TimeZone -> 0] <> "_" <> StringTake[CreateUUID[], 8] <> "Z"}];
directoryCreated = Quiet[Check[
  CreateDirectory[recordDirectory, CreateIntermediateDirectories -> True], $Failed]];
recordReady = StringQ[directoryCreated] && DirectoryQ[recordDirectory];
checkRecords = {}; textLog = {};
sourceHash = If[loadedFromFile,
  Quiet[Check[FileHash[scriptPath, "SHA256", "HexString"], Null]], Null];
metadata = <|"started_utc" -> utcString[], "script" -> "verify_mixed_depth.wl",
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
  "scope" -> "Exact finite recurrence examples and algebraic ingredients. Universal recurrence and arbitrary-resolution statements are analytic, not numerically enumerated."|>;
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

 variants=<|"original"->original,"improved"->improved,"anchor_preserving"->anchored|>;
 permPower[p_,n_Integer] := Module[{ans=Range[Length[p]],b0=p,k0=n},
  While[k0>0,If[OddQ[k0],ans=b0[[ans]]];b0=b0[[b0]];k0=Quotient[k0,2]];ans];
 cycleTally[p_] := Module[{seen=ConstantArray[False,Length[p]],lens={},y,k0},
  Do[If[!seen[[j]],y=j;k0=0;While[!seen[[y]],seen[[y]]=True;k0++;y=p[[y]]];
   assert[y==j,"permutation orbit closes at its starting label"];AppendTo[lens,k0]],{j,Length[p]}];Sort[Tally[lens]]];
 blocks=<||>;cycles=<||>;orders=<||>;hashes=<||>;
 expectedCycles=<|"original"->{{16,2880}},"improved"->{{8,2880},{12,1920}},"anchor_preserving"->{{16,2880}}|>;
 Do[tables=variants[model];digests=<||>;
  Do[p=tables[name];assert[Sort[p]==Range[4nn],model<>" "<>name<>" full permutation"];
   If[recordReady && KeyExistsQ[expectedHashes,model],
    hp=FileNameJoin[{recordDirectory,model<>"_"<>name<>"_canonical.csv"}];stream=OpenWrite[hp,BinaryFormat->True];
    BinaryWrite[stream,ToCharacterCode[StringRiffle[ToString /@ (p-1),","],"ASCII"],"Byte"];Close[stream];
    digest=FileHash[hp,"SHA256","HexString"];AssociateTo[digests,name->digest];
    assert[digest==expectedHashes[model][name],model<>" "<>name<>" frozen hash"]],{name,controls}];
  AssociateTo[hashes,model->digests];
  p=tables["C"][[tables["HA"][[tables["TA"]]]]];
  assert[Sort[p]==Range[4nn],model<>" mixed block is a full permutation"];
  ct=cycleTally[p];ord=Apply[LCM,ct[[All,1]]];
  assert[ct==expectedCycles[model],model<>" exact cycle multiplicities"];
  assert[permPower[p,ord]==Range[4nn],model<>" all frames return at the block order"];
  assert[permPower[p,96]==Range[4nn],model<>" all frames return at 96 blocks"];
  AssociateTo[blocks,model->p];AssociateTo[cycles,model->ct];AssociateTo[orders,model->ord],{model,Keys[variants]}];
 currentCase="ideal two-qubit probability and irrational eigenphase ingredients";
 t=DiagonalMatrix[{1,(1+I)/Sqrt[2]}];vp=simp[h . t];vm=simp[sz . vp];nominal=simp[cm . KroneckerProduct[vp,eye]];
 xb=KroneckerProduct[eye,sx];assert[simp[nominal . xb-xb . nominal]==ConstantArray[0,{4,4}],"target X sectors are invariant"];
 Do[embedding=KroneckerProduct[eye,{{1},{sign}}/Sqrt[2]];
  assert[simp[ConjugateTranspose[embedding] . nominal . embedding]==If[sign==1,vp,vm],"exact target-X sector restriction"],{sign,{1,-1}}];
 assert[nominal[[All,1]]=={1,0,0,1}/Sqrt[2],"one block maps 00 to a Bell state"];
 taup=simp[Tr[vp] Conjugate[Tr[vp]]];taum=simp[Tr[vm] Conjugate[Tr[vm]]];
 assert[taup==1-1/Sqrt[2] && taum==1+1/Sqrt[2],"exact sector squared traces"];
 Do[assert[simp[2tau^2-4tau+1]==0,"primitive trace polynomial"],{tau,{taup,taum}}];
 assert[simp[taup taum]==1/2 && !IntegerQ[Sqrt[8]],"nonintegral norm and irreducible discriminant"];
 xp=-(2+Sqrt[2])/4;xm=-(2-Sqrt[2])/4;
 cp=ConstantArray[0,145];cn=ConstantArray[0,145];cp[[1]]=1;cn[[1]]=1;cp[[2]]=xp;cn[[2]]=xm;
 Do[cp[[k+1]]=simp[2xp cp[[k]]-cp[[k-1]]];cn[[k+1]]=simp[2xm cn[[k]]-cn[[k-1]]],{k,2,144}];
 ideal=ConstantArray[0,145];mat=IdentityMatrix[4];
 Do[prob=simp[((1-cp[[k+1]])(6-Sqrt[2])+(1-cn[[k+1]])(6+Sqrt[2]))/34];
  assert[IntegerQ[prob] || Head[prob]===Rational,"ideal probability is rational by conjugate cancellation"];
  assert[0<=prob<=1,"ideal probability is valid"];ideal[[k+1]]=prob;
  If[k<=96,direct=simp[Sum[Conjugate[mat[[i,1]]] mat[[i,1]],{i,3,4}]];
   assert[direct==prob,"4x4 matrix power agrees with Chebyshev formula"];mat=simp[nominal . mat]],{k,0,144}];
 assert[simp[1/(6+Sqrt[2])+1/(6-Sqrt[2])]==6/17,"Cesaro coefficient is 6/17"];
 assert[6/17-1/3==1/51 && 51(12/17)==36,"explicit finite-depth threshold coefficients"];
 Do[sp=(1-cp[[q+1]])/2;sm=(1-cn[[q+1]])/2;nm=simp[sp sm];scaled=simp[4*16^q nm];
  assert[nm>0 && (IntegerQ[nm] || Head[nm]===Rational),"small-Q sine-square norm is positive rational"];
  assert[IntegerQ[scaled] && scaled>=1,"norm denominator bound (finite examples)"],{q,1,16}];
 initial[k_,l_] := enc[Mod[k,2],Mod[l,2],Lookup[index,key[compose[pow[sa,Quotient[k,2]],pow[sb,Quotient[l,2]]]]]];
 zero=state[{4->1,13->1,16->1}];rows={};
 currentCase="calibrated, uniform and every-context return probabilities";
 Do[modelRows=<||>;
  Do[pn=permPower[blocks[model],count];values={};
   Do[{a,b,g}=dec[pn[[initial[k,l]]]];rr=phase[act[group[[g]],phase[zero,-k,-l]],a,b];prob=simp[(1-rr[[13]])/2];
    assert[IntegerQ[prob] || Head[prob]===Rational,model<>" context probability is rational"];
    assert[0<=prob<=1,model<>" context probability is valid"];AppendTo[values,prob],{k,0,7},{l,0,7}];
   AssociateTo[modelRows,model-><|"calibrated_A1"->ToString[First[values],InputForm],"uniform_A1"->ToString[Mean[values],InputForm],
    "minimum_context_A1"->ToString[Min[values],InputForm],"maximum_context_A1"->ToString[Max[values],InputForm]|>];
   If[Mod[count,orders[model]]==0,assert[AllTrue[values,#==0&],model<>" every phase context returns to zero A1"]],{model,Keys[variants]}];
  AppendTo[rows,<|"blocks"->count,"elementary_commands"->3count,"ideal_A1"->ToString[ideal[[count+1]],InputForm],"models"->modelRows|>],
  {count,{1,2,4,8,16,24,32,48,72,96,144}}];
 assert[ideal[[17]]==157/256,"original return discrepancy at 16 blocks"];
 assert[ideal[[25]]==1813/4096,"improved return discrepancy at 24 blocks"];
 assert[ideal[[97]]==124202613707725/281474976710656 && ideal[[97]]>1/3,"common 96-block discrepancy exceeds one third"];
 assert[ideal[[49]]==2421629/16777216 && ideal[[49]]<1/3,"not every return depth exceeds one third"];
 resultTable=<|"schema"->"raqm-mixed-depth-v1","chronological_block"->{"TA","HA","C"},"frames"->4nn,
  "cycle_counts_length_multiplicity"->cycles,"block_orders"->orders,"rows"->rows,
  "scope"->"Three coarse models checked. Universal recurrence and arbitrary-resolution results are analytic."|>;
 summary=<|"frames"->4nn,"cycle_counts_length_multiplicity"->cycles,"block_orders"->orders,
  "common_example_blocks"->96,"common_example_commands"->288,"finite_A1_at_common_example"->"0 in every context for all three tables",
  "ideal_A1_at_common_example"->ToString[ideal[[97]],InputForm],"universal_return_subsequence_mean"->"6/17",
  "uniform_positive_threshold"->"1/3","scope"->resultTable["scope"]|>;
 log[summary];log["All finite mixed-depth checks passed."];"passed","verification"],"aborted"];
AssociateTo[metadata,{"finished_utc"->utcString[],"status"->status,"checks"->checkRecords,
 "passed_assertions"->Count[Lookup[checkRecords,"passed"],True],"total_assertions"->Length[checkRecords],"last_case"->currentCase}];
If[ValueQ[summary],AssociateTo[metadata,"summary"->summary]];
jsonPath=FileNameJoin[{recordDirectory,"mixed_depth_run_record.json"}];logPath=FileNameJoin[{recordDirectory,"mixed_depth_assertion_log.txt"}];
jsonExport=$Failed;logExport=$Failed;
If[recordReady,jsonExport=Quiet[Check[Export[jsonPath,metadata,"RawJSON"],$Failed]];
 logExport=Quiet[Check[Export[logPath,StringRiffle[textLog,"\n"]<>"\n"<>
 StringRiffle[Map[If[TrueQ[#["passed"]],"PASS: ","FAIL: "]<>#["label"]&,checkRecords],"\n"],"Text"],$Failed]]];
exportSucceeded=StringQ[jsonExport] && StringQ[logExport] && FileExistsQ[jsonPath] && FileExistsQ[logPath];
Print["Mathematical verification: ",status];Print["Record export: ",If[exportSucceeded,"saved","FAILED"]];Print["Local record: ",recordDirectory];
If[!StringQ[sourceHash],Print["Source hash unavailable. Load the saved .wl file using Get to capture it."]];
If[ValueQ[resultTable] && recordReady,
 tablePath=FileNameJoin[{recordDirectory,"mixed_depth_results.json"}];
 tableExport=Quiet[Check[Export[tablePath,resultTable,"RawJSON"],$Failed]];
 Print["Result table export: ",If[StringQ[tableExport] && FileExistsQ[tablePath],"saved","FAILED"]]];
End[];
EndPackage[];
