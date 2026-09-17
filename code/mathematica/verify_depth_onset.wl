(* ::Package:: *)

(* Depth-onset certificates v0.1. Mathematica 13.3.
   Standalone: enumerates 11520 Clifford frames using exact signed Pauli actions.
   With[{f=SystemDialogInput["FileOpen"]},If[StringQ[f],Get[f]]]
   Saves JSON and a text log under reproducibility/ beside the saved script.
   All verification arithmetic is exact. *)
BeginPackage["RaQMDepthOnsetVerification`"];
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
  "depth_onset_mathematica_" <> DateString[{"Year", "Month", "Day", "T", "Hour", "Minute", "Second"},
    TimeZone -> 0] <> "_" <> StringTake[CreateUUID[], 8] <> "Z"}];
directoryCreated = Quiet[Check[
  CreateDirectory[recordDirectory, CreateIntermediateDirectories -> True], $Failed]];
recordReady = StringQ[directoryCreated] && DirectoryQ[recordDirectory];
checkRecords = {}; textLog = {};
sourceHash = If[loadedFromFile,
  Quiet[Check[FileHash[scriptPath, "SHA256", "HexString"], Null]], Null];
metadata = <|"started_utc" -> utcString[], "script" -> "verify_depth_onset.wl",
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
  "scope" -> "Exact spectral-budget certificates and earliest witnesses for three frozen tables. Universal finite-window bounds are analytic; no near-sqrt(L) sufficiency claim."|>;
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

 currentCase="exact spectral-budget certificates";
 budgetData={{8,27,7},{16,480,11},{32,480,22},{48,480,44},{64,481,55},{128,6970,128},
  {256,30220,139},{464,44293,406},{1024,20885827,684},{4096,20885827,3420},{46080,20885844,45828}};
 radicalSign[a_Integer,b_Integer] := Which[b==0,Sign[a],a==0,Sign[b],a>0 && b>0,1,a<0 && b<0,-1,True,Sign[a^2-2b^2] Sign[a]];
 assert[OrderedQ[budgetData[[All,2]]],"budget bounds are monotone"];
 aa=2;bb=0;a=-2;b=1;saved=<||>;budgetRows={};
 Do[If[q>1,{aa,bb,a,b}={a,b,-2a+2b-4aa,a-2b-4bb}];
  If[MemberQ[budgetData[[All,3]],q],AssociateTo[saved,q->{a,b}]];
  bd=SelectFirst[budgetData,#[[1]]>=q&];kk=bd[[2]];
  assert[radicalSign[kk^2(2^(q+1)-a)-q^2 2^(q+2),-kk^2 b]>=0,"spectral certificate at period "<>ToString[q]];
  If[MemberQ[budgetData[[All,1]],q],bd=SelectFirst[budgetData,#[[1]]==q&];{mb,kb,wb}=bd;{wa,wbb}=saved[wb];
   assert[radicalSign[(kb-1)^2(2^(wb+1)-wa)-wb^2 2^(wb+2),-(kb-1)^2 wbb]<0,"minimal integer bound at budget "<>ToString[q]];
   AppendTo[budgetRows,<|"frame_label_bound_M"->mb,"ceil_kappa_M"->kb,"minimality_witness_period"->wb,
    "window_length_D"->20(kb+1),"maximum_blocks"->20(kb+1)-1,"maximum_elementary_commands"->3(20(kb+1)-1)|>];
   Print["Certified spectral budget: ",q]],{q,1,46080}];
 assert[1<2<9/4 && 16<17<(17/4)^2,"rational enclosures for sqrt(2) and sqrt(17)"];
 assert[7+4>16(4/5)^2 && 7-17/4>16(2/5)^2,"cross-frequency sine lower bounds"];
 assert[7/8>(9/10)^2,"double-angle sine lower bound"];
 assert[8/17+15/56+10/81<1,"rational C0 remainder certificate"];
 assert[(1/5)/2-1/20==1/20,"universal contrast threshold arithmetic"];
 currentCase="earliest robust single-depth scan for the three tables";
 blocks=<||>;
 Do[tables=variants[model];p=tables["C"][[tables["HA"][[tables["TA"]]]]];
  assert[Sort[p]==Range[4nn],model<>" full mixed-block permutation"];AssociateTo[blocks,model->p],{model,Keys[variants]}];
 initial[k_,l_] := enc[Mod[k,2],Mod[l,2],Lookup[index,key[compose[pow[sa,Quotient[k,2]],pow[sb,Quotient[l,2]]]]]];
 starts=Flatten[Table[initial[k,l],{k,0,7},{l,0,7}]];
 positions=AssociationThread[Keys[variants],ConstantArray[starts,Length[variants]]];zero=state[{4->1,13->1,16->1}];
 xp=-(2+Sqrt[2])/4;xm=-(2-Sqrt[2])/4;cp0=1;cp1=xp;cn0=1;cn1=xm;
 firstSmall=Null;firstLarge=Null;rows={};
 Do[cp=If[depth==0,1,cp1];cn=If[depth==0,1,cn1];
  idealProb=simp[((1-cp)(6-Sqrt[2])+(1-cn)(6+Sqrt[2]))/34];modelRows=<||>;gaps=<||>;allValues=<||>;
  Do[values={};poss=positions[model];
   Do[{a,b,g}=dec[poss[[8k+l+1]]];rr=phase[act[group[[g]],phase[zero,-k,-l]],a,b];prob=simp[(1-rr[[13]])/2];
    assert[IntegerQ[prob] || Head[prob]===Rational,model<>" context probability is rational"];
    assert[0<=prob<=1,model<>" context probability is valid"];AppendTo[values,prob],{k,0,7},{l,0,7}];
   lo=Min[values];hi=Max[values];gap=Max[lo-idealProb,idealProb-hi,0];
   AssociateTo[gaps,model->gap];AssociateTo[allValues,model->values];
   AssociateTo[modelRows,model-><|"minimum_A1"->ToString[lo,InputForm],"maximum_A1"->ToString[hi,InputForm],
    "calibrated_A1"->ToString[First[values],InputForm],"uniform_A1"->ToString[Mean[values],InputForm],
    "minimum_prior_discrepancy"->ToString[gap,InputForm]|>],{model,Keys[variants]}];
  robust=Min[Values[gaps]];wmodel=SelectFirst[Keys[gaps],gaps[#]==robust&];values=allValues[wmodel];lo=Min[values];hi=Max[values];
  target=Min[Max[idealProb,lo],hi];weight=If[lo==hi,0,(target-lo)/(hi-lo)];
  assert[0<=weight<=1 && Abs[(1-weight)lo+weight hi-idealProb]==robust,"two-context prior attains the benchmark infimum"];
  If[firstSmall===Null && robust>1/20,firstSmall=depth];If[firstLarge===Null && robust>1/3,firstLarge=depth];
  AppendTo[rows,<|"blocks"->depth,"ideal_A1"->ToString[idealProb,InputForm],"minimum_over_three_tables_and_all_priors"->ToString[robust,InputForm],
   "models"->modelRows,"attaining_prior"-><|"model"->wmodel,"low_context_index"->First[FirstPosition[values,lo]]-1,
    "high_context_index"->First[FirstPosition[values,hi]]-1,"weight_on_high"->ToString[weight,InputForm],"context_index_encoding"->"8*k+ell; phases 0..7"|>|>];
  Do[AssociateTo[positions,model->blocks[model][[positions[model]]]],{model,Keys[variants]}];
  If[depth>=1,{cp0,cp1}={cp1,simp[2xp cp1-cp0]};{cn0,cn1}={cn1,simp[2xm cn1-cn0]}],{depth,0,24}];
 assert[firstSmall==5 && firstLarge==7,"earliest strict thresholds are five and seven blocks"];
 assert[rows[[6]]["ideal_A1"]=="3/8" && rows[[6]]["minimum_over_three_tables_and_all_priors"]=="1/8","five-block exact witness"];
 assert[rows[[8]]["ideal_A1"]=="1/16" && rows[[8]]["minimum_over_three_tables_and_all_priors"]=="7/16","seven-block exact witness"];
 resultTable=<|"schema"->"raqm-depth-onset-v1","spectral_budget_certificates"->budgetRows,"three_table_single_depth_scan"->rows,
  "earliest_three_table_strict_thresholds"-><|"1/20"->firstSmall,"1/3"->firstLarge|>,
  "scope"->"Earliest times concern only three named tables. Universal spectral result is a multi-depth contrast, not a single-depth optimum or near-sqrt(L) onset."|>;
 summary=<|"spectral_periods_certified"->46080,"largest_budget_certificate"->Last[budgetRows],
  "earliest_three_table_strict_thresholds"->resultTable["earliest_three_table_strict_thresholds"],"five_block_discrepancy"->"1/8",
  "seven_block_discrepancy"->"7/16","scope"->resultTable["scope"]|>;
 log[summary];log["All finite depth-onset checks passed."];"passed","verification"],"aborted"];
AssociateTo[metadata,{"finished_utc"->utcString[],"status"->status,"checks"->checkRecords,
 "passed_assertions"->Count[Lookup[checkRecords,"passed"],True],"total_assertions"->Length[checkRecords],"last_case"->currentCase}];
If[ValueQ[summary],AssociateTo[metadata,"summary"->summary]];
jsonPath=FileNameJoin[{recordDirectory,"depth_onset_run_record.json"}];logPath=FileNameJoin[{recordDirectory,"depth_onset_assertion_log.txt"}];
jsonExport=$Failed;logExport=$Failed;
If[recordReady,jsonExport=Quiet[Check[Export[jsonPath,metadata,"RawJSON"],$Failed]];
 logExport=Quiet[Check[Export[logPath,StringRiffle[textLog,"\n"]<>"\n"<>
 StringRiffle[Map[If[TrueQ[#["passed"]],"PASS: ","FAIL: "]<>#["label"]&,checkRecords],"\n"],"Text"],$Failed]]];
exportSucceeded=StringQ[jsonExport] && StringQ[logExport] && FileExistsQ[jsonPath] && FileExistsQ[logPath];
Print["Mathematical verification: ",status];Print["Record export: ",If[exportSucceeded,"saved","FAILED"]];Print["Local record: ",recordDirectory];
If[!StringQ[sourceHash],Print["Source hash unavailable. Load the saved .wl file using Get to capture it."]];
If[ValueQ[resultTable] && recordReady,
 tablePath=FileNameJoin[{recordDirectory,"depth_onset_results.json"}];
 tableExport=Quiet[Check[Export[tablePath,resultTable,"RawJSON"],$Failed]];
 Print["Result table export: ",If[StringQ[tableExport] && FileExistsQ[tablePath],"saved","FAILED"]]];
End[];
EndPackage[];
