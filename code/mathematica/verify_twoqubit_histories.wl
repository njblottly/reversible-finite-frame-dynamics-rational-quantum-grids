(* ::Package:: *)

(* Standalone two-qubit phase-lift and whole-history verifier.
Load the saved file using Get for source hash and local JSON/text records.
With[{f=SystemDialogInput["FileOpen"]},If[StringQ[f],Get[f]]]
Exact integer, rational and algebraic arithmetic. The 184320 objects are retained labels, not distinct physical matrices.
*)
BeginPackage["RaQMTwoQubitHistoryVerification`"];
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
  "twoqubit_histories_mathematica_" <> DateString[{"Year", "Month", "Day", "T", "Hour", "Minute", "Second"},
    TimeZone -> 0] <> "_" <> StringTake[CreateUUID[], 8] <> "Z"}];
directoryCreated = Quiet[Check[
  CreateDirectory[recordDirectory, CreateIntermediateDirectories -> True], $Failed]];
recordReady = StringQ[directoryCreated] && DirectoryQ[recordDirectory];
checkRecords = {}; textLog = {};
sourceHash = If[loadedFromFile,
  Quiet[Check[FileHash[scriptPath, "SHA256", "HexString"], Null]], Null];
metadata = <|"started_utc" -> utcString[], "script" -> "verify_twoqubit_histories.wl",
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
  "scope" -> "Exact labelled phase lift and classified shared local-H family with fixed CNOT; no unrestricted matching or refinement optimality claim."|>;
log[x_] := (AppendTo[textLog, ToString[x, InputForm]]; Print[x]);
assert[test_, label_] := (AppendTo[checkRecords, <|"label" -> label, "passed" -> TrueQ[test]|>];
  If[!TrueQ[test], log["FAILED: " <> label]; log["Case: " <> ToString[currentCase,InputForm]]; Throw["failed", "verification"]]);


(* Compare exact algebraic entries after subtracting, not whole lists
   with differently represented radicals. No numerical tolerance is used. *)
exactEqual[x_,y_] := If[Dimensions[x]===Dimensions[y],
 And@@(TrueQ[RootReduce[#]==0]& /@ Flatten[{x-y}]),False];
compose[p_,q_] := Sign[q] p[[Abs[q]]];
key[p_] := FromDigits[p+16,33];
act[p_,r_] := Module[{v=ConstantArray[0,16]},v[[Abs[p]]]=Sign[p] r;v];
local[p_,side_] := Flatten[Table[If[side==1,Sign[p[[a]]](4(Abs[p[[a]]]-1)+b),Sign[p[[b]]](4(a-1)+Abs[p[[b]]])],{a,4},{b,4}]];
enumerate[generators_,identity_] := Module[{out=ConstantArray[0,{12000,Length[identity]}],seen=<||>,nn=1,head=1,x,y,ky},
 out[[1]]=identity;AssociateTo[seen,key[identity]->1];
 While[head<=nn,x=out[[head]];
 Do[y=compose[g,x];ky=key[y];If[!KeyExistsQ[seen,ky],nn++;If[nn>12000,Throw["group overflow","verification"]];out[[nn]]=y;AssociateTo[seen,ky->nn]],{g,generators}];head++];{Take[out,nn],seen}];
state[rules_] := Normal[SparseArray[Join[{1->1},rules],16]];
phaseCoefficients[rr_,a_,b_] := Module[{out=rr,co,si,xx,yy,rx,ry},
 Do[co=Cos[Pi pair[[2]]/4];si=Sin[Pi pair[[2]]/4];
 Do[{xx,yy}=If[pair[[1]]==1,{5+j,9+j},{4 j+2,4 j+3}];rx=out[[xx]];ry=out[[yy]];
 out[[xx]]=Expand[co rx-si ry];out[[yy]]=Expand[si rx+co ry],{j,0,3}],{pair,{{1,a},{2,b}}}];out];
computational[r_] := Flatten[Table[Expand[(1+(-1)^a r[[13]]+(-1)^b r[[4]]+(-1)^(a+b) r[[16]])/4],{a,0,1},{b,0,1}]];
status=CheckAbort[Catch[Block[{$IterationLimit=Infinity},
 currentCase="exact algebraic comparison regression checks";
 assert[exactEqual[RootReduce[{{(1+I)/Sqrt[2],0},{0,(1-I)/Sqrt[2]}}],
 {{(1+I)/Sqrt[2],0},{0,(1-I)/Sqrt[2]}}],"equivalent radical matrix forms compare exactly"];
 assert[!exactEqual[{{1,0},{0,1}},{{1,0},{0,-1}}],"unequal exact matrices are rejected"];
 currentCase="exact phase-closure obstruction";
 eye=IdentityMatrix[2];sx={{0,1},{1,0}};sy={{0,-I},{I,0}};sz=DiagonalMatrix[{1,-1}];
 h={{1,1},{1,-1}}/Sqrt[2];t=DiagonalMatrix[{1,(1+I)/Sqrt[2]}];s=t . t;
 cm={{1,0,0,0},{0,1,0,0},{0,0,0,1},{0,0,1,0}};
 paulis=Flatten[Table[KroneckerProduct[a,b],{a,{eye,sx,sy,sz}},{b,{eye,sx,sy,sz}}],1];nonclosure={};
 Do[ee=KroneckerProduct[If[a==0,eye,t],If[b==0,eye,t]];
 kk=ConjugateTranspose[ee] . cm . KroneckerProduct[eye,t];im=kk . paulis[[2]] . ConjugateTranspose[kk];
 coeff=RootReduce[Table[Tr[pp . im]/4,{pp,paulis}]];support=Flatten[Position[coeff,Except[0],{1},Heads->False]];
 assert[Length[support]==If[b==0,2,4],"C T_B absent from original phase-Clifford translate"];
 AppendTo[nonclosure,<|"a"->a,"b"->b,"Pauli_support"->support-1|>];
 rel=RootReduce[ConjugateTranspose[cm] . ee . cm . ConjugateTranspose[ee]];
 assert[exactEqual[rel,DiagonalMatrix[If[b==0,{1,1,1,1},{1,1,(1+I)/Sqrt[2],(1-I)/Sqrt[2]}]]],"exact unchanged CNOT relative spectrum"],{a,0,1},{b,0,1}];
 currentCase="complete calibrated local-H classification";
 hid={1,4,-3,2};sid={1,3,-2,4};lid=Range[4];lacts={lid};ldata=<|key[lid]->{eye,""}|>;head=1;
 While[head<=Length[lacts],gg=lacts[[head]];
 Do[yk=compose[item[[2]],gg];ky=key[yk];If[!KeyExistsQ[ldata,ky],
 AppendTo[lacts,yk];AssociateTo[ldata,ky->{RootReduce[item[[3]] . ldata[key[gg]][[1]]],ldata[key[gg]][[2]]<>item[[1]]}]],
 {item,{{"H",hid,h},{"S",sid,s}}}];head++];
 assert[Length[lacts]==24,"complete local Clifford group"];
 threshold=3/2+Sqrt[2];candidates={};
 Do[mat=ldata[key[gg]][[1]];ww=ldata[key[gg]][[2]];physical=t . mat . ConjugateTranspose[t];zz=Tr[ConjugateTranspose[h] . physical];
 ss=Sign[RootReduce[Conjugate[zz] zz-threshold]];assert[MemberQ[{-1,0,1},ss],"exact local-H edge comparison"];
 If[ss>=0,AppendTo[candidates,{gg,mat,ww}]],{gg,lacts}];
 assert[candidates[[All,3]]=={"H","SH","HSSS","SHSSS"},"exactly four calibrated odd-parity H choices"];
 expected={h,h . s,ConjugateTranspose[s] . h,ConjugateTranspose[s] . h . s};
 Do[zz=Tr[ConjugateTranspose[candidates[[k,2]]] . expected[[k]]];assert[RootReduce[Conjugate[zz] zz]==4,"closed-form candidate matrix"];
 zz=Tr[ConjugateTranspose[h] . t . candidates[[k,2]] . ConjugateTranspose[t]];
 assert[RootReduce[Conjugate[zz] zz]==threshold,"sharp Hadamard bound for each choice"],{k,4}];
 zz=Tr[ConjugateTranspose[t . h] . h . t];assert[RootReduce[Conjugate[zz] zz]==threshold,"colliding physical frame labels have distinct H successors"];
 currentCase="exact two-qubit Clifford base";
 names={"HA0","HB0","SA","SB","C"};
 mats={KroneckerProduct[h,eye],KroneckerProduct[eye,h],KroneckerProduct[s,eye],KroneckerProduct[eye,s],cm};
 actions=Table[Table[coeff=RootReduce[Table[Tr[paulis[[j]] . mats[[k]] . paulis[[i]] . ConjugateTranspose[mats[[k]]]]/4,{j,16}]];
 pos=Flatten[Position[coeff,Except[0],{1},Heads->False]];assert[Length[pos]==1&&MemberQ[{-1,1},coeff[[First[pos]]]],"exact signed Pauli generator image"];First[pos] coeff[[First[pos]]],{i,16}],{k,5}];
 {group,index}=enumerate[actions,Range[16]];nn=Length[group];labelCount=16 nn;
 assert[nn==11520&&labelCount==184320,"base group and retained label counts"];
 left0=Table[Lookup[index,key[compose[gg,v]]],{gg,actions},{v,group}];
 leftA=Table[Lookup[index,key[compose[local[candidates[[k,1]],1],v]]],{k,4},{v,group}];
 leftB=Table[Lookup[index,key[compose[local[candidates[[k,1]],2],v]]],{k,4},{v,group}];
 rightA=Table[Lookup[index,key[compose[v,actions[[3]]]]],{v,group}];
 rightB=Table[Lookup[index,key[compose[v,actions[[4]]]]],{v,group}];
 enc[a_,b_,g_,r_,u_] := ((2 a+b) nn+g-1)4+2 r+u+1;
 dec[i_] := Module[{q=Quotient[i-1,4],rs=Mod[i-1,4],ab},ab=Quotient[q,nn];{Quotient[ab,2],Mod[ab,2],Mod[q,nn]+1,Quotient[rs,2],Mod[rs,2]}];
 Print["Constructing shared tables on ",labelCount," retained labels"];
 rhoA=Flatten[Table[enc[a,b,If[r==0,g,rightA[[g]]],1-r,u],{a,0,1},{b,0,1},{g,nn},{r,0,1},{u,0,1}]];
 rhoB=Flatten[Table[enc[a,b,If[u==0,g,rightB[[g]]],r,1-u],{a,0,1},{b,0,1},{g,nn},{r,0,1},{u,0,1}]];
 ta=Flatten[Table[enc[1-a,b,If[a==0,g,left0[[3,g]]],r,u],{a,0,1},{b,0,1},{g,nn},{r,0,1},{u,0,1}]];
 tb=Flatten[Table[enc[a,1-b,If[b==0,g,left0[[4,g]]],r,u],{a,0,1},{b,0,1},{g,nn},{r,0,1},{u,0,1}]];
 pc=Flatten[Table[enc[a,b,left0[[5,g]],r,u],{a,0,1},{b,0,1},{g,nn},{r,0,1},{u,0,1}]];
 has=Table[Flatten[Table[enc[a,b,If[a==0,left0[[1,g]],leftA[[k,g]]],r,u],{a,0,1},{b,0,1},{g,nn},{r,0,1},{u,0,1}]],{k,4}];
 hbs=Table[Flatten[Table[enc[a,b,If[b==0,left0[[2,g]],leftB[[k,g]]],r,u],{a,0,1},{b,0,1},{g,nn},{r,0,1},{u,0,1}]],{k,4}];
 allTables=Join[{ta,tb,pc},has,hbs];
 Do[assert[Sort[pp]==Range[labelCount],"full lifted table is a permutation"];
 assert[pp[[rhoA]]==rhoA[[pp]]&&pp[[rhoB]]==rhoB[[pp]],"full table commutes with both right-phase actions"],{pp,allTables}];
 assert[rhoA[[rhoB]]==rhoB[[rhoA]],"right-phase generators commute"];
 Do[assert[Nest[pp[[#]]&,Range[labelCount],8]==Range[labelCount],"right-phase generator eighth power identity"],{pp,{rhoA,rhoB}}];
 Do[assert[pa[[pb]]==pb[[pa]],"all opposite-party local tables commute"],{pa,Join[{ta},has]},{pb,Join[{tb},hbs]}];
 seen=ConstantArray[False,labelCount];orbitCount=0;
 Do[If[!seen[[i]],oa=NestList[rhoA[[#]]&,i,7];os=Union[Flatten[NestList[rhoB[[#]]&,#,7]& /@ oa]];
 assert[Length[os]==64,"free 64-label phase orbit"];seen[[os]]=ConstantArray[True,Length[os]];orbitCount++],{i,labelCount}];
 assert[orbitCount==2880,"phase quotient label count"];
 origin=enc[0,0,1,0,0];contexts=Flatten[Table[Nest[rhoB[[#]]&,Nest[rhoA[[#]]&,origin,k],l],{k,0,7},{l,0,7}]];
 Do[assert[pp[[pp[[origin]]]]==origin,"Hadamard calibrated identity round trip"],{pp,Join[has,hbs]}];
 assert[pc[[pc[[origin]]]]==origin,"CNOT calibrated identity round trip"];
 Do[assert[pc[[pc[[pp[[origin]]]]]]==pp[[origin]],"CNOT calibration at H_A and C H_A"],{pp,has}];
 currentCase="exact entangling histories for all fixed phase priors";
 xp=-(2+Sqrt[2])/4;xm=-(2-Sqrt[2])/4;cp={1,xp};cn={1,xm};
 Do[AppendTo[cp,Expand[2 xp cp[[-1]]-cp[[-2]]]];AppendTo[cn,Expand[2 xm cn[[-1]]-cn[[-2]]]],{63}];
 ideal=Table[RootReduce[((1-cp[[n+1]])(6-Sqrt[2])+(1-cn[[n+1]])(6+Sqrt[2]))/34],{n,0,64}];
 nominal=cm . KroneckerProduct[h . t,eye];mat=IdentityMatrix[4];
 Do[assert[RootReduce[Conjugate[mat[[{3,4},1]]] . mat[[{3,4},1]]]==ideal[[n+1]],"direct ideal matrix probability matches formula"];
 mat=RootReduce[nominal . mat],{n,0,64}];
 zero=state[{4->1,13->1,16->1}];rows={};
 Do[ends=contexts;probs={};dists={};first=-1;
 Do[values={};joint={};
 Do[{a,b,g,r,u}=dec[ends[[j]]];k=Quotient[j-1,8];l=Mod[j-1,8];
 rr=phaseCoefficients[act[group[[g]],phaseCoefficients[phaseCoefficients[zero,-k,-l],r,u]],a,b];
 value=RootReduce[(1-rr[[13]])/2];assert[IntegerQ[value]||Head[value]===Rational,"exact rational entangling probability"];
 AppendTo[values,value];AppendTo[joint,computational[rr]],{j,64}];
 assert[(SameQ@@values)&&(SameQ@@joint),"all 64 contexts have identical joint probabilities"];
 AppendTo[probs,First[values]];AppendTo[dists,First[joint]];
 If[first<0&&Abs[First[values]-ideal[[n+1]]]>1/3,first=n];ends=pc[[has[[choice,ta[[ends]]]]]],{n,0,64}];
 assert[first=={7,6,7,7}[[choice]],"first strict one-third discrepancy for local choice"];
 assert[probs[[8]]==1/2&&ideal[[8]]==1/16,"common seven-block discrepancy is 7/16"];
 AppendTo[rows,<|"matrix_form"->{"H","H S","S^-1 H","S^-1 H S"}[[choice]],"first_gap_above_1_over_3"->first,
 "A1_probabilities"->(ToString[#,InputForm]& /@ probs),"joint_probabilities"->Map[ToString[#,InputForm]&,dists,{2}]|>],{choice,4}];
 assert[And@@Table[Abs[ToExpression[rows[[1]]["A1_probabilities"][[n+1]]]-ideal[[n+1]]]<=1/3,{n,0,6}],"one fixed table survives the whole six-block history"];
 end=pc[[has[[1,origin]]]];{a,b,g,r,u}=dec[end];rr=phaseCoefficients[act[group[[g]],phaseCoefficients[zero,r,u]],a,b];
 assert[rr==state[{6->1,11->-1,16->1}],"H_A followed by C prepares an exact Bell state"];
 summary=<|"base_frames"->46080,"retained_labels"->labelCount,"phase_quotient_labels"->orbitCount,"additional_label_bits"->2,
 "labels_need_not_be_distinct_matrices"->True,"nonclosure_witness"->nonclosure,"history_rows"->rows,"common_point_depth_above_one_third"->7,
 "D_curve"->8,"common_point_gap"->"7/16","scope"->"Classified inherited local-H family, original C fixed; not unrestricted matching or refinement optimality."|>;
 "passed"],"verification"],"aborted"];
If[!StringQ[status],status="failed"];If[!AssociationQ[summary],summary=<||>];
record=Join[metadata,<|"finished_utc"->utcString[],"status"->status,"last_case"->currentCase,"total_assertions"->Length[checkRecords],
 "passed_assertions"->Count[Lookup[checkRecords,"passed"],True],"results"->summary,"checks"->checkRecords|>];
If[recordReady,jp=FileNameJoin[{recordDirectory,"twoqubit_histories_run_record.json"}];tp=FileNameJoin[{recordDirectory,"twoqubit_histories_assertion_log.txt"}];
 jr=Quiet[Check[Export[jp,record,"RawJSON"],$Failed]];
 tr=Quiet[Check[Export[tp,StringRiffle[(If[TrueQ[#["passed"]],"PASS: ","FAIL: "]<>#["label"]& /@ checkRecords),"\n"],"Text"],$Failed]];
 exportStatus=If[StringQ[jr]&&StringQ[tr],"saved","failed"],exportStatus="failed"];
Print["Mathematical verification: ",status];Print["Record export: ",exportStatus];Print["Local record: ",recordDirectory];
If[!loadedFromFile,Print["Source hash unavailable. Load the saved .wl file using Get to capture it."]];
End[];EndPackage[];
