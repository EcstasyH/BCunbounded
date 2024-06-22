#!wolframscript

myVerify [varSet_, flowVec_, initBoundary_, unsafeBoundary_, bc_, 
  deg_] := Module[
  {lambda, unsafeCons, initCons, bcLie, bc1, bc2, bc3, timev, range},
  (*Verify a BC candidate returned by the numerical solver satisfies \
the BC conditions*)
  lambda = -1;
  unsafeCons = And @@ (GreaterEqual[#, 0] & /@ unsafeBoundary);
  initCons = And @@ (GreaterEqual[#, 0] & /@ initBoundary);
  If[Length[varSet] == 2,
   bcLie = 
    D[bc, varSet[[1]]]*flowVec[[1]] + 
     D[bc, varSet[[2]]]*flowVec[[2]],
   bcLie = 
     D[bc, varSet[[1]]]*flowVec[[1]] + 
      D[bc, varSet[[2]]]*flowVec[[2]] + 
      D[bc, varSet[[3]]]*flowVec[[3]];
   ];
  
  bc1 = FindInstance[bc > 0 && initCons, varSet];
  bc2 = FindInstance[bc <= 0 && unsafeCons, varSet];
  bc3 = FindInstance[bc == 0 && bcLie > 0, varSet];
  If[Length[bc1] + Length[bc2] + Length[bc3] == 0,
   Return[1]];
  Return[0]
  ]

verifyResults[benchmark_, plot_, plotRange_, verifyTime_] := Module[
  {systemPath, varSet, flowVec, initCons, unsafeCons, bcPath, bcData, 
   ebc, hbc, ebcDeg, hbcDeg, totaltime, totalsdptime, time, etime, 
   htime, ptime, result, i, sbc, sbcDeg},
  (*Read the results from files and verify each results*)
  
  systemPath = 
   FileNameJoin[{".", "Results", "systems", 
     StringJoin[benchmark, ".txt"]}];
  {varSet, flowVec, initCons, unsafeCons} = ReadList[systemPath];
  ebcDeg = 0;
  hbcDeg = 0;
  Print["Variables: ", varSet];
  Print["Flow: ", flowVec];
  Print["Initial: ", initCons];
  Print["Unsafe: ", unsafeCons];
  
  (*Using the original SOS characterization from CAV20*)
  (*Print[Style[
  "Using CAV20 incomplete Characterization for Polynomial BC: ", 
  Blue]];*)
  Print["=== Thm. 3 (CAV20) Polynomial BC ==="];
  bcPath = 
   FileNameJoin[{".", "Results", "sound", 
     StringJoin[benchmark, ".txt"]}];
  bcData = 
   ReadList[bcPath, {Expression, Expression}, 
    RecordSeparators -> "\n"];
  totaltime = 0;
  totalsdptime = 0;
  For[i = 1, i <= Length[bcData], i++,
   ebc = bcData[[i]][[1]];
   time = bcData[[i]][[2]];
   result = 
    Timing[TimeConstrained[
      myVerify[varSet, flowVec, initCons, unsafeCons, ebc, i], 
      verifyTime, 2]];
   totaltime = totaltime + time + result[[1]];
   totalsdptime = totalsdptime + time;
   If[result[[2]] == 1,
    ebcDeg = i;
    Print["Verified at deg ", i];
    (*Print[ebc];*)
    Break[]
    ];
   If[result[[2]] == 2,
    Print["Verify out of time at deg ", i];
    ebcDeg = 0;
    (*Print[ebc];*)
    Break[]
    ];
   ];
  etime = totaltime;
  Print["SDP time: ", Round[totalsdptime, 0.01] ];
  If[ebcDeg > 0,
   Print["Verified, SDP+verify time: ", Round[etime, 0.01]];,
   Print["No valid soltion, SDP+verify time: ", Round[etime, 0.01]];
   ];
  
  (*Using our SOS characterization for polynomial BC*)
  (*Print[Style["Using Complete Characterization for Polynomial BC:", 
  Blue]];*)
  Print["=== Thm. 5 Polynomial BC ==="];
  bcPath = 
   FileNameJoin[{".", "Results", "complete", 
     StringJoin[benchmark, ".txt"]}];
  bcData = 
   ReadList[bcPath, {Expression, Expression}, 
    RecordSeparators -> "\n"];
  totaltime = 0;
  totalsdptime = 0;
  For[i = 1, i <= Length[bcData], i++,
   hbc = bcData[[i]][[1]];
   time = bcData[[i]][[2]];
   result = 
    Timing[TimeConstrained[
      myVerify[varSet, flowVec, initCons, unsafeCons, hbc, i], 
      verifyTime, 2]];
   totaltime = totaltime + time + result[[1]];
   totalsdptime = totalsdptime + time;
   If[result[[2]] == 1,
    hbcDeg = i;
    Print["Verified at deg ", i];
    (*Print[hbc];*)
    Break[]
    ];
   If[result[[2]] == 2,
    Print["Verify out of time at deg ", i];
    hbcDeg = 0;
    (*Print[hbc];*)
    Break[]
    ];
   ];
  htime = totaltime;
  Print["SDP time: ", Round[totalsdptime, 0.01] ];
  If[hbcDeg > 0,
   Print["Verified, SDP+verify time: ", Round[htime, 0.01]];,
   Print["No valid soltion, SDP+verify time: ", Round[htime, 0.01]];
   ];
  
  (*Using our SOS characterization for semialgebraic BC*)
  (*Print[Style[
  "Using Complete Characterization for Semialgebraic BC:", Blue]];*)
  Print["=== Thm. 7 Semialgebraic BC ==="];
  bcPath = 
   FileNameJoin[{".", "Results", "completesemi", 
     StringJoin[benchmark, ".txt"]}];
  bcData = 
   ReadList[bcPath, {Expression, Expression}, 
    RecordSeparators -> "\n"];
  totaltime = 0;
  totalsdptime = 0;
  For[i = 1, i <= Length[bcData], i++,
   sbc = bcData[[i]][[1]];
   time = bcData[[i]][[2]];
   sbc = sbc /. u -> Sqrt[Total[varSet^2] + 1];
   result = 
    Timing[TimeConstrained[
      myVerify[varSet, flowVec, initCons, unsafeCons, sbc, i], 
      verifyTime, 2]];
   totaltime = totaltime + time + result[[1]];
   totalsdptime = totalsdptime + time;
   If[result[[2]] == 1,
    sbcDeg = i;
    Print["Verified at deg ", i];
    (*Print[sbc];*)
    Break[]
    ];
   If[result[[2]] == 2,
    Print["Verify out of time at deg ", i];
    (*Print[sbc];*)
    sbcDeg = 0;
    Break[]
    ];
   ];
  htime = totaltime;
  Print["SDP time: ", Round[totalsdptime, 0.01] ];
  
  If[sbcDeg > 0,
   Print["Verified, SDP+verify time: ", Round[htime, 0.01]],
   Print["No valid soltion, SDP+verify time: ", Round[htime, 0.01]]
   ];
  If[ plot == 1,
   If[Length[varSet] == 2, 
    myPlot2D[benchmark, varSet, plotRange, flowVec, unsafeCons, 
     initCons, If[ebcDeg == 0, 1, ebc], If[hbcDeg == 0, 1, hbc], 
     If[sbcDeg == 0, 1, sbc]]];
   If[Length[varSet] == 3, 
    myPlot3D[benchmark, varSet, plotRange, flowVec, unsafeCons, 
     initCons, If[ebcDeg == 0, 1, ebc], If[hbcDeg == 0, 1, hbc], 
     If[sbcDeg == 0, 1, sbc]]];
   ]
  ]
myPlot2D [benchmark_, varSet_, varRange_, flowVec_, unsafeBoundary_, 
  initBoundary_, bc_, bch_, sbc_] := Module[
  {ranges, flowPlot, initCons, unsafeCons, initPlot, unsafePlot, 
   initSamples, flow, trajGen, trajTable, trajPlot, sampleNum, 
   solution, trajPlotAll, bcPlot, bchPlot, sbcPlot},
  Print["Plot2d..."];
  ranges = Sequence @@ MapThread[Prepend, {varRange, varSet}];
  flowPlot = 
   VectorPlot[flowVec, Evaluate[ranges], VectorScaling -> Automatic, 
    VectorSizes -> Automatic, VectorColorFunction -> None, 
    VectorStyle -> Gray];
  initCons = And @@ Thread[GreaterEqualThan[0][initBoundary]];
  unsafeCons = And @@ Thread[GreaterEqualThan[0][unsafeBoundary]];
  initPlot = 
   RegionPlot[initCons, Evaluate[ranges], 
    PlotStyle -> {Opacity[0.65], RGBColor[
      0.28026441037696703`, 0.715, 0.4292089322474965]}, 
    BoundaryStyle -> None , PlotPoints -> 100];
  unsafePlot = 
   RegionPlot[unsafeCons, Evaluate[ranges], 
    PlotStyle -> {Opacity[0.65], RGBColor[0.915, 0.3325, 0.2125]}, 
    BoundaryStyle -> None, PlotPoints -> 100];
  sampleNum = 5;
  initSamples = 
   RandomPoint[
    ImplicitRegion[
     initCons && 
      varRange[[1]][[1]] <= varSet[[1]] <= varRange[[1]][[2]] && 
      varRange[[2]][[1]] <= varSet[[2]] <= varRange[[2]][[2]], 
     varSet], sampleNum];
  Off[NDSolve::ndsz];
  solution[s0_] := 
   NDSolve[Join[
     Thread[#'[t] & /@ 
        varSet == (flowVec /. 
         Table[varSet[[i]] -> varSet[[i]][t], {i, 1, 
           Length[varSet]}])], Thread[#[0] & /@ varSet == s0]], 
    varSet, {t, 0, 20}, Method -> "StiffnessSwitching"]; 
  trajPlot = 
   ParametricPlot[
    Evaluate[(#[t] & /@ varSet) /. solution[#] & /@ initSamples], {t, 
     0, 20}, PlotStyle -> Directive[Black, Thickness[Large]]];
  bcPlot = 
   RegionPlot[bc <= 0, Evaluate[ranges], 
    PlotStyle -> {Opacity[0.15], RGBColor[1, 0.75, 0]}, 
    BoundaryStyle -> {
RGBColor[1, 0.75, 0]}, PlotPoints -> 100];
  bchPlot = 
   RegionPlot[bch <= 0, Evaluate[ranges], 
    PlotStyle -> {Opacity[0.15], RGBColor[
      0.363898, 0.618501, 0.782349]}, BoundaryStyle -> {
RGBColor[0.363898, 0.618501, 0.782349]}, PlotPoints -> 100];
  sbcPlot = 
   RegionPlot[sbc <= 0, Evaluate[ranges], 
    PlotStyle -> {Opacity[0.15], RGBColor[
      0.6900000000000001, 0.52, 0.6900000000000001]}, 
    BoundaryStyle -> {
RGBColor[0.6900000000000001, 0.52, 0.6900000000000001]}, 
    PlotPoints -> 100];
  Export[
   FileNameJoin[{".", "Results", "plots", 
     StringJoin[benchmark, ".pdf"]}], 
   Show[flowPlot, initPlot, unsafePlot, trajPlot, bcPlot, bchPlot, 
    sbcPlot, FrameLabel -> varSet]]
  ]

myPlot3D [benchmark_, varSet_, varRange_, flowVec_, unsafeBoundary_, 
   initBoundary_, bc_, bch_, sbc_] := Module[
   {ranges, flowPlot, initCons, unsafeCons, initPlot, unsafePlot, 
    initSamples, flow, trajGen, trajTable, trajPlot, sampleNum, 
    solution, trajPlotAll, bcPlot, bchPlot, sbcPlot},
   Print["Plot3d..."];
   ranges = Sequence @@ MapThread[Prepend, {varRange, varSet}];
   (*flowPlot =VectorPlot3D[flowVec,Evaluate[ranges],
   VectorScaling\[Rule]Automatic,VectorSizes\[Rule]Automatic,
   VectorColorFunction\[Rule]None,VectorStyle\[Rule]{Opacity[0.4],
   Gray}];*)
   (*plot init and unsafe regions*)
   initCons = And @@ Thread[GreaterEqualThan[0][initBoundary]];
   unsafeCons = And @@ Thread[GreaterEqualThan[0][unsafeBoundary]];
   initPlot = 
    RegionPlot3D[initCons, Evaluate[ranges], 
     PlotStyle -> {Opacity[0.65], RGBColor[
       0.28026441037696703`, 0.715, 0.4292089322474965]}, 
     BoundaryStyle -> None , PlotPoints -> 100, Mesh -> None];
   unsafePlot = 
    RegionPlot3D[unsafeCons, Evaluate[ranges], 
     PlotStyle -> {Opacity[0.65], RGBColor[0.915, 0.3325, 0.2125]}, 
     BoundaryStyle -> None, PlotPoints -> 100, Mesh -> None];
   (*plot sample trajectories*)
   sampleNum = 5;
   initSamples = 
    RandomPoint[
     ImplicitRegion[
      initCons && 
       varRange[[1]][[1]] <= varSet[[1]] <= varRange[[1]][[2]] && 
       varRange[[2]][[1]] <= varSet[[2]] <= varRange[[2]][[2]] && 
       varRange[[3]][[1]] <= varSet[[3]] <= varRange[[3]][[2]], 
      varSet], sampleNum];
   Off[NDSolve::ndsz];
   solution[s0_] := 
    NDSolve[Join[
      Thread[#'[t] & /@ 
         varSet == (flowVec /. 
          Table[varSet[[i]] -> varSet[[i]][t], {i, 1, 
            Length[varSet]}])], Thread[#[0] & /@ varSet == s0]], 
     varSet, {t, 0, 20}, Method -> "StiffnessSwitching"]; 
   trajPlot = 
    ParametricPlot3D[
     Evaluate[(#[t] & /@ varSet) /. solution[#] & /@ initSamples], {t,
       0, 20}, PlotStyle -> Directive[Black, Thickness[Large]]];
   bcPlot = 
    RegionPlot3D[bc <= 0, Evaluate[ranges], 
     PlotStyle -> {Opacity[0.3], RGBColor[1, 0.75, 0]}, 
     BoundaryStyle -> {
RGBColor[1, 0.75, 0]}, PlotPoints -> 100, Mesh -> None];
   bchPlot = 
    RegionPlot3D[bch <= 0, Evaluate[ranges], 
     PlotStyle -> {Opacity[0.3], RGBColor[
       0.363898, 0.618501, 0.782349]}, BoundaryStyle -> {
RGBColor[0.363898, 0.618501, 0.782349]}, PlotPoints -> 100, 
     Mesh -> None];
   sbcPlot = 
    RegionPlot3D[sbc <= 0, Evaluate[ranges], 
     PlotStyle -> {Opacity[0.3], RGBColor[
       0.6900000000000001, 0.52, 0.6900000000000001]}, 
     BoundaryStyle -> {
RGBColor[0.6900000000000001, 0.52, 0.6900000000000001]}, 
     PlotPoints -> 100, Mesh -> None];
   Export[
    FileNameJoin[{".", "Results", "plots", 
      StringJoin[benchmark, ".pdf"]}], 
    Show[initPlot, unsafePlot, trajPlot, bcPlot, bchPlot, sbcPlot, 
     FrameLabel -> varSet]]
   ];

benchmarks = {"lorenz-1", "lorenz-2"};
For[i = 1, i <= 2, i++,
 benchmarkname = benchmarks[[i]];
 Print["============", benchmarkname, "============"];
 verifytime = 600;
 plot = 1;
 plotRange = {{-20, 20}, {-20, 20}, {-20, 20}};
 verifyResults[benchmarkname, plot, plotRange, verifytime];
]