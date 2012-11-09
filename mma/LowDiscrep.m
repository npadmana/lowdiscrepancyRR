(* Collection of code for generating low discrepancy sequences *)

(* Reshape arrays *)
Clear[reshape];
reshape[ll_, {d__?((IntegerQ[#] && Positive[#]) &)}] := 
  With[{e = Flatten[ll]}, 
   Fold[Partition, e, 
     Take[{d}, {-1, 2, -1}]] /; (Length[e] === Times[d])];

(* Wrapper function for low discrepancy sequences in dimension d *)
Clear[generateLDSGrid];
generateLDSGrid[num_, dim_:2] := BlockRandom[
      SeedRandom[
       Method -> {"MKL", 
         Method -> {"Niederreiter", "Dimension" -> dim}}];
      RandomReal[1, {num, dim}]];
      

(* Generate a shifted low discrepancy sequence in dimension dim *)
Clear[generateLDSshiftedGrid];
generateLDSshiftedGrid[num_, dim_:2] := With[{x0 = RandomReal[{0, 1}, dim]}, 
    Plus[x0, #] & /@  generateLDSGrid[num, dim] // FractionalPart // reshape[#, {num, 2, dim/2}] &
   ];
   
(* Generate a purely random grid *)   
Clear[generateRandomGrid];
generateRandomGrid[snum_, dim_:1] := With[{xp = RandomReal[{0, 1}, {snum, dim}]}, 
   Distribute[{xp, xp}, List]
   ];   
   

(* Run func nsamp times, storing mean and standard deviation 

getmeanerror has attribute HoldFirst so the function is re-evaluated every time.
*)
Clear[getmeanerror];
getmeanerror[func_, nsamp_] := Module[{ll},
   ll = Reap[Do[Sow[func], {nsamp}]][[2, 1]];
   {Mean[ll], StandardDeviation[ll]}
   ];
SetAttributes[getmeanerror, HoldFirst];

(* Basic Monte Carlo Integrator

Applies a function to a list of points, and sums them, and divides by the numberof points
 *)
Clear[mcIntegrate];
mcIntegrate[func_, ll_] := Total[func @@@ ll]/Dimensions[ll][[1]];
SetAttributes[mcIntegrate, HoldFirst];

