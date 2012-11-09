(* Collection of code for generating low discrepancy sequences *)

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
    Plus[x0, #] & /@  generateLDSGrid[num, dim] // FractionalPart
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