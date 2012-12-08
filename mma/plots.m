(* Codes for plots... *)

(* Some shapes etc *)
pSquare = Graphics[{Rectangle[]}];
pDisk = Graphics[{Disk[]}];

Clear[npLogTicks];
npLogTicks[lnstart_, lnend_, label_:True] :=
    Module[ {bigticks, smallticks},
        smallticks = 
         Table[{10.0^i*j, "", {0.01, 0}}, {i, lnstart, lnend}, {j, 2., 
            9}] // Flatten[#, {{1, 2}, {3}}] &;
        bigticks = 
         Table[{10.0^i, If[ label,
                            Superscript[10, i],
                            ""
                        ], {0.02, 0}}, {i, 
           lnstart, lnend}];
        Join[bigticks, smallticks]
    ];


(* Code to read in the binary files we generate *)

readSims[fn_] :=
    Module[ {hdr, ff, out},
        ff = OpenRead[fn, BinaryFormat->True];
        hdr = BinaryRead[ff, {"Integer32", "Integer32", "Real64", "Real64"}];
        out = BinaryReadList[ff, "Real64", hdr[[1]]];
        Close[ff];
        {hdr, out}
    ];

stats[ll_List] := {Mean[ll], StandardDeviation[ll], StandardDeviation[ll]/Mean[ll]};

genfilenames[hdr_, ibin_, npts_] := hdr <> "_" <> IntegerString[ibin, 10, 2]  \
	<> "_" <> IntegerString[npts, 10, 8] <> ".dat";

processbin[hdr_, ibin_, nptlist_] := With[{tmp=readSims[genfilenames[hdr, ibin, #]]}, {tmp[[1]],stats[tmp[[2]]]}] & /@ nptlist;