(* Codes for plots... *)

Needs["Misc`npPlotting`"];

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

(* Read in John's RR counts 

Three functions :
  readjp : reads in one file and dumps out a normalized array
  readalljp : reads in all the files
  processjp : data in a similar format as previously.
*)

readjp[fn_] := Module[{dat, norm}, 
	dat = Import[fn,"Table"];
	norm = dat[[1,1]];
	dat[[3;;,3]]/norm^2
]

readalljp[dir_, npts_] := Module[{fnlist},
	fnlist = FileNames[dir<>IntegerString[npts]<>"*.dat"];
	readjp /@ fnlist
];