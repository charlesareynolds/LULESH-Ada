package LULESH.Viz is
   --- // lulesh-viz

   --x void DumpToVisit(Domain& domain, int numFiles, int myRank, int numRanks);
   procedure DumpToVisit
     (domainn  : not null access Domain_Record;
      numFiles : in IC.int;
      myRank   : in IC.int;
      numRanks : in IC.int);

end LULESH.Viz;
