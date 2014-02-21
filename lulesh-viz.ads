package LULESH.Viz is
   --- // lulesh-viz

   --x void DumpToVisit(Domain& domain, int numFiles, int myRank, int numRanks);
   procedure DumpToVisit
     (domainn  : not null access Domain_Record;
      numFiles : in Int_t;
      myRank   : in Int_t;
      numRanks : in Int_t);

end LULESH.Viz;
