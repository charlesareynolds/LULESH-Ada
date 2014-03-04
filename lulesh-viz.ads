package LULESH.Viz is
   --- // lulesh-viz

   --x void DumpToVisit(Domain& domain, int numFiles, int myRank, int numRanks);
   procedure DumpToVisit
     (domainn  : in out Domain_Record;
      numFiles : in Int_t;
      myRank   : in Rank_Type;
      numRanks : in Rank_Type);

end LULESH.Viz;
