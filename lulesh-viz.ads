package LULESH.Viz is
   --- // lulesh-viz

   --x void DumpToVisit(Domain& domain, int numFiles, int myRank, int numRanks);
   procedure DumpToVisit
     (domain  : in out Domain_Record;
      numFiles : in Integer;
      myRank   : in Rank_Type;
      numRanks : in Rank_Type);

end LULESH.Viz;
