package LULESH.Viz is
   --- // lulesh-viz

   --x void DumpToVisit(Domain& domain, int numFiles, int myRank, int numRanks);
   procedure DumpToVisit
     (domain  : in out Domain_Record;
      numFiles : in Int_t;
      myRank   : in MPI.Rank_Type;
      numRanks : in MPI.Rank_Type);

end LULESH.Viz;
