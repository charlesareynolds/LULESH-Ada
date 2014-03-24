package LULESH.Init is

   --x    // Constructor
   --x    Domain(Int_t numRanks, Index_t colLoc,
   --x           Index_t rowLoc, Index_t planeLoc,
   --x           Index_t nx, Int_t tp, Int_t nr, Int_t balance, Int_t cost);
   function Create
     (numRanks    : in Rank_Count_Range;
      colLoc      : in Domain_Index;
      rowLoc      : in Domain_Index;
      planeLoc    : in Domain_Index;
      side_length : in Element_Index;
      tp          : in Domain_Index;
      nr          : in Region_Index;
      balance     : in Balance_Type;
      cost        : in Cost_Type)
      return Domain_Record;

   --- // lulesh-init

   --x void InitMeshDecomp(Int_t numRanks, Int_t myRank,
   --x                     Int_t *col, Int_t *row, Int_t *plane, Int_t *side);
   procedure InitMeshDecomp
     (numRanks         : in Rank_Count_Range;
      myRank           : in Rank_Type;
      domain_column    : out Domain_Index;
      domain_row       : out Domain_Index;
      domain_plane     : out Domain_Index;
      domains_per_side : out Domain_Index);

   --x    void BuildMesh(Int_t nx, Int_t edgeNodes, Int_t edgeElems);
   procedure BuildMesh
     (this        : in out Domain_Record;
      side_length : in Element_Index;
      edgeNodes   : in Node_Index;
      edgeElems   : in Element_Index);

   --x    void SetupThreadSupportStructures();
   procedure SetupThreadSupportStructures
     (this : in out Domain_Record);

   --x    void CreateRegionIndexSets(Int_t nreg, Int_t balance);
   procedure CreateRegionIndexSets
     (this    : in out Domain_Record;
      nreg    : in Region_Index;
      balance : in Balance_Type);

   --x    void SetupCommBuffers(Int_t edgeNodes);
   procedure SetupCommBuffers
     (This      : in out Domain_Record;
      EdgeNodes : in Node_Index);

   --x    void SetupSymmetryPlanes(Int_t edgeNodes);
   procedure SetupSymmetryPlanes
     (this      : in out Domain_Record;
      edgeNodes : in Node_Index);

   --x    void SetupElementConnectivities(Int_t edgeElems);
   procedure SetupElementConnectivities
     (this      : in out Domain_Record;
      edgeElems : in Element_Index);

   --x    void SetupBoundaryConditions(Int_t edgeElems);
   procedure SetupBoundaryConditions
     (this      : in out Domain_Record;
      edgeElems : in Element_Index);

end LULESH.Init;
