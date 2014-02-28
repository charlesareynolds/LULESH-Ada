package LULESH.Init is

   --x    // Constructor
   --x    Domain(Int_t numRanks, Index_t colLoc,
   --x           Index_t rowLoc, Index_t planeLoc,
   --x           Index_t nx, Int_t tp, Int_t nr, Int_t balance, Int_t cost);
   function Create
     (numRanks : in Int_t;
      colLoc   : in Index_Type;
      rowLoc   : in Index_Type;
      planeLoc : in Index_Type;
      nx       : in Index_Type;
      tp       : in Index_Type;
      nr       : in Int_t;
      balance  : in Int_t;
      cost     : in Int_t)
      return Domain_Record;

   --- // lulesh-init

   --x void InitMeshDecomp(Int_t numRanks, Int_t myRank,
   --x                     Int_t *col, Int_t *row, Int_t *plane, Int_t *side);
   procedure InitMeshDecomp
     (numRanks : in Int_t;
      myRank   : in Int_t;
      col      : out Int_t;
      row      : out Int_t;
      plane    : out Int_t;
      side     : out Int_t);

   --x    void BuildMesh(Int_t nx, Int_t edgeNodes, Int_t edgeElems);
   procedure BuildMesh
     (this      : in out Domain_Record;
      nx        : in Index_Type;
      edgeNodes : in Node_Index_Type;
      edgeElems : in Element_Index_Type);

   --x    void SetupThreadSupportStructures();
   procedure SetupThreadSupportStructures;

   --x    void CreateRegionIndexSets(Int_t nreg, Int_t balance);
   procedure CreateRegionIndexSets
     (nreg    : in Int_t;
      balance : in Int_t);

   --x    void SetupCommBuffers(Int_t edgeNodes);
   procedure SetupCommBuffers
     (edgeNodes : in Node_Index_Type);

   --x    void SetupSymmetryPlanes(Int_t edgeNodes);
   procedure SetupSymmetryPlanes
     (edgeNodes : in Node_Index_Type);

   --x    void SetupElementConnectivities(Int_t edgeElems);
   procedure SetupElementConnectivities
     (edgeElems : in Element_Index_Type);

   --x    void SetupBoundaryConditions(Int_t edgeElems);
   procedure SetupBoundaryConditions
     (this      : in out Domain_Record;
      edgeElems : in Element_Index_Type);

end LULESH.Init;
