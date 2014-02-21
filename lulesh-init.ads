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
      tp       : in Int_t;
      nr       : in Int_t;
      balance  : in Int_t;
      cost     : in Int_t)
      return Access_Domain;

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
     (nx        : in Int_t;
      edgeNodes : in Int_t;
      edgeElems : in Int_t);

   --x    void SetupThreadSupportStructures();
   procedure SetupThreadSupportStructures;

   --x    void CreateRegionIndexSets(Int_t nreg, Int_t balance);
   procedure CreateRegionIndexSets
     (nreg    : in Int_t;
      balance : in Int_t);

   --x    void SetupCommBuffers(Int_t edgeNodes);
   procedure SetupCommBuffers
     (edgeNodes : in Int_t);

   --x    void SetupSymmetryPlanes(Int_t edgeNodes);
   procedure SetupSymmetryPlanes
     (edgeNodes : in Int_t);

   --x    void SetupElementConnectivities(Int_t edgeElems);
   procedure SetupElementConnectivities
     (edgeElems : in Int_t);

   --x    void SetupBoundaryConditions(Int_t edgeElems);
   procedure SetupBoundaryConditions
     (edgeElems : in Int_t);

end LULESH.Init;
