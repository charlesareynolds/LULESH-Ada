--x #include <math.h>
--x #if USE_MPI
--x # include <mpi.h>
--x #endif
--x #if _OPENMP
--x #include <omp.h>
--x #endif
--x #include <stdio.h>
--x #include <stdlib.h>
--x #include <string.h>
--x #include <limits.h>
--x #include <cstdlib>
--x #include "lulesh.h"

with Ada.Exceptions;
with Ada.Numerics.Float_Random;
with LULESH.Par;

package body LULESH.Init is

   package AEX renames Ada.Exceptions;

   package Random_Selection is

      --- Returns an integer between Min and Max:
      function Choose
        (Minimum : in Integer;
         Maximum : in Integer)
         return Integer
        with
          Post => (Choose'Result in Minimum .. Maximum);

      --- Returns an integer between 0 and Modulo - 1:
      function Choose_Rem
        (Modulo : in Integer)
         return Integer is
        (Choose (0, Modulo - 1));

      --- Intialzes the generator. Initializing with the same seed will produce
      --- the same sequence of outputs.
      procedure Initialize (Seed : in Integer);

   end Random_Selection;

   package body Random_Selection is
      package ANFR renames Ada.Numerics.Float_Random;
      Generator : ANFR.Generator;

      procedure Initialize (Seed : in Integer) is
      begin
         ANFR.Reset (Generator, Seed);
      end Initialize;

      function Choose
        (Minimum : in Integer;
         Maximum : in Integer)
         return Integer is
      begin
         return Integer (Float (ANFR.Random (Generator)) * Float (Maximum - Minimum)) + Minimum;
      end Choose;
   end Random_Selection;

   function Choose_Rem
     (Modulo : in Integer)
      return Element_Index is
     (Element_Index (Random_Selection.Choose_Rem (Modulo)));

   function Choose_Rem
     (Modulo : in Integer)
      return Int_T is
     (Int_T (Random_Selection.Choose_Rem (Modulo)));

   function Choose_Rem
     (Modulo : in Cost_Type)
      return Cost_Type is
     (Cost_Type (Random_Selection.Choose_Rem (Integer (Modulo))));

   procedure Abort_Or_Raise
     (X       : in AEX.Exception_Id;
      Message : in String) is
   begin
      --z #if USE_MPI
      --z       MPI_Abort(MPI_COMM_WORLD, -1) ;
      --z #else
      --x       exit(-1);
      --z #endif
      if USE_MPI then
         MPI.Abortt (MPI.COMM_WORLD, MPI.Errorcode_Type (-1));
      else
         AEX.Raise_Exception (X, Message);
      end if;
   end Abort_Or_Raise;

   ---    public:

   ---    //
   ---    // ALLOCATION
   ---    //

   --x    void AllocateNodePersistent(Int_t numNode) // Node-centered
   --x    {
   --x       m_x.resize(numNode);  // coordinates
   --x       m_y.resize(numNode);
   --x       m_z.resize(numNode);

   --x       m_xd.resize(numNode); // velocities
   --x       m_yd.resize(numNode);
   --x       m_zd.resize(numNode);

   --x       m_xdd.resize(numNode); // accelerations
   --x       m_ydd.resize(numNode);
   --x       m_zdd.resize(numNode);

   --x       m_fx.resize(numNode);  // forces
   --x       m_fy.resize(numNode);
   --x       m_fz.resize(numNode);

   --x       m_nodalMass.resize(numNode);  // mass
   --x    }

   --x    void AllocateElemPersistent(Int_t numElem) // Elem-centered
   --x    {
   --x       m_nodelist.resize(8*numElem);

   ---       // elem connectivities through face
   --x       m_lxim.resize(numElem);
   --x       m_lxip.resize(numElem);
   --x       m_letam.resize(numElem);
   --x       m_letap.resize(numElem);
   --x       m_lzetam.resize(numElem);
   --x       m_lzetap.resize(numElem);

   --x       m_elemBC.resize(numElem);

   --x       m_e.resize(numElem);
   --x       m_p.resize(numElem);

   --x       m_q.resize(numElem);
   --x       m_ql.resize(numElem);
   --x       m_qq.resize(numElem);

   --x       m_v.resize(numElem);

   --x       m_volo.resize(numElem);
   --x       m_delv.resize(numElem);
   --x       m_vdov.resize(numElem);

   --x       m_arealg.resize(numElem);

   --x       m_ss.resize(numElem);

   --x       m_elemMass.resize(numElem);
   --x    }

   --- AllocateGradients and DeallocateGradients are not needed if gradient is in
   --- Element_Record:
   --x    void AllocateGradients(Int_t numElem, Int_t allElem)
   --x    {
   --x       // Position gradients
   --x       m_delx_xi.resize(numElem) ;
   --x       m_delx_eta.resize(numElem) ;
   --x       m_delx_zeta.resize(numElem) ;
   --x       // Velocity gradients
   --x       m_delv_xi.resize(allElem) ;
   --x       m_delv_eta.resize(allElem);
   --x       m_delv_zeta.resize(allElem) ;
   --x    }
   --x    void DeallocateGradients()
   --x    {
   --x       m_delx_zeta.clear() ;
   --x       m_delx_eta.clear() ;
   --x       m_delx_xi.clear() ;
   --x       m_delv_zeta.clear() ;
   --x       m_delv_eta.clear() ;
   --x       m_delv_xi.clear() ;
   --x    }
   --- AllocateStrains and DeallocateStrains are not needed if strain is in
   --- Element_Record:
   --x    void AllocateStrains(Int_t numElem)
   --x    {
   --x       m_dxx.resize(numElem) ;
   --x       m_dyy.resize(numElem) ;
   --x       m_dzz.resize(numElem) ;
   --x    }
   --x    void DeallocateStrains()
   --x    {
   --x       m_dzz.clear() ;
   --x       m_dyy.clear() ;
   --x       m_dxx.clear() ;
   --x    }
   --- /////////////////////////////////////////////////////////////////////
   --x Domain::Domain(Int_t numRanks, Index_t colLoc,
   --x                Index_t rowLoc, Index_t planeLoc,
   --x                Index_t nx, int tp, int nr, int balance, Int_t cost)
   --x    :
   -------------
   --- EXPORTED:
   -------------
   function Create
     (NumRanks    : in Rank_Count_Range;
      ColLoc      : in Domain_Index;
      RowLoc      : in Domain_Index;
      PlaneLoc    : in Domain_Index;
      Side_Length : in Element_Index;
      Tp          : in Domain_Index;
      Nr          : in Region_Index;
      Balance     : in Balance_Type;
      Cost        : in Cost_Type)
      return Domain_Record
   is
      This : Domain_Record;
      --x    Index_t edgeElems = nx ;
      --x    Index_t edgeNodes = edgeElems+1 ;
      EdgeElems : constant Element_Index := Side_Length;
      EdgeNodes : constant Node_Index := Node_Index (Side_Length) + 1;
   begin

      --x    m_e_cut(Real_t(1.0e-7)),
      --x    m_p_cut(Real_t(1.0e-7)),
      --x    m_q_cut(Real_t(1.0e-7)),
      --x    m_v_cut(Real_t(1.0e-10)),
      --x    m_u_cut(Real_t(1.0e-7)),
      --x    m_hgcoef(Real_t(3.0)),
      --x    m_ss4o3(Real_t(4.0)/Real_t(3.0)),
      --x    m_qstop(Real_t(1.0e+12)),
      --x    m_monoq_max_slope(Real_t(1.0)),
      --x    m_monoq_limiter_mult(Real_t(2.0)),
      --x    m_qlc_monoq(Real_t(0.5)),
      --x    m_qqc_monoq(Real_t(2.0)/Real_t(3.0)),
      --x    m_qqc(Real_t(2.0)),
      --x    m_eosvmax(Real_t(1.0e+9)),
      --x    m_eosvmin(Real_t(1.0e-9)),
      --x    m_pmin(Real_t(0.)),
      --x    m_emin(Real_t(-1.0e+15)),
      --x    m_dvovmax(Real_t(0.1)),
      --x    m_refdens(Real_t(1.0))
      --x {
      This.Parameters :=
        (Energy_Tolerance           => 1.0e-7 * J,
         Pressure_Tolerance         => 1.0e-7 * Pa,
         Artificial_Viscosity_Tolerance => 1.0e-7,
         Volume_Relative_Tolerance  => 1.0e-10,
         Velocity_Tolerance         => 1.0e-7 * Mps,
         Hgcoef                     => 3.0,
         Four_Thirds                => 4.0 / 3.0,
         Qstop                      => 1.0e+12,
         Monoq_Max_Slope            => 1.0,
         Monoq_Limiter_Mult         => 2.0,
         Qlc_Monoq                  => 0.5,
         Qqc_Monoq                  => 2.0 / 3.0,
         Qqc                        => 2.0,
         Eosvmax                    => 1.0e+9 * M3,
         Eosvmin                    => 1.0e-9 * M3,
         Pressure_Floor             => 0.0 * Pa,
         Energy_Floor               => -1.0e+15 * J,
         Dvovmax                    => Compression (0.1),
         Reference_Density          => 1.0 * Kgpm3,
         --x    this->cost() = cost;
         Imbalance_Cost             => Cost,
         --x    m_sizeX = edgeElems ;
         --x    m_sizeY = edgeElems ;
         --x    m_sizeZ = edgeElems ;
         Size                       => (others => EdgeElems));
      --x    m_tp       = tp ;
      This.Variables.Tp := Tp;
      --x    m_numRanks = numRanks ;
      This.Variables.NumRanks := NumRanks;

      ---    ///////////////////////////////
      ---    //   Initialize Sedov Mesh
      ---    ///////////////////////////////

      ---    // construct a uniform box for this processor

      --x    m_colLoc   =   colLoc ;
      --x    m_rowLoc   =   rowLoc ;
      --x    m_planeLoc = planeLoc ;
      This.Variables.ColLoc   := ColLoc;
      This.Variables.RowLoc   := RowLoc;
      This.Variables.PlaneLoc := PlaneLoc;
      --x    m_numElem = edgeElems*edgeElems*edgeElems ;
      --x    m_numNode = edgeNodes*edgeNodes*edgeNodes ;
      This.NumElem := EdgeElems ** 3;
      This.NumNode := EdgeNodes ** 3;

      --x    m_regNumList = new Index_t[numElem()] ;  // material indexset

      --x    // Elem-centered
      --x    AllocateElemPersistent(numElem()) ;
      This.Elements := new Element_Array (0 .. This.NumElem - 1);

      --x    // Node-centered
      --x    AllocateNodePersistent(numNode()) ;
      This.Nodes := new Node_Array (0 .. This.NumNode - 1);

      --x    SetupCommBuffers(edgeNodes);
      SetupCommBuffers (This, EdgeNodes);

      --x    // Basic Field Initialization
      --x    for (Index_t i=0; i<numElem(); ++i) {
      --x       e(i) =  Real_t(0.0) ;
      --x       p(i) =  Real_t(0.0) ;
      --x       q(i) =  Real_t(0.0) ;
      --x       ss(i) = Real_t(0.0) ;
      --x    }
      ---    // Note - v initializes to 1.0, not 0.0!
      --x    for (Index_t i=0; i<numElem(); ++i) {
      --x       v(i) = Real_t(1.0) ;
      --x    }
      for Index in 0 .. This.NumElem - 1 loop
         This.Elements (Index).Energy               := 0.0 * J;
         This.Elements (Index).Pressure             := 0.0 * Pa;
         This.Elements (Index).Artificial_Viscosity := 0.0;
         This.Elements (Index).Sound_Speed          := 0.0 * Mps;
         This.Elements (Index).Volume               := 1.0;
      end loop;
      --x    for (Index_t i=0; i<numNode(); ++i) {
      --x       xd(i) = Real_t(0.0) ;
      --x       yd(i) = Real_t(0.0) ;
      --x       zd(i) = Real_t(0.0) ;
      --x    }
      --x    for (Index_t i=0; i<numNode(); ++i) {
      --x       xdd(i) = Real_t(0.0) ;
      --x       ydd(i) = Real_t(0.0) ;
      --x       zdd(i) = Real_t(0.0) ;
      --x    }
      --x    for (Index_t i=0; i<numNode(); ++i) {
      --x       nodalMass(i) = Real_t(0.0) ;
      --x    }
      for Index in 0 .. This.NumNode - 1 loop
         This.Nodes (Index).Velocity     := (others => 0.0 * Mps);
         This.Nodes (Index).Acceleration := (others => 0.0 * Mps2);
         This.Nodes (Index).Mass         := 0.0 * Kg;
      end loop;

      --x    BuildMesh(nx, edgeNodes, edgeElems);
      BuildMesh (This, Side_Length, EdgeNodes, EdgeElems);

      --x #if _OPENMP
      --x    SetupThreadSupportStructures();
      --x #else
      --x    // These arrays are not used if we're not threaded
      --x    m_nodeElemStart = NULL;
      --x    m_nodeElemCornerList = NULL;
      --x #endif
      if COMPILER_SUPPORTS_OPENMP then
         SetupThreadSupportStructures (This);
      else
         This.Variables.NodeElemStart := null;
         This.Variables.NodeElemCornerList := null;
      end if;

      ---    // Setup region index sets. For now, these are constant sized
      ---    // throughout the run, but could be changed every cycle to
      ---    // simulate effects of ALE on the lagrange solver
      --x    CreateRegionIndexSets(nr, balance);
      ---    // Setup symmetry nodesets
      --x    SetupSymmetryPlanes(edgeNodes);
      ---    // Setup element connectivities
      --x    SetupElementConnectivities(edgeElems);
      ---    // Setup symmetry planes and free surface boundary arrays
      --x    SetupBoundaryConditions(edgeElems);
      CreateRegionIndexSets (This, Nr, Balance);
      SetupSymmetryPlanes (This, EdgeNodes);
      SetupElementConnectivities (This, EdgeElems);
      SetupBoundaryConditions (This, EdgeElems);

      ---    // Setup defaults

      ---    // These can be changed (requires recompile) if you want to run
      ---    // with a fixed timestep, or to a different end time, but it's
      ---    // probably easier/better to just run a fixed number of timesteps
      ---    // using the -i flag in 2.x

      --x    dtfixed() = Real_t(-1.0e-6) ; // Negative means use courant condition
      --x    stoptime()  = Real_t(1.0e-2); // *Real_t(edgeElems*tp/45.0) ;
      This.Variables.Dtfixed  := Time_Span_Last;
      This.Variables.Use_Courant_Condition  := True;
      This.Variables.Stoptime := Time_First + 1.0e-2 * S;

      ---    // Initial conditions
      --x    deltatimemultlb() = Real_t(1.1) ;
      --x    deltatimemultub() = Real_t(1.2) ;
      --x    dtcourant() = Real_t(1.0e+20) ;
      --x    dthydro()   = Real_t(1.0e+20) ;
      --x    dtmax()     = Real_t(1.0e-2) ;
      --x    time()    = Real_t(0.) ;
      --x    cycle()   = Int_t(0) ;
      This.Variables.Delta_Time_Multiplier_Lower_Bound := 1.1;
      This.Variables.Delta_Time_Multiplier_Upper_Bound := 1.2;
      This.Variables.Dtcourant       := Time_Span_Last;
      This.Variables.Dthydro         := Time_Span_Last;
      This.Variables.Dtmax           := 1.0e-2 * S;
      This.Variables.Current_Time    := Time_First;
      This.Variables.Cycle           := 0;

      ---    // initialize field data
      --x    for (Index_t i=0; i<numElem(); ++i) {
      --x       Real_t x_local[8], y_local[8], z_local[8] ;
      --x       Index_t *elemToNode = nodelist(i) ;
      --x       for( Index_t lnode=0 ; lnode<8 ; ++lnode )
      --x       {
      --x         Index_t gnode = elemToNode[lnode];
      --x         x_local[lnode] = x(gnode);
      --x         y_local[lnode] = y(gnode);
      --x         z_local[lnode] = z(gnode);
      --x       }
      ---       // volume calculations
      --x       Real_t volume = CalcElemVolume(x_local, y_local, z_local );
      --x       volo(i) = volume ;
      --x       elemMass(i) = volume ;
      --x       for (Index_t j=0; j<8; ++j) {
      --x          Index_t idx = elemToNode[j] ;
      --x          nodalMass(idx) += volume / Real_t(8.0) ;
      --x       }
      --x    }
      for Element in 0 .. This.NumElem - 1 loop
         declare
            Local_Coords : NodesPerElement_Coordinate_Array;
            ElemToNode   : constant NodesPerElement_Element_Index_Array :=
              This.Elements (Element).Node_Indexes;
         begin
            for Node in NodesPerElement_Range loop
               Local_Coords (Node) :=
                 This.Nodes (ElemToNode (Node)).Coordinate;
            end loop;
            declare
               Element_Volume : constant Volume :=
                 LULESH.Par.CalcElemVolume (Local_Coords);
               Mass_Share     : constant Mass := Mass
                 (Element_Volume / Dimensionless (NODES_PER_ELEMENT));
            begin
               This.Elements (Element).Volume_Reference := Element_Volume;
               This.Elements (Element).Mass := Mass_Type (Element_Volume);
               for Node in NodesPerElement_Range loop
                  declare
                     Node_Mass : Mass renames
                       This.Nodes (ElemToNode (Node)).Mass;
                  begin
                     Node_Mass := Node_Mass + Mass_Share;
                  end;
               end loop;
            end;
         end;
      end loop;

      ---    // deposit initial energy
      ---    // An energy of 3.948746e+7 is correct for a problem with
      ---    // 45 zones along a side - we need to scale it

      --x    const Real_t ebase = Real_t(3.948746e+7);
      --x    Real_t scale = (nx*m_tp)/Real_t(45.0);
      --x    Real_t einit = ebase*scale*scale*scale;
      --x    if (m_rowLoc + m_colLoc + m_planeLoc == 0) {
      ---       // Dump into the first zone (which we know is in the corner)
      ---       // of the domain that sits at the origin
      --x       e(0) = einit;
      --x    }
      declare
         Ebase : constant Energy := 3.948746e+7 * J;
         Scale : constant Dimensionless   :=
           Dimensionless (Side_Length) * Dimensionless (This.Variables.Tp) / 45.0;
         Einit : constant Energy := Ebase * Scale ** 3;
      begin
         if (This.Variables.RowLoc
             + This.Variables.ColLoc
             + This.Variables.PlaneLoc = 0) then
            This.Elements (0).Energy := Einit;
         end if;
         ---    //set initial deltatime base on analytic CFL calculation
         --x    deltatime() = (Real_t(.5)*cbrt(volo(0)))/sqrt(Real_t(2.0)*einit);
         This.Variables.Deltatime := Time (
                                           ((0.5 *  Cbrt (Dimensionless (This.Elements (0).Volume_Reference))) /
                                             Sqrt (2.0 * Dimensionless (Einit))));
         --             ((0.5 * Cbrt (Real10 (This.Elements (0).Volume_Reference))) /
         --                Sqrt (2.0 * Real10 (Einit))) * S;
      end;
      --x } // End constructor
      return This;
   end Create;


   --- ////////////////////////////////////////////////////////////////////////////////
   --x void
   --x Domain::BuildMesh(Int_t nx, Int_t edgeNodes, Int_t edgeElems)
   --x {
   procedure BuildMesh
     (This        : in out Domain_Record;
      Side_Length : in Element_Index;
      EdgeNodes   : in Node_Index;
      EdgeElems   : in Element_Index)
   is
      subtype Edge_Elements_Range is Element_Index range 0 .. EdgeElems - 1;
      subtype Edge_Nodes_Range is Node_Index range 0 .. EdgeNodes - 1;
      --x   Index_t meshEdgeElems = m_tp*nx ;
      ---   // initialize nodal coordinates
      --x   Index_t nidx = 0 ;
      MeshEdgeElems : constant Element_Index :=
        Element_Index (This.Variables.Tp) * Side_Length ;
      Node          : Node_Index;
      T             : Coordinate_Vector;
      Zidx          : Element_Index;

      function Calc_T_Part
        (Loc   : in Domain_Index;
         Nodee : in Node_Index)
        --- Coordinate vector components aren't dimensioned, so:
         return Dimensionless is
        (1.125 * Dimensionless
           ((Element_Index (Loc) * Side_Length)+
                Element_Index (Nodee) / MeshEdgeElems))
          with Inline;
   begin
      --x   Real_t tz = Real_t(1.125)*Real_t(m_planeLoc*nx)/Real_t(meshEdgeElems) ;
      --x   for (Index_t plane=0; plane<edgeNodes; ++plane) {
      --x     Real_t ty = Real_t(1.125)*Real_t(m_rowLoc*nx)/Real_t(meshEdgeElems) ;
      --x     for (Index_t row=0; row<edgeNodes; ++row) {
      --x       Real_t tx = Real_t(1.125)*Real_t(m_colLoc*nx)/Real_t(meshEdgeElems) ;
      --x       for (Index_t col=0; col<edgeNodes; ++col) {
      --x 	x(nidx) = tx ;
      --x 	y(nidx) = ty ;
      --x 	z(nidx) = tz ;
      --x 	++nidx ;
      --x 	// tx += ds ; // may accumulate roundoff...
      --x 	tx = Real_t(1.125)*Real_t(m_colLoc*nx+col+1)/Real_t(meshEdgeElems) ;
      --x       }
      --x       // ty += ds ;  // may accumulate roundoff...
      --x       ty = Real_t(1.125)*Real_t(m_rowLoc*nx+row+1)/Real_t(meshEdgeElems) ;
      --x     }
      --x     // tz += ds ;  // may accumulate roundoff...
      --x     tz = Real_t(1.125)*Real_t(m_planeLoc*nx+plane+1)/Real_t(meshEdgeElems) ;
      --x   }
      Node := 0;
      for This_Plane in Edge_Nodes_Range loop
         T (Z) := Calc_T_Part (This.Variables.PlaneLoc, This_Plane);
         for This_Row in Edge_Nodes_Range loop
            T (Y) := Calc_T_Part (This.Variables.RowLoc, This_Row);
            for This_Col in Edge_Nodes_Range loop
               T (X) := Calc_T_Part (This.Variables.ColLoc, This_Col);
               This.Nodes (Node).Coordinate := T;
               Node := Node + 1;
            end loop;
         end loop;
      end loop;

      ---   // embed hexehedral elements in nodal point lattice
      --x   Index_t zidx = 0 ;
      --x   nidx = 0 ;
      --x   for (Index_t plane=0; plane<edgeElems; ++plane) {
      --x     for (Index_t row=0; row<edgeElems; ++row) {
      --x       for (Index_t col=0; col<edgeElems; ++col) {
      --x 	Index_t *localNode = nodelist(zidx) ;
      --x 	localNode[0] = nidx                                       ;
      --x 	localNode[1] = nidx                                   + 1 ;
      --x 	localNode[2] = nidx                       + edgeNodes + 1 ;
      --x 	localNode[3] = nidx                       + edgeNodes     ;
      --x 	localNode[4] = nidx + edgeNodes*edgeNodes                 ;
      --x 	localNode[5] = nidx + edgeNodes*edgeNodes             + 1 ;
      --x 	localNode[6] = nidx + edgeNodes*edgeNodes + edgeNodes + 1 ;
      --x 	localNode[7] = nidx + edgeNodes*edgeNodes + edgeNodes     ;
      --x 	++zidx ;
      --x 	++nidx ;
      --x       }
      --x       ++nidx ;
      --x     }
      --x     nidx += edgeNodes ;
      --x   }
      Zidx := 0 ;
      Node := 0;
      for This_Plane in Edge_Elements_Range loop
         for This_Row in Edge_Elements_Range loop
            for This_Col in Edge_Elements_Range loop
               This.Elements (Zidx).Node_Indexes :=
                 (0 => Node                                       ,
                  1 => Node                                   + 1 ,
                  2 => Node                       + EdgeNodes + 1 ,
                  3 => Node                       + EdgeNodes     ,
                  4 => Node + EdgeNodes ** 2                 ,
                  5 => Node + EdgeNodes ** 2             + 1 ,
                  6 => Node + EdgeNodes ** 2 + EdgeNodes + 1 ,
                  7 => Node + EdgeNodes ** 2 + EdgeNodes     );
               Zidx := Zidx + 1;
               Node := Node + 1;
            end loop;
            Node := Node + 1;
         end loop;
         Node := Node + EdgeNodes ;
      end loop;
      --x }
   end BuildMesh;

   --x ////////////////////////////////////////////////////////////////////////////////
   --x void
   --x Domain::SetupThreadSupportStructures()
   --x {
   procedure SetupThreadSupportStructures
     (This : in out Domain_Record)
   is
      --x #if _OPENMP
      --x    Index_t numthreads = omp_get_max_threads();
      --x #else
      --x    Index_t numthreads = 1;
      --x #endif
      Numthreads         : constant Index_Type :=
        (if COMPILER_SUPPORTS_OPENMP then OMP.Get_Max_Threads else 1);
      --- Number of elements each node is in. Interior nodes are in 8 elements,
      --- faces 4, edges 2, and corners 1:
      NodeElemCounts     : Node_Element_Index_Array_Access;
      NodeElemStart      : Node_Element_Index_Array_Access renames
        This.Variables.NodeElemStart;
      NodeElemCornerList : Element_Element_Index_Array_Access renames
        This.Variables.NodeElemCornerList;

      procedure Count_Elements_For_Each_Node is
      begin
         NodeElemCounts := new Node_Element_Index_Array (0 .. This.NumNode - 1);
         NodeElemCounts.all := (others => 0);
         --x     for (Index_t i=0; i<numElem(); ++i) {
         --x       Index_t *nl = nodelist(i) ;
         --x       for (Index_t j=0; j < 8; ++j) {
         --x 	++(nodeElemCount[nl[j]] );
         --x       }
         --x     }
         for Element in This.Elements'Range loop
            for Enode in NodesPerElement_Range loop
               declare
                  Necount : Element_Index renames
                    NodeElemCounts (This.Elements (Element).Node_Indexes (Enode));
               begin
                  Necount := Necount + 1;
               end;
            end loop;
         end loop;
      end Count_Elements_For_Each_Node;

      procedure Calc_First_Element_For_Each_Node is
      begin
         --x     m_nodeElemStart = new Index_t[numNode()+1] ;
         --x     m_nodeElemStart[0] = 0;
         --x     for (Index_t i=1; i <= numNode(); ++i) {
         --x       m_nodeElemStart[i] =
         --x 	m_nodeElemStart[i-1] + nodeElemCount[i-1] ;
         --x     }
         NodeElemStart := new Node_Element_Index_Array (0 .. This.NumNode);
         NodeElemStart (0) := 0;
         for Node in 1 .. This.NumNode loop
            NodeElemStart (Node) :=
              NodeElemStart (Node - 1) + NodeElemCounts (Node - 1);
         end loop;
      end Calc_First_Element_For_Each_Node;

      procedure Calc_Nodes_Per_Element_Offset_For_Each_Corner is
      begin
         --x     m_nodeElemCornerList = new Index_t[m_nodeElemStart[numNode()]];
         NodeElemCornerList := new Element_Element_Index_Array
           (0 .. NodeElemStart (This.NumNode)-1);
         --x     for (Index_t i=0; i < numNode(); ++i) {
         --x       nodeElemCount[i] = 0;
         --x     }
         ---!! Again?
         NodeElemCounts.all := (others => 0);
         --x     for (Index_t i=0; i < numElem(); ++i) {
         --x       Index_t *nl = nodelist(i) ;
         --x       for (Index_t j=0; j < 8; ++j) {
         --x 	Index_t m = nl[j];
         --x 	Index_t k = i*8 + j ;
         --x 	Index_t offset = m_nodeElemStart[m] + nodeElemCount[m] ;
         --x 	m_nodeElemCornerList[offset] = k;
         --x 	++(nodeElemCount[m]) ;
         --x       }
         --x     }
         for Element in 0 .. This.NumElem - 1 loop
            for Enode in NodesPerElement_Range loop
               declare
                  Node   : constant Node_Index :=
                    This.Elements (Element).Node_Indexes (Enode);
                  Corner_Offset : constant Element_Index :=
                    NodeElemStart (Node) + NodeElemCounts (Node) ;
                  Nodes_Per_Element_Offset      : constant Element_Index :=
                    Element * NODES_PER_ELEMENT + Element_Index (Enode) ;
               begin
                  NodeElemCornerList (Corner_Offset) := Nodes_Per_Element_Offset;
                  ---!! Is it ok that this changes from one reference to ther next?
                  NodeElemCounts (Node) := NodeElemCounts (Node) + 1;
               end;
            end loop;
         end loop;
      end Calc_Nodes_Per_Element_Offset_For_Each_Corner;

      procedure Check_Each_Nodes_Per_Element_Offset is
      begin
         --x     Index_t clSize = m_nodeElemStart[numNode()] ;
         --x     for (Index_t i=0; i < clSize; ++i) {
         --x       Index_t clv = m_nodeElemCornerList[i] ;
         --x       if ((clv < 0) || (clv > numElem()*8)) {
         --x 	fprintf(stderr,
         --x 		"AllocateNodeElemIndexes(): nodeElemCornerList entry out of range!\n");
         --x #if USE_MPI
         --x 	MPI_Abort(MPI_COMM_WORLD, -1);
         --x #else
         --x 	exit(-1);
         --x #endif
         --x       }
         --x     }
         declare
            ClSize     : constant Element_Index :=
              NodeElemStart (This.NumNode);
            Max_ClSize : constant Element_Index :=
              This.NumElem * NODES_PER_ELEMENT;
         begin
            for Element in 0 .. ClSize - 1 loop
               declare
                  Clv : constant Element_Index :=
                    NodeElemCornerList (Element);
               begin
                  if Clv > Max_ClSize then
                     if USE_MPI then
                        MPI.Abortt (MPI.COMM_WORLD, MPI.Errorcode_Type (-1));
                     else
                        raise Coding_Error with
                          "AllocateNodeElemIndexes(): nodeElemCornerList entry out of range!" &
                          "  i:" & Element'Img & " clv:" & Clv'Img;
                     end if;
                  end if;
               end;
            end loop;
         end;
      end Check_Each_Nodes_Per_Element_Offset;

   begin
      --x   if (numthreads > 1) {
      ---     // set up node-centered indexing of elements
      --x     Index_t *nodeElemCount = new Index_t[numNode()] ;
      --x     for (Index_t i=0; i<numNode(); ++i) {
      --x       nodeElemCount[i] = 0 ;
      --x     }
      if Numthreads > 1 then
         Count_Elements_For_Each_Node;
         Calc_First_Element_For_Each_Node;
         Calc_Nodes_Per_Element_Offset_For_Each_Corner;
         Check_Each_Nodes_Per_Element_Offset;
         --x     delete [] nodeElemCount ;
         Release (NodeElemCounts);
         --x   }
         --x   else {
      else
         ---     // These arrays are not used if we're not threaded
         --x     m_nodeElemStart = NULL;
         --x     m_nodeElemCornerList = NULL;
         NodeElemStart      := null;
         NodeElemCornerList := null;
         --x   }
      end if;
      --x }
   end SetupThreadSupportStructures;


   --- ////////////////////////////////////////////////////////////////////////////////
   --x void
   --x Domain::SetupCommBuffers(Int_t edgeNodes)
   --x {
   ---   // allocate a buffer large enough for nodal ghost data
   procedure SetupCommBuffers
     (This      : in out Domain_Record;
      EdgeNodes : in Node_Index)
   is
      TV : Variables_Record  renames This.Variables;
      --x   Index_t maxEdgeSize = MAX(this->sizeX(), MAX(this->sizeY(), this->sizeZ()))+1 ;
      MaxEdgeSize : constant Node_Index :=
        Node_Index (MAX (This.Parameters.Size (X),
                    MAX (This.Parameters.Size (Y),
                      This.Parameters.Size (Z))) + 1);
   begin
      --x   m_maxPlaneSize = CACHE_ALIGN_REAL(maxEdgeSize*maxEdgeSize) ;
      --x   m_maxEdgeSize = CACHE_ALIGN_REAL(maxEdgeSize) ;
      TV.MaxPlaneSize := CACHE_ALIGN_REAL (MaxEdgeSize ** 2) ;
      TV.MaxEdgeSize := CACHE_ALIGN_REAL (MaxEdgeSize) ;
      ---   // assume communication to 6 neighbors by default
      --x   m_rowMin = (m_rowLoc == 0)        ? 0 : 1;
      --x   m_rowMax = (m_rowLoc == m_tp-1)     ? 0 : 1;
      --x   m_colMin = (m_colLoc == 0)        ? 0 : 1;
      --x   m_colMax = (m_colLoc == m_tp-1)     ? 0 : 1;
      --x   m_planeMin = (m_planeLoc == 0)    ? 0 : 1;
      --x   m_planeMax = (m_planeLoc == m_tp-1) ? 0 : 1;
      TV.At_Limit (Row, Min)   := TV.RowLoc = 0;
      TV.At_Limit (Row, Max)   := TV.RowLoc = TV.Tp - 1;
      TV.At_Limit (Col, Min)   := TV.ColLoc = 0;
      TV.At_Limit (Col, Max)   := TV.ColLoc = TV.Tp - 1;
      TV.At_Limit (Plane, Min) := TV.PlaneLoc = 0;
      TV.At_Limit (Plane, Max) := TV.PlaneLoc = TV.Tp - 1;

      --x #if USE_MPI
      if USE_MPI then
         declare
            ---   // account for face communication
            --x   Index_t comBufSize =
            --x     (m_rowMin + m_rowMax + m_colMin + m_colMax + m_planeMin + m_planeMax) *
            --x     m_maxPlaneSize * MAX_FIELDS_PER_MPI_COMM ;
            Plane_Buffer_Size : constant Index_Type :=
              Index_Type (TV.MaxPlaneSize * MAX_FIELDS_PER_MPI_COMM);
            Edge_Buffer_Size : constant Index_Type :=
              Index_Type (TV.MaxEdgeSize * MAX_FIELDS_PER_MPI_COMM);
            Corner_Buffer_Size : constant Index_Type :=
              Index_Type (1 * MAX_FIELDS_PER_MPI_COMM);
            ComBufSize : Index_Type := 0;
         begin
            -- 6 planes:
            for RCP in Row_Col_Plane loop
               for M in Min_Max loop
                  if TV.At_Limit (RCP, M) then
                     ComBufSize := ComBufSize + Plane_Buffer_Size;
                  end if;
               end loop;
            end loop;
            ---   // account for edge communication
            --x   comBufSize +=
            --x     ((m_rowMin & m_colMin) + (m_rowMin & m_planeMin) + (m_colMin & m_planeMin) +
            --x      (m_rowMax & m_colMax) + (m_rowMax & m_planeMax) + (m_colMax & m_planeMax) +
            --x      (m_rowMax & m_colMin) + (m_rowMin & m_planeMax) + (m_colMin & m_planeMax) +
            --x      (m_rowMin & m_colMax) + (m_rowMax & m_planeMin) + (m_colMax & m_planeMin)) *
            --x     m_maxPlaneSize * MAX_FIELDS_PER_MPI_COMM ;
            -- 12 edges:
            for RCP in Row_Col_Plane loop
               for This_Min_Max in Min_Max loop
                  for Succ_Min_Max in Min_Max loop
                     if TV.At_Limit (RCP, This_Min_Max) and then
                       TV.At_Limit (Succ_Wrap (RCP), Succ_Min_Max) then
                        ComBufSize := ComBufSize + Edge_Buffer_Size;
                     end if;
                  end loop;
               end loop;
            end loop;
            ---   // account for corner communication
            ---   // factor of 16 is so each buffer has its own cache line
            --x   comBufSize += ((m_rowMin & m_colMin & m_planeMin) +
            --x 		 (m_rowMin & m_colMin & m_planeMax) +
            --x 		 (m_rowMin & m_colMax & m_planeMin) +
            --x 		 (m_rowMin & m_colMax & m_planeMax) +
            --x 		 (m_rowMax & m_colMin & m_planeMin) +
            --x 		 (m_rowMax & m_colMin & m_planeMax) +
            --x 		 (m_rowMax & m_colMax & m_planeMin) +
            --x 		 (m_rowMax & m_colMax & m_planeMax)) * CACHE_COHERENCE_PAD_REAL ;
            --x   this->commDataSend = new Real_t[comBufSize] ;
            --x   this->commDataRecv = new Real_t[comBufSize] ;
            -- 8 corners:
            for Row_Min_Max in Min_Max loop
               for Col_Min_Max in Min_Max loop
                  for Plane_Min_Max in Min_Max loop
                     if TV.At_Limit (Row, Row_Min_Max) and then
                       TV.At_Limit (Col, Col_Min_Max) and then
                       TV.At_Limit (Plane, Plane_Min_Max) then
                        ComBufSize := ComBufSize + Edge_Buffer_Size;
                     end if;
                  end loop;
               end loop;
            end loop;

            TV.CommDataSend := new MPI.Comm_Buffer (0 .. Natural (ComBufSize) - 1);
            TV.CommDataRecv := new MPI.Comm_Buffer (0 .. Natural (ComBufSize) - 1);
            ---   // prevent floating point exceptions
            --x   memset(this->commDataSend, 0, comBufSize*sizeof(Real_t)) ;
            --x   memset(this->commDataRecv, 0, comBufSize*sizeof(Real_t)) ;
            TV.CommDataSend.all := (others => 0.0);
            TV.CommDataRecv.all := (others => 0.0);
            --x #endif
         end;
      end if;

      ---   // Boundary nodesets
      --x   if (m_colLoc == 0)
      --x     m_symmX.resize(edgeNodes*edgeNodes);
      --x   if (m_rowLoc == 0)
      --x     m_symmY.resize(edgeNodes*edgeNodes);
      --x   if (m_planeLoc == 0)
      --x     m_symmZ.resize(edgeNodes*edgeNodes);
      if TV.ColLoc = 0 then
         TV.Symmetry_Plane_Nodes (X) := new Node_Index_Array (0 .. EdgeNodes ** 2);
      end if;
      if TV.RowLoc = 0 then
         TV.Symmetry_Plane_Nodes (Y) := new Node_Index_Array (0 .. EdgeNodes ** 2);
      end if;
      if TV.PlaneLoc = 0 then
         TV.Symmetry_Plane_Nodes (Z) := new Node_Index_Array (0 .. EdgeNodes ** 2);
      end if;
      --x }
   end SetupCommBuffers;

   --- ////////////////////////////////////////////////////////////////////////////////
   --x void
   --x Domain::CreateRegionIndexSets(Int_t nr, Int_t balance)
   --x {
   procedure CreateRegionIndexSets
     (This    : in out Domain_Record;
      Nreg    : in Region_Index;
      Balance : in Balance_Type)
   is
      NextIndex : Element_Index := 0;
      --x #if USE_MPI
      --x    Index_t myRank;
      --x    MPI_Comm_rank(MPI_COMM_WORLD, &myRank) ;
      --x    srand(myRank);
      --x #else
      --x    srand(0);
      --x    Index_t myRank = 0;
      --x #endif
      MyRank    : MPI.Rank_Type := 0;
   begin
      if USE_MPI then
         MPI.Comm_Rank (MPI.COMM_WORLD, MyRank);
      end if;
      Random_Selection.Initialize (MyRank);

      --x    this->numReg() = nr;
      This.NumReg := Nreg;
      --x    m_regElemSize = new Index_t[numReg()];
      --x    m_regElemlist = new Index_t*[numReg()];
      This.Regions := new Region_Array (0 .. This.NumReg - 1);
      --x    Index_t nextIndex = 0;
      ---    //if we only have one region just fill it
      ---    // Fill out the regNumList with material numbers, which are always
      ---    // the region index plus one
      --x    if(numReg() == 1) {
      --x       while (nextIndex < numElem()) {
      --x 	 this->regNumList(nextIndex) = 1;
      --x          nextIndex++;
      --x       }
      --x       regElemSize(0) = 0;
      --x    }
      if This.NumReg = 1 then
         while (NextIndex < This.NumElem) loop
            This.Elements (NextIndex).Region := 1;
            NextIndex := NextIndex + 1;
         end loop;
         This.Regions (0).Size := 0;
         ---    //If we have more than one region distribute the elements.
         --x    else {
      else
         declare
            --x       Int_t regionNum;
            --x       Int_t regionVar;
            --x       Int_t lastReg = -1;
            --x       Int_t binSize;
            --x       Index_t elements;
            --x       Index_t runto = 0;
            --x       Int_t costDenominator = 0;
            --x       Int_t* regBinEnd = new Int_t[numReg()];
            RegionNum       : Region_Index;
            RegionVar       : Cost_Type;
            LastReg         : Region_Index := Region_Index'Last;
            BinSize         : Int_T;
            Elements        : Element_Index;
            Runto           : Element_Index := 0;
            CostDenominator : Cost_Type := 0;
            RegBinEnd       : constant Region_Bin_End_Array_Access :=
              new Region_Bin_End_Array (0 .. This.NumReg - 1);
         begin
            ---       //Determine the relative weights of all the regions.  This is based off the -b flag.  Balance is the value passed into b.
            --x       for (Index_t i=0 ; i<numReg() ; ++i) {
            --x          regElemSize(i) = 0;
            --x 	 costDenominator += pow((i+1), balance);  //Total sum of all regions weights
            --x 	 regBinEnd[i] = costDenominator;  //Chance of hitting a given region is (regBinEnd[i] - regBinEdn[i-1])/costDenominator
            --x       }
            for Region in This.Regions'Range loop
               This.Regions (Region).Size := 0;
               CostDenominator := CostDenominator + Cost_Type ((Region + 1) ** Natural (Balance));  --- //Total sum of all regions weights
               RegBinEnd (Region) := CostDenominator;  --- //Chance of hitting a given region is (regBinEnd[i] - regBinEdn[i-1])/costDenominator
            end loop;

            ---       //Until all elements are assigned
            --x       while (nextIndex < numElem()) {
            while NextIndex < This.NumElem loop
               --- 	 //pick the region
               --- 	 regionVar = rand() % costDenominator;
               RegionVar := Choose_Rem (CostDenominator);
               declare
                  --x 	 Index_t i = 0;
                  I : Region_Index := 0;
               begin
                  --x          while(regionVar >= regBinEnd[i])
                  --x 	    i++;
                  while RegionVar >= RegBinEnd (I) loop
                     I := I + 1;
                  end loop;
                  ---          //rotate the regions based on MPI rank.
                  --- Rotation is Rank % NumRegions this makes each domain
                  --- have a different region with
                  ---          //the highest representation
                  --x 	 regionNum = ((i + myRank) % numReg()) + 1;
                  --- 	 // make sure we don't pick the same region twice in a row
                  --x          while(regionNum == lastReg) {
                  --x 	    regionVar = rand() % costDenominator;
                  --x 	    i = 0;
                  --x             while(regionVar >= regBinEnd[i])
                  --x 	       i++;
                  --x 	    regionNum = ((i + myRank) % numReg()) + 1;
                  --x          }
                  RegionNum := ((I + Region_Index (MyRank)) rem This.NumReg) + 1;
                  while RegionNum = LastReg loop
                     RegionVar := Choose_Rem (CostDenominator);
                     I := 0;
                     while RegionVar >= RegBinEnd (I) loop
                        I := I + 1;
                     end loop;
                     RegionNum := ((I + Region_Index (MyRank)) rem This.NumReg) + 1;
                  end loop;
               end;
               --- 	 //Pick the bin size of the region and determine the number of elements.
               --x          binSize = rand() % 1000;
               --x 	 if(binSize < 773) {
               --x 	   elements = rand() % 15 + 1;
               --x 	 }
               --x 	 else if(binSize < 937) {
               --x 	   elements = rand() % 16 + 16;
               --x 	 }
               --x 	 else if(binSize < 970) {
               --x 	   elements = rand() % 32 + 32;
               --x 	 }
               --x 	 else if(binSize < 974) {
               --x 	   elements = rand() % 64 + 64;
               --x 	 }
               --x 	 else if(binSize < 978) {
               --x 	   elements = rand() % 128 + 128;
               --x 	 }
               --x 	 else if(binSize < 981) {
               --x 	   elements = rand() % 256 + 256;
               --x 	 }
               --x 	 else
               --x 	    elements = rand() % 1537 + 512;
               BinSize := Choose_Rem (1000);
               if BinSize < 773 then
                  Elements := Choose_Rem (15 + 1);
               elsif (BinSize < 937) then
                  Elements := Choose_Rem (16 + 16);
               elsif (BinSize < 970) then
                  Elements := Choose_Rem (32 + 32);
               elsif (BinSize < 974) then
                  Elements := Choose_Rem (64 + 64);
               elsif (BinSize < 978) then
                  Elements := Choose_Rem (128 + 128);
               elsif (BinSize < 981) then
                  Elements := Choose_Rem (256 + 256);
               else
                  Elements := Choose_Rem (1537 + 512);
               end if;
               --x 	 runto = elements + nextIndex;
               --x 	 //Store the elements.  If we hit the end before we run out of elements then just stop.
               --x          while (nextIndex < runto && nextIndex < numElem()) {
               --x 	    this->regNumList(nextIndex) = regionNum;
               --x 	    nextIndex++;
               --x 	 }
               Runto := Elements + NextIndex;
               while NextIndex < Runto and NextIndex < This.NumElem loop
                  This.Elements (NextIndex).Region := RegionNum;
                  NextIndex := NextIndex + 1;
               end loop;
               --x 	 lastReg = regionNum;
               LastReg := RegionNum;
               --x       }
            end loop;
            --x    }
         end;
      end if;
      ---    // Convert regNumList to region index sets
      ---    // First, count size of each region
      --x    for (Index_t i=0 ; i<numElem() ; ++i) {
      --x       int r = this->regNumList(i)-1; // region index == regnum-1
      --x       regElemSize(r)++;
      --x    }
      for Element in This.Elements'Range loop
         declare
            --// region index == regnum-1
            Region : constant Region_Index :=
              This.Elements (Element).Region - 1;
            Size : Element_Index renames
              This.Regions (Region).Size;
         begin
            Size := Size + 1;
         end;
      end loop;
      ---    // Second, allocate each region index set
      --x    for (Index_t i=0 ; i<numReg() ; ++i) {
      --x       m_regElemlist[i] = new Index_t[regElemSize(i)];
      --x       regElemSize(i) = 0;
      --x    }
      for Region in This.Regions'Range loop
         This.Regions (Region).Elements := new
           Element_Element_Index_Array (0 .. This.Regions (Region).Size - 1);
         This.Regions (Region).Size := 0;
      end loop;
      ---    // Third, fill index sets
      --x    for (Index_t i=0 ; i<numElem() ; ++i) {
      --x       Index_t r = regNumList(i)-1;       // region index == regnum-1
      --x       Index_t regndx = regElemSize(r)++; // Note increment
      --x       regElemlist(r,regndx) = i;
      --x    }
      for Element in This.Elements'Range loop
         declare
            Region              : constant Region_Index :=
              This.Elements (Element).Region - 1;
            Last_Region_Element : Element_Index renames
              This.Regions (Region).Size;
         begin
            Last_Region_Element := Last_Region_Element + 1;
            This.Regions (Region).Elements (Last_Region_Element) := Element;
         end;
      end loop;
      --x }
   end CreateRegionIndexSets;

   --- /////////////////////////////////////////////////////////////
   --x void
   --x Domain::SetupSymmetryPlanes(Int_t edgeNodes)
   --x {
   procedure SetupSymmetryPlanes
     (This      : in out Domain_Record;
      EdgeNodes : in Node_Index)
   is
      --x   Index_t nidx = 0 ;
      Node : Node_Index := 0;
   begin
      --x   for (Index_t i=0; i<edgeNodes; ++i) {
      --x     Index_t planeInc = i*edgeNodes*edgeNodes ;
      --x     Index_t rowInc   = i*edgeNodes ;
      --x     for (Index_t j=0; j<edgeNodes; ++j) {
      --x       if (m_planeLoc == 0) {
      --x 	m_symmZ[nidx] = rowInc   + j ;
      --x       }
      --x       if (m_rowLoc == 0) {
      --x 	m_symmY[nidx] = planeInc + j ;
      --x       }
      --x       if (m_colLoc == 0) {
      --x 	m_symmX[nidx] = planeInc + j*edgeNodes ;
      --x       }
      --x       ++nidx ;
      --x     }
      --x   }
      for I in 0 .. EdgeNodes - 1 loop
         declare
            PlaneInc : constant Node_Index := I * EdgeNodes ** 2;
            RowInc   : constant Node_Index := I * EdgeNodes;
         begin
            for J in 0 .. EdgeNodes - 1 loop
               if This.Variables.PlaneLoc = 0 then
                  This.Variables.Symmetry_Plane_Nodes (Z) (Node) := RowInc + J;
               end if;
               if This.Variables.RowLoc = 0 then
                  This.Variables.Symmetry_Plane_Nodes (Y) (Node) := PlaneInc + J;
               end if;
               if This.Variables.ColLoc = 0 then
                  This.Variables.Symmetry_Plane_Nodes (X) (Node) := PlaneInc + J * EdgeNodes;
               end if;
               Node := Node + 1;
            end loop;
         end;
      end loop;
      --x }
   end SetupSymmetryPlanes;

   --- /////////////////////////////////////////////////////////////
   --x void
   --x Domain::SetupElementConnectivities(Int_t edgeElems)
   --x {
   procedure SetupElementConnectivities
     (This      : in out Domain_Record;
      EdgeElems : in Element_Index) is
   begin
      --x    lxim(0) = 0 ;
      --x    for (Index_t i=1; i<numElem(); ++i) {
      --x       lxim(i)   = i-1 ;
      --x       lxip(i-1) = i ;
      --x    }
      --x    lxip(numElem()-1) = numElem()-1 ;
      This.Elements (0).Connections (Xi, M) := 0;
      for I in 0 .. This.NumElem - 1 loop
         This.Elements (I).Connections (Xi, M) := I - 1;
         This.Elements (I - 1).Connections (Xi, P) := I;
      end loop;
      This.Elements (This.NumElem - 1).Connections (Xi, P) := This.NumElem - 1;
      --x    for (Index_t i=0; i<edgeElems; ++i) {
      --x       letam(i) = i ;
      --x       letap(numElem()-edgeElems+i) = numElem()-edgeElems+i ;
      --x    }
      for I in 0 .. EdgeElems - 1 loop
         This.Elements (I).Connections (Eta, M) := I;
         This.Elements (This.NumElem - EdgeElems + I).Connections (Eta, P) :=
           This.NumElem - EdgeElems + I;
      end loop;
      --x    for (Index_t i=edgeElems; i<numElem(); ++i) {
      --x       letam(i) = i-edgeElems ;
      --x       letap(i-edgeElems) = i ;
      --x    }
      for I in EdgeElems .. This.NumElem - 1 loop
         This.Elements (I).Connections (Eta, M) := I - EdgeElems;
         This.Elements (I - EdgeElems).Connections (Eta, P) := I;
      end loop;
      --x    for (Index_t i=0; i<edgeElems*edgeElems; ++i) {
      --x       lzetam(i) = i ;
      --x       lzetap(numElem()-edgeElems*edgeElems+i) = numElem()-edgeElems*edgeElems+i ;
      --x    }
      for I in 0 .. EdgeElems * EdgeElems - 1 loop
         This.Elements (I).Connections (Zeta, M) := I;
         This.Elements (This.NumElem - EdgeElems * EdgeElems + I).Connections (Zeta, P) :=
           This.NumElem - EdgeElems * EdgeElems + I;
      end loop;
      --x    for (Index_t i=edgeElems*edgeElems; i<numElem(); ++i) {
      --x       lzetam(i) = i - edgeElems*edgeElems ;
      --x       lzetap(i-edgeElems*edgeElems) = i ;
      --x    }
      for I in EdgeElems * EdgeElems .. This.NumElem - 1 loop
         This.Elements (I).Connections (Zeta, M) := I - EdgeElems * EdgeElems;
         This.Elements (I - EdgeElems * EdgeElems).Connections (Zeta, P) := I;
      end loop;
      --x }
   end SetupElementConnectivities;

   --- /////////////////////////////////////////////////////////////
   --x void
   --x Domain::SetupBoundaryConditions(Int_t edgeElems)
   --x {
   --x   Index_t ghostIdx[6] ;  // offsets to ghost locations
   procedure SetupBoundaryConditions
     (This      : in out Domain_Record;
      EdgeElems : in Element_Index)
   is
      GhostIdx : FacesPerNode_Element_Index_Array;
      Pidx     : Element_Index;
      Size     : constant Cartesian_Size_Array := This.Parameters.Size;
   begin
      ---   // set up boundary condition information
      --x   for (Index_t i=0; i<numElem(); ++i) {
      --x      elemBC(i) = Int_t(0) ;
      --x   }
      for I in 0 .. This.NumElem - 1 loop
         This.Elements (I).ElemBC := (others => (others => (others => False)));
      end loop;
      --x   for (Index_t i=0; i<6; ++i) {
      --x     ghostIdx[i] = INT_MIN ;
      --x   }
      GhostIdx := (others => Element_Index'Last);
      --x   Int_t pidx = numElem() ;
      --x   if (m_planeMin != 0) {
      --x     ghostIdx[0] = pidx ;
      --x     pidx += sizeX()*sizeY() ;
      --x   }
      Pidx := This.NumElem;
      if This.Variables.At_Limit (Plane, Min) then
         GhostIdx (0) := Pidx;
         Pidx := Pidx + Size (X) * Size (Y);
      end if;
      --x   if (m_planeMax != 0) {
      --x     ghostIdx[1] = pidx ;
      --x     pidx += sizeX()*sizeY() ;
      --x   }
      if This.Variables.At_Limit (Plane, Max) then
         GhostIdx (1) := Pidx;
         Pidx := Pidx + Size (X) * Size (Y);
      end if;
      --x   if (m_rowMin != 0) {
      --x     ghostIdx[2] = pidx ;
      --x     pidx += sizeX()*sizeZ() ;
      --x   }
      if This.Variables.At_Limit (Row, Min) then
         GhostIdx (2) := Pidx;
         Pidx := Pidx + Size (Z) * Size (X);
      end if;
      --x   if (m_rowMax != 0) {
      --x     ghostIdx[3] = pidx ;
      --x     pidx += sizeX()*sizeZ() ;
      --x   }
      if This.Variables.At_Limit (Row, Max) then
         GhostIdx (3) := Pidx;
         Pidx := Pidx + Size (Z) * Size (X);
      end if;
      --x   if (m_colMin != 0) {
      --x     ghostIdx[4] = pidx ;
      --x     pidx += sizeY()*sizeZ() ;
      --x   }
      if This.Variables.At_Limit (Col, Min) then
         GhostIdx (4) := Pidx;
         Pidx := Pidx + Size (Y) * Size (Z);
      end if;
      --x   if (m_colMax != 0) {
      --x     ghostIdx[5] = pidx ;
      --x   }
      if This.Variables.At_Limit (Col, Max) then
         GhostIdx (5) := Pidx;
      end if;

      ---   // symmetry plane or free surface BCs
      --x   for (Index_t i=0; i<edgeElems; ++i) {
      --x     Index_t planeInc = i*edgeElems*edgeElems ;
      --x     Index_t rowInc   = i*edgeElems ;
      --x     for (Index_t j=0; j<edgeElems; ++j) {
      for I in 0 .. EdgeElems - 1 loop
         declare
            PlaneInc : constant Element_Index := I * EdgeElems * EdgeElems ;
            RowInc   : constant Element_Index := I * EdgeElems ;
            function At_Beginning (Domain : in Domain_Index) return Boolean is
              (Domain = 0);
            function At_End (Domain : in Domain_Index) return Boolean is
              (Domain = This.Variables.Tp - 1);
         begin
            for J in 0 .. EdgeElems - 1 loop
               --x       if (m_planeLoc == 0) {
               --x 	elemBC(rowInc+j) |= ZETA_M_SYMM ;
               --x       }
               --x       else {
               --x 	elemBC(rowInc+j) |= ZETA_M_COMM ;
               --x 	lzetam(rowInc+j) = ghostIdx[0] + rowInc + j ;
               --x       }
               declare
                  Element : Element_Record renames
                    This.Elements (RowInc + J);
               begin
                  if At_Beginning (This.Variables.PlaneLoc) then
                     Element.ElemBC (Zeta, M, Symm) := True;
                  else
                     Element.ElemBC (Zeta, M, Common) := True;
                     Element.Connections (Zeta, M) := GhostIdx (0) + RowInc + J ;
                  end if;
               end;
               --x       if (m_planeLoc == m_tp-1) {
               --x 	elemBC(rowInc+j+numElem()-edgeElems*edgeElems) |=
               --x 	  ZETA_P_FREE;
               --x       }
               --x       else {
               --x 	elemBC(rowInc+j+numElem()-edgeElems*edgeElems) |=
               --x 	  ZETA_P_COMM ;
               --x 	lzetap(rowInc+j+numElem()-edgeElems*edgeElems) =
               --x 	  ghostIdx[1] + rowInc + j ;
               --x       }
               declare
                  Element : Element_Record renames
                    This.Elements (RowInc + J + This.NumElem - EdgeElems * EdgeElems);
               begin
                  if At_End (This.Variables.PlaneLoc) then
                     Element.ElemBC (Zeta, P, Free) := True;
                  else
                     Element.ElemBC (Zeta, P, Common) := True;
                     Element.Connections (Zeta, P) := GhostIdx (1) + RowInc + J ;
                  end if;
               end;
               --x       if (m_rowLoc == 0) {
               --x 	elemBC(planeInc+j) |= ETA_M_SYMM ;
               --x       }
               --x       else {
               --x 	elemBC(planeInc+j) |= ETA_M_COMM ;
               --x 	letam(planeInc+j) = ghostIdx[2] + rowInc + j ;
               --x       }
               declare
                  Element : Element_Record renames
                    This.Elements (PlaneInc + J);
               begin
                  if At_Beginning (This.Variables.RowLoc) then
                     Element.ElemBC (Eta, M, Symm) := True;
                  else
                     Element.ElemBC (Eta, M, Common) := True;
                     Element.Connections (Eta, M) := GhostIdx (2) + RowInc + J ;
                  end if;
               end;
               --x       if (m_rowLoc == m_tp-1) {
               --x 	elemBC(planeInc+j+edgeElems*edgeElems-edgeElems) |=
               --x 	  ETA_P_FREE ;
               --x       }
               --x       else {
               --x 	elemBC(planeInc+j+edgeElems*edgeElems-edgeElems) |=
               --x 	  ETA_P_COMM ;
               --x 	letap(planeInc+j+edgeElems*edgeElems-edgeElems) =
               --x 	  ghostIdx[3] +  rowInc + j ;
               --x       }
               declare
                  Element : Element_Record renames
                    This.Elements (PlaneInc + J + EdgeElems * EdgeElems - EdgeElems);
               begin
                  if At_End (This.Variables.RowLoc) then
                     Element.ElemBC (Eta, P, Free) := True;
                  else
                     Element.ElemBC (Eta, P, Common) := True;
                     Element.Connections (Eta, P) := GhostIdx (3) + RowInc + J ;
                  end if;
               end;
               --x       if (m_colLoc == 0) {
               --x 	elemBC(planeInc+j*edgeElems) |= XI_M_SYMM ;
               --x       }
               --x       else {
               --x 	elemBC(planeInc+j*edgeElems) |= XI_M_COMM ;
               --x 	lxim(planeInc+j*edgeElems) = ghostIdx[4] + rowInc + j ;
               --x       }
               declare
                  Element : Element_Record renames
                    This.Elements (PlaneInc + J * EdgeElems);
               begin
                  if At_Beginning (This.Variables.ColLoc) then
                     Element.ElemBC (Xi, M, Symm) := True;
                  else
                     Element.ElemBC (Xi, M, Common) := True;
                     Element.Connections (Eta, M) := GhostIdx (4) + RowInc + J ;
                  end if;
               end;
               --x       if (m_colLoc == m_tp-1) {
               --x 	elemBC(planeInc+j*edgeElems+edgeElems-1) |= XI_P_FREE ;
               --x       }
               --x       else {
               --x 	elemBC(planeInc+j*edgeElems+edgeElems-1) |= XI_P_COMM ;
               --x 	lxip(planeInc+j*edgeElems+edgeElems-1) =
               --x 	  ghostIdx[5] + rowInc + j ;
               --x       }
               declare
                  Element : Element_Record renames
                    This.Elements (PlaneInc + J * EdgeElems + EdgeElems - 1);
               begin
                  if At_End (This.Variables.ColLoc) then
                     Element.ElemBC (Xi, P, Free) := True;
                  else
                     Element.ElemBC (Xi, P, Common) := True;
                     Element.Connections (Xi, P) := GhostIdx (5) + RowInc + J ;
                  end if;
               end;
               --x     }
            end loop;
            --x   }
         end;
      end loop;
      --x }
   end SetupBoundaryConditions;

   --- ///////////////////////////////////////////////////////////////////////////
   --x void InitMeshDecomp(Int_t numRanks, Int_t myRank,
   --x                     Int_t *col, Int_t *row, Int_t *plane, Int_t *side)
   --x {
   procedure InitMeshDecomp
     (NumRanks         : in Rank_Count_Range;
      MyRank           : in Rank_Type;
      Domain_Column    : out Domain_Index;
      Domain_Row       : out Domain_Index;
      Domain_Plane     : out Domain_Index;
      Domains_Per_Side : out Domain_Index)
   is
      --x    Int_t testProcs;
      --x    Int_t dx, dy, dz;
      --x    Int_t myDom;
      Ranks_Root : Rank_Count_Range;
      Dx         : Domain_Index;
      Dy         : Domain_Index;
      Dz         : Domain_Index;
      Dxyz       : Domain_Index;
      My_Domain  : Domain_Index;
   begin
      ---    // Assume cube processor layout for now
      --x    testProcs = Int_t(cbrt(Real_t(numRanks))+0.5) ;
      --x    if (testProcs*testProcs*testProcs != numRanks) {
      --x       printf("Num processors must be a cube of an integer (1, 8, 27, ...)\n") ;
      --z #if USE_MPI
      --z       MPI_Abort(MPI_COMM_WORLD, -1) ;
      --z #else
      --x       exit(-1);
      --z #endif
      --x    }
      --!! same as c++?:
      Ranks_Root := Rank_Type (Cbrt (Real10 (NumRanks)));
      if Ranks_Root ** 3 /= NumRanks then
         Abort_Or_Raise
           (Usage_Error'Identity,
            "Num processors (" & NumRanks'Img &
              ") must be a cube of an integer (1, 8, 27, ...)");
      end if;
      --x    if (sizeof(Real_t) != 4 && sizeof(Real_t) != 8) {
      --x       printf("MPI operations only support float and double right now...\n");
      --x #if USE_MPI
      --x       MPI_Abort(MPI_COMM_WORLD, -1) ;
      --x #else
      --x       exit(-1);
      --x #endif
      --x    }
      if Real_Type'Size /= 4 and Real_Type'Size /= 8 then
         Abort_Or_Raise
           (Usage_Error'Identity,
            "MPI operations only support float and double right now." &
              "  Size of Real_t is " & Real_Type'Size'Img);
      end if;
      --x    if (MAX_FIELDS_PER_MPI_COMM > CACHE_COHERENCE_PAD_REAL) {
      --x       printf("corner element comm buffers too small.  Fix code.\n") ;
      -- #if USE_MPI
      --       MPI_Abort(MPI_COMM_WORLD, -1) ;
      -- #else
      --x       exit(-1);
      -- #endif
      --x    }
      if MAX_FIELDS_PER_MPI_COMM > CACHE_COHERENCE_PAD_REAL then
         raise Coding_Error
           with "corner element comm buffers too small.  Fix code.";
      end if;
      --x    dx = testProcs ;
      --x    dy = testProcs ;
      --x    dz = testProcs ;
      ---    // temporary test
      --x    if (dx*dy*dz != numRanks) {
      --x       printf("error -- must have as many domains as procs\n") ;
      -- #if USE_MPI
      --       MPI_Abort(MPI_COMM_WORLD, -1) ;
      -- #else
      --x       exit(-1);
      -- #endif
      --x    }
      Dx := Domain_Index (Ranks_Root);
      Dy := Domain_Index (Ranks_Root);
      Dz := Domain_Index (Ranks_Root);
      Dxyz := Dx * Dy * Dz;
      if Process_Count_Range (Dxyz) /= NumRanks then
         raise Usage_Error with
           "Domain count:" & Dxyz'Img &
           " does nor equal proc count:" & NumRanks'Img;
      end if;
      --x    Int_t remainder = dx*dy*dz % numRanks ;
      --x    if (myRank < remainder) {
      --x       myDom = myRank*( 1+ (dx*dy*dz / numRanks)) ;
      --x    }
      --x    else {
      --x       myDom = remainder*( 1+ (dx*dy*dz / numRanks)) +
      --x          (myRank - remainder)*(dx*dy*dz/numRanks) ;
      --x    }
      declare
         -- Always 0, since dxyz = numRanks:
         Remainder : constant Rank_Count_Range := Rank_Count_Range (Dxyz) rem NumRanks;
         -- Always 1, since dxyz = numRanks:
         Quotient  : constant Rank_Count_Range := Rank_Count_Range (Dxyz) / NumRanks;
      begin
         -- Always false:
         if MyRank < Remainder then
            My_Domain := Domain_Index (MyRank * (1 + Quotient)) ;
         else
            -- my_domain := 0 * (1 + 1) + (myRank - 0) * 1:
            -- my_domain := myRank:
            My_Domain := Domain_Index (Remainder * (1 + Quotient) +
                                       (MyRank - Remainder) * Quotient);
         end if;
      end;
      --x    *col = myDom % dx ;
      --x    *row = (myDom / dx) % dy ;
      --x    *plane = myDom / (dx*dy) ;
      --x    *side = testProcs;
      Domain_Column    := My_Domain rem Dx;
      Domain_Row       := (My_Domain / Dx) rem Dy;
      Domain_Plane     := My_Domain / (Dx * Dy);
      Domains_Per_Side := Domain_Index (Ranks_Root);
      --x    return;
      --x }
   end InitMeshDecomp;

end LULESH.Init;
