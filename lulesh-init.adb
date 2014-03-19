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

with LULESH.Par;
with Ada.Numerics.Float_Random;

package body LULESH.Init is

   package Random_Selection is

      -- Returns an integer between Min and Max:
      function Choose
        (Min : in Integer;
         Max : in integer)
         return Integer
        with
          post => (Choose'Result in Min..Max);

      -- Returns an integer between 0 and Modulo - 1:
      function Choose_Rem
        (Modulo : in Integer)
         return Integer is
        (Choose (0, Modulo-1));

   end Random_Selection;

   package body Random_Selection is
      package ANFR renames Ada.Numerics.Float_Random;
      Generator : ANFR.Generator;
      function Choose
        (Min : in Integer;
         Max : in integer)
         return Integer is
      begin
         return Integer(Float(ANFR.Random (Generator)) * Float (Max - Min)) + Min;
      end Choose;
   begin
      ANFR.Reset(Generator);
   end Random_Selection;

   function Choose_Rem
     (Modulo : in Integer)
      return Element_Index is
     (Element_Index(Random_Selection.Choose_Rem(Modulo)));

   function Choose_Rem
     (Modulo : in Integer)
      return Int_t is
     (Int_t(Random_Selection.Choose_Rem(Modulo)));

   function Choose_Rem
     (Modulo : in Cost_Type)
      return Cost_Type is
     (Cost_Type(Random_Selection.Choose_Rem(Integer(Modulo))));

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

   --    void AllocateGradients(Int_t numElem, Int_t allElem)
   --    {
   --       // Position gradients
   --       m_delx_xi.resize(numElem) ;
   --       m_delx_eta.resize(numElem) ;
   --       m_delx_zeta.resize(numElem) ;

   --       // Velocity gradients
   --       m_delv_xi.resize(allElem) ;
   --       m_delv_eta.resize(allElem);
   --       m_delv_zeta.resize(allElem) ;
   --    }

   --    void DeallocateGradients()
   --    {
   --       m_delx_zeta.clear() ;
   --       m_delx_eta.clear() ;
   --       m_delx_xi.clear() ;

   --       m_delv_zeta.clear() ;
   --       m_delv_eta.clear() ;
   --       m_delv_xi.clear() ;
   --    }

   --    void AllocateStrains(Int_t numElem)
   --    {
   --       m_dxx.resize(numElem) ;
   --       m_dyy.resize(numElem) ;
   --       m_dzz.resize(numElem) ;
   --    }

   --    void DeallocateStrains()
   --    {
   --       m_dzz.clear() ;
   --       m_dyy.clear() ;
   --       m_dxx.clear() ;
   --    }
   --- /////////////////////////////////////////////////////////////////////
   --x Domain::Domain(Int_t numRanks, Index_t colLoc,
   --x                Index_t rowLoc, Index_t planeLoc,
   --x                Index_t nx, int tp, int nr, int balance, Int_t cost)
   --x    :
   ------------
   -- EXPORTED:
   ------------
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
      return Domain_Record
   is
      this : Domain_Record;
      --x    Index_t edgeElems = nx ;
      --x    Index_t edgeNodes = edgeElems+1 ;
      edgeElems : constant Element_Index := side_length;
      edgeNodes : constant Node_Index := Node_Index(side_length) + 1;
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
      this.parameters :=
        (energy_tolerance           => 1.0e-7,
         pressure_tolerance         => 1.0e-7,
         artificial_viscosity_tolerance => 1.0e-7,
         volume_relative_tolerance  => 1.0e-10,
         velocity_tolerance         => 1.0e-7,
         hgcoef                     => 3.0,
         four_thirds                => 4.0/3.0,
         qstop                      => 1.0e+12,
         monoq_max_slope            => 1.0,
         monoq_limiter_mult         => 2.0,
         qlc_monoq                  => 0.5,
         qqc_monoq                  => 2.0/3.0,
         qqc                        => 2.0,
         eosvmax                    => 1.0e+9,
         eosvmin                    => 1.0e-9,
         pressure_floor             => 0.0,
         energy_floor               => -1.0e+15,
         dvovmax                     => 0.1,
         reference_density          => 1.0,
         --x    this->cost() = cost;
         imbalance_cost             => cost,
         --x    m_sizeX = edgeElems ;
         --x    m_sizeY = edgeElems ;
         --x    m_sizeZ = edgeElems ;
         size                       => (others => edgeElems));
      --x    m_tp       = tp ;
      this.variables.tp := tp;
      --x    m_numRanks = numRanks ;
      this.variables.numRanks := numRanks;

      ---    ///////////////////////////////
      ---    //   Initialize Sedov Mesh
      ---    ///////////////////////////////

      ---    // construct a uniform box for this processor

      --x    m_colLoc   =   colLoc ;
      --x    m_rowLoc   =   rowLoc ;
      --x    m_planeLoc = planeLoc ;
      this.variables.colLoc   := colLoc;
      this.variables.rowLoc   := rowLoc;
      this.variables.planeLoc := planeLoc;
      --x    m_numElem = edgeElems*edgeElems*edgeElems ;
      --x    m_numNode = edgeNodes*edgeNodes*edgeNodes ;
      this.numElem := edgeElems**3;
      this.numNode := edgeNodes**3;

      --x    m_regNumList = new Index_t[numElem()] ;  // material indexset

      --x    // Elem-centered
      --x    AllocateElemPersistent(numElem()) ;
      this.elements := new Element_Array (0..this.numElem-1);

      --x    // Node-centered
      --x    AllocateNodePersistent(numNode()) ;
      this.nodes := new Node_Array (0..this.numNode-1);

      --x    SetupCommBuffers(edgeNodes);
      SetupCommBuffers(this, edgeNodes);

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
      for index in 0..this.numElem-1 loop
         this.elements(index).eenergy              := 0.0;
         this.elements(index).epressure             := 0.0;
         this.elements(index).artificial_viscosity := 0.0;
         this.elements(index).sound_speed          := 0.0;
         this.elements(index).volume_relative      := 1.0;
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
      for index in 0..this.numNode-1 loop
         this.nodes(index).velocity     := (others => 0.0);
         this.nodes(index).acceleration := (others => 0.0);
         this.nodes(index).nmass := 0.0;
      end loop;

      --x    BuildMesh(nx, edgeNodes, edgeElems);
      BuildMesh(this, side_length, edgeNodes, edgeElems);

      -- #if _OPENMP
      --    SetupThreadSupportStructures();
      -- #else
      --    // These arrays are not used if we're not threaded
      --    m_nodeElemStart = NULL;
      --    m_nodeElemCornerList = NULL;
      -- #endif
      ---    // Setup region index sets. For now, these are constant sized
      ---    // throughout the run, but could be changed every cycle to
      ---    // simulate effects of ALE on the lagrange solver
      --x    CreateRegionIndexSets(nr, balance);
      CreateRegionIndexSets(this, nr, balance);
      ---    // Setup symmetry nodesets
      --x    SetupSymmetryPlanes(edgeNodes);
      SetupSymmetryPlanes(this, edgeNodes);
      ---    // Setup element connectivities
      --x    SetupElementConnectivities(edgeElems);
      SetupElementConnectivities(this,edgeElems);
      ---    // Setup symmetry planes and free surface boundary arrays
      --x    SetupBoundaryConditions(edgeElems);
      SetupBoundaryConditions(this, edgeElems);

      ---    // Setup defaults

      ---    // These can be changed (requires recompile) if you want to run
      ---    // with a fixed timestep, or to a different end time, but it's
      ---    // probably easier/better to just run a fixed number of timesteps
      ---    // using the -i flag in 2.x

      --x    dtfixed() = Real_t(-1.0e-6) ; // Negative means use courant condition
      --x    stoptime()  = Real_t(1.0e-2); // *Real_t(edgeElems*tp/45.0) ;
      this.variables.dtfixed  := ART.Time_Span_Last;
      this.variables.use_courant_condition  := True;
      this.variables.stoptime := ART."+"
        (ART.Time_First, ART.To_Time_Span(Duration(1.0e-2)));

      ---    // Initial conditions
      --x    deltatimemultlb() = Real_t(1.1) ;
      --x    deltatimemultub() = Real_t(1.2) ;
      --x    dtcourant() = Real_t(1.0e+20) ;
      --x    dthydro()   = Real_t(1.0e+20) ;
      --x    dtmax()     = Real_t(1.0e-2) ;
      --x    time()    = Real_t(0.) ;
      --x    cycle()   = Int_t(0) ;
      this.variables.delta_time_multiplier_lower_bound := 1.1;
      this.variables.delta_time_multiplier_upper_bound := 1.2;
      this.variables.dtcourant       := ART.Time_Span_Last;
      this.variables.dthydro         := ART.Time_Span_Last;
      this.variables.dtmax           := ART.To_Time_Span(Duration(1.0e-2));
      this.variables.current_time    := ART.Time_First;
      this.variables.cycle           := 0;

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
      for element in 0..this.numElem-1 loop
         declare
            local_coords : NodesPerElement_Coordinate_Array;
            elemToNode   : constant NodesPerElement_Element_Index_Array :=
              this.elements(element).node_indexes;
         begin
            for node in NodesPerElement_Range loop
               local_coords(node) :=
                 this.nodes(elemToNode(node)).coordinate;
            end loop;
            declare
               element_volume : constant Volume :=
                 LULESH.Par.CalcElemVolume (local_coords);
               mass_share     : constant Mass :=
                 Mass(element_volume) / Mass(NODES_PER_ELEMENT);
            begin
               this.elements(element).volume_reference := element_volume;
               this.elements(element).emass := Mass(element_volume);
               for node in NodesPerElement_Range loop
                  declare
                     node_mass : Mass renames
                       this.nodes(elemToNode(node)).nmass;
                  begin
                     node_mass := node_mass + mass_share;
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
         ebase : constant Energy := 3.948746e+7;
         scale : constant Real_Type   :=
           Real_Type(side_length)*Real_Type(this.variables.tp)/45.0;
         einit : constant Energy := ebase*Energy(scale**3);
      begin
         if (this.variables.rowLoc
             + this.variables.colLoc
             + this.variables.planeLoc = 0) then
            this.elements(0).eenergy := einit;
         end if;
         ---    //set initial deltatime base on analytic CFL calculation
         --x    deltatime() = (Real_t(.5)*cbrt(volo(0)))/sqrt(Real_t(2.0)*einit);
         this.variables.deltatime := ART.To_Time_Span
           (Duration
              ((0.5*cbrt(real10(this.elements(0).volume_reference)))/
                   sqrt(2.0*real10(einit))));
      end;
      --x } // End constructor
      return this;
   end Create;


   -- ////////////////////////////////////////////////////////////////////////////////
   -- void
   -- Domain::BuildMesh(Int_t nx, Int_t edgeNodes, Int_t edgeElems)
   -- {
   procedure BuildMesh
     (this        : in out Domain_Record;
      side_length : in Element_Index;
      edgeNodes   : in Node_Index;
      edgeElems   : in Element_Index)
   is
      subtype Edge_Nodes_Range is Node_Index range 0..edgeNodes-1;
      --x   Index_t meshEdgeElems = m_tp*nx ;
      ---   // initialize nodal coordinates
      --x   Index_t nidx = 0 ;
      meshEdgeElems : constant Element_Index :=
        Element_Index(this.variables.tp)*side_length ;
      node          : Node_Index;
      t             : Coordinate_Vector;
      zidx          : Element_Index;

      function Calc_T_Part
        (loc  : in Domain_Index;
         node : in Node_Index)
         return Length is
        (1.125 *
           Length
             ((Element_Index(loc)*side_length)+
                  Element_Index(node)/meshEdgeElems))
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
      node := 0;
      for plane in Edge_Nodes_Range loop
         t(Z) := Calc_T_Part(this.variables.planeLoc, plane);
         for row in Edge_Nodes_Range loop
            t(Y) := Calc_T_Part(this.variables.rowLoc, row);
            for col in Edge_Nodes_Range loop
               t(X) := Calc_T_Part(this.variables.colLoc, col);
               this.nodes(node).coordinate := T;
               node := node + 1;
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
      zidx := 0 ;
      node := 0;
      for plane in Edge_Nodes_Range loop
         for row in Edge_Nodes_Range loop
            for col in Edge_Nodes_Range loop
               this.elements(zidx).node_indexes :=
                 (0 => node                                       ,
                  1 => node                                   + 1 ,
                  2 => node                       + edgeNodes + 1 ,
                  3 => node                       + edgeNodes     ,
                  4 => node + edgeNodes**2                 ,
                  5 => node + edgeNodes**2             + 1 ,
                  6 => node + edgeNodes**2 + edgeNodes + 1 ,
                  7 => node + edgeNodes**2 + edgeNodes     );
               zidx := zidx + 1;
               node := node + 1;
            end loop;
            node := node + 1;
         end loop;
         node := node + edgeNodes ;
      end loop;
      -- }
   end BuildMesh;

   --x ////////////////////////////////////////////////////////////////////////////////
   --x void
   --x Domain::SetupThreadSupportStructures()
   --x {
   procedure SetupThreadSupportStructures
     (this : in out Domain_Record)
   is
   -- #if _OPENMP
   --    Index_t numthreads = omp_get_max_threads();
   -- #else
   --x    Index_t numthreads = 1;
   -- #endif
      numthreads         : Index_Type := 1;
      --- Number of elements each node is in. Interior nodes are in 8 elements,
      --- faces 4, edges 2, and corners 1:
      nodeElemCounts     : Node_Element_Index_Array_Access;
      nodeElemStart      : Node_Element_Index_Array_Access renames
        this.variables.nodeElemStart;
      nodeElemCornerList : Element_Element_Index_Array_Access renames
        this.variables.nodeElemCornerList;

      procedure Count_Elements_For_Each_Node is
      begin
         nodeElemCounts := new Node_Element_Index_Array (0..this.numNode-1);
         nodeElemCounts.all := (others =>0);
         --x     for (Index_t i=0; i<numElem(); ++i) {
         --x       Index_t *nl = nodelist(i) ;
         --x       for (Index_t j=0; j < 8; ++j) {
         --x 	++(nodeElemCount[nl[j]] );
         --x       }
         --x     }
         for element in this.elements'Range loop
            for enode in NodesPerElement_Range loop
               declare
                  necount : Element_Index renames
                    nodeElemCounts(this.elements(element).node_indexes(enode));
               begin
                  necount := necount + 1;
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
         nodeElemStart := new Node_Element_Index_Array (0..this.numNode);
         nodeElemStart (0) := 0;
         for node in 1..this.numNode loop
            nodeElemStart (node) :=
              nodeElemStart (node-1) + nodeElemCounts (node-1);
         end loop;
      end Calc_First_Element_For_Each_Node;

      procedure Calc_Nodes_Per_Element_Offset_For_Each_Corner is
      begin
         --x     m_nodeElemCornerList = new Index_t[m_nodeElemStart[numNode()]];
         nodeElemCornerList := new Element_Element_Index_Array
           (0..nodeElemStart(this.numNode)-1);
         --x     for (Index_t i=0; i < numNode(); ++i) {
         --x       nodeElemCount[i] = 0;
         --x     }
         ---!! Again?
         nodeElemCounts.all := (others =>0);
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
         for element in 0..this.numElem-1 loop
            for enode in NodesPerElement_Range loop
               declare
                  node   : constant Node_Index :=
                    this.elements(element).node_indexes(enode);
                  corner_offset : constant Element_Index :=
                    nodeElemStart(node) + nodeElemCounts(node) ;
                  Nodes_Per_Element_Offset      : constant Element_Index :=
                    element * NODES_PER_ELEMENT + Element_Index(enode) ;
               begin
                  nodeElemCornerList(corner_offset) := Nodes_Per_Element_Offset;
                  ---!! Is it ok that this changes from one reference to ther next?
                  nodeElemCounts(node) := nodeElemCounts(node) + 1;
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
         -- #if USE_MPI
         -- 	MPI_Abort(MPI_COMM_WORLD, -1);
         -- #else
         --x 	exit(-1);
         -- #endif
         --x       }
         --x     }
         declare
            clSize     : constant Element_Index :=
              nodeElemStart(this.numNode);
            max_clSize : constant Element_Index :=
              this.numElem * NODES_PER_ELEMENT;
         begin
            for element in 0..clSize-1 loop
               declare
                  clv : constant Element_Index :=
                    nodeElemCornerList(element);
               begin
                  if not (clv in 0 .. max_clSize) then
                     raise Coding_Error with
                       "AllocateNodeElemIndexes(): nodeElemCornerList entry out of range!"&
                       "  i:" & element'Img & " clv:" & clv'Img;
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
      if numthreads > 1 then
         Count_Elements_For_Each_Node;
         Calc_First_Element_For_Each_Node;
         Calc_Nodes_Per_Element_Offset_For_Each_Corner;
         Check_Each_Nodes_Per_Element_Offset;
         --x     delete [] nodeElemCount ;
         Release (nodeElemCounts);
         --x   }
         --x   else {
      else
         ---     // These arrays are not used if we're not threaded
         --x     m_nodeElemStart = NULL;
         --x     m_nodeElemCornerList = NULL;
         nodeElemStart      := null;
         nodeElemCornerList := null;
         --x   }
      end if;
      --x }
   end SetupThreadSupportStructures;


   -- ////////////////////////////////////////////////////////////////////////////////
   -- void
   -- Domain::SetupCommBuffers(Int_t edgeNodes)
   -- {
   --   // allocate a buffer large enough for nodal ghost data
   --   Index_t maxEdgeSize = MAX(this->sizeX(), MAX(this->sizeY(), this->sizeZ()))+1 ;
   --   m_maxPlaneSize = CACHE_ALIGN_REAL(maxEdgeSize*maxEdgeSize) ;
   --   m_maxEdgeSize = CACHE_ALIGN_REAL(maxEdgeSize) ;

   --   // assume communication to 6 neighbors by default
   --   m_rowMin = (m_rowLoc == 0)        ? 0 : 1;
   --   m_rowMax = (m_rowLoc == m_tp-1)     ? 0 : 1;
   --   m_colMin = (m_colLoc == 0)        ? 0 : 1;
   --   m_colMax = (m_colLoc == m_tp-1)     ? 0 : 1;
   --   m_planeMin = (m_planeLoc == 0)    ? 0 : 1;
   --   m_planeMax = (m_planeLoc == m_tp-1) ? 0 : 1;

   -- #if USE_MPI
   --   // account for face communication
   --   Index_t comBufSize =
   --     (m_rowMin + m_rowMax + m_colMin + m_colMax + m_planeMin + m_planeMax) *
   --     m_maxPlaneSize * MAX_FIELDS_PER_MPI_COMM ;

   --   // account for edge communication
   --   comBufSize +=
   --     ((m_rowMin & m_colMin) + (m_rowMin & m_planeMin) + (m_colMin & m_planeMin) +
   --      (m_rowMax & m_colMax) + (m_rowMax & m_planeMax) + (m_colMax & m_planeMax) +
   --      (m_rowMax & m_colMin) + (m_rowMin & m_planeMax) + (m_colMin & m_planeMax) +
   --      (m_rowMin & m_colMax) + (m_rowMax & m_planeMin) + (m_colMax & m_planeMin)) *
   --     m_maxPlaneSize * MAX_FIELDS_PER_MPI_COMM ;

   --   // account for corner communication
   --   // factor of 16 is so each buffer has its own cache line
   --   comBufSize += ((m_rowMin & m_colMin & m_planeMin) +
   -- 		 (m_rowMin & m_colMin & m_planeMax) +
   -- 		 (m_rowMin & m_colMax & m_planeMin) +
   -- 		 (m_rowMin & m_colMax & m_planeMax) +
   -- 		 (m_rowMax & m_colMin & m_planeMin) +
   -- 		 (m_rowMax & m_colMin & m_planeMax) +
   -- 		 (m_rowMax & m_colMax & m_planeMin) +
   -- 		 (m_rowMax & m_colMax & m_planeMax)) * CACHE_COHERENCE_PAD_REAL ;

   --   this->commDataSend = new Real_t[comBufSize] ;
   --   this->commDataRecv = new Real_t[comBufSize] ;
   --   // prevent floating point exceptions
   --   memset(this->commDataSend, 0, comBufSize*sizeof(Real_t)) ;
   --   memset(this->commDataRecv, 0, comBufSize*sizeof(Real_t)) ;
   -- #endif

   --   // Boundary nodesets
   --   if (m_colLoc == 0)
   --     m_symmX.resize(edgeNodes*edgeNodes);
   --   if (m_rowLoc == 0)
   --     m_symmY.resize(edgeNodes*edgeNodes);
   --   if (m_planeLoc == 0)
   --     m_symmZ.resize(edgeNodes*edgeNodes);
   -- }

   -- ////////////////////////////////////////////////////////////////////////////////
   -- void
   -- Domain::CreateRegionIndexSets(Int_t nr, Int_t balance)
   -- {
   procedure CreateRegionIndexSets
     (this    : in out Domain_Record;
      nreg    : in Region_Index;
      balance : in Balance_Type)
   is
      -- #if USE_MPI
      --    Index_t myRank;
      --    MPI_Comm_rank(MPI_COMM_WORLD, &myRank) ;
      --    srand(myRank);
      -- #else
      --    srand(0);
      --x    Index_t myRank = 0;
      -- #endif
      myRank   : Index_Type := 0;
      nextIndex : Element_Index := 0;
   begin
      --x    this->numReg() = nr;
      this.numReg := nreg;
      --x    m_regElemSize = new Index_t[numReg()];
      --x    m_regElemlist = new Index_t*[numReg()];
      this.regions := new Region_Array (0..this.numReg-1);
      --x    Index_t nextIndex = 0;
      ---    //if we only have one region just fill it
      ---    // Fill out the regNumList with material numbers, which are always
      ---    // the region index plus one
      --x    if(numReg() == 1) {
      --x       while (nextIndex < numElem()) {
      --x 	 this->regNumList(nextIndex) = 1;
      --x          nextIndex++;
      --x       }
      --       regElemSize(0) = 0;
      --x    }
      if this.numReg = 1 then
         while (nextIndex < this.numElem) loop
            this.elements (nextIndex).region_number := 1;
            nextIndex := nextIndex + 1;
         end loop;
         this.regions(0).size :=0;
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
            regionNum       : Region_Index;
            regionVar       : Cost_Type;
            lastReg         : Region_Index := Region_Index'Last;
            binSize         : Int_t;
            elements        : Element_Index;
            runto           : Element_Index := 0;
            costDenominator : Cost_Type := 0;
            regBinEnd       : Region_Bin_End_Array_Access :=
              new Region_Bin_End_Array(0..this.numReg-1);
         begin
            ---       //Determine the relative weights of all the regions.  This is based off the -b flag.  Balance is the value passed into b.
            --x       for (Index_t i=0 ; i<numReg() ; ++i) {
            --x          regElemSize(i) = 0;
            --x 	 costDenominator += pow((i+1), balance);  //Total sum of all regions weights
            --x 	 regBinEnd[i] = costDenominator;  //Chance of hitting a given region is (regBinEnd[i] - regBinEdn[i-1])/costDenominator
            --x       }
            for region_index in this.regions'Range loop
               this.regions(region_index).size := 0;
               costDenominator := costDenominator + Cost_Type((region_index+1)**Natural(balance));  -- //Total sum of all regions weights
               regBinEnd (region_index) := costDenominator;  -- //Chance of hitting a given region is (regBinEnd[i] - regBinEdn[i-1])/costDenominator
            end loop;

            ---       //Until all elements are assigned
            --x       while (nextIndex < numElem()) {
            while nextIndex < this.numElem loop
               --- 	 //pick the region
               --- 	 regionVar = rand() % costDenominator;
               regionVar := Choose_Rem(costDenominator);
               declare
                  --x 	 Index_t i = 0;
                  i : Region_Index := 0;
               begin
                  --x          while(regionVar >= regBinEnd[i])
                  --x 	    i++;
                  while regionVar >= regBinEnd(i) loop
                     i := i+1;
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
                  regionNum := ((i + Region_Index(myRank)) rem this.numReg) + 1;
                  while regionNum = lastReg loop
                     regionVar := Choose_Rem(costDenominator);
                     i := 0;
                     while regionVar >= regBinEnd(i) loop
                        i := i+1;
                     end loop;
                     regionNum := ((i + Region_Index(myRank)) rem this.numReg) + 1;
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
               binSize := Choose_Rem(1000);
               if binSize < 773 then
                  elements := Choose_Rem(15 + 1);
               elsif(binSize < 937) then
                  elements := Choose_Rem(16 + 16);
               elsif(binSize < 970) then
                  elements := Choose_Rem(32 + 32);
               elsif(binSize < 974) then
                  elements := Choose_Rem(64 + 64);
               elsif(binSize < 978) then
                  elements := Choose_Rem(128 + 128);
               elsif(binSize < 981) then
                  elements := Choose_Rem(256 + 256);
               else
                  elements := Choose_Rem(1537 + 512);
               end if;
               --x 	 runto = elements + nextIndex;
               --x 	 //Store the elements.  If we hit the end before we run out of elements then just stop.
               --x          while (nextIndex < runto && nextIndex < numElem()) {
               --x 	    this->regNumList(nextIndex) = regionNum;
               --x 	    nextIndex++;
               --x 	 }
               runto := elements + nextIndex;
               while nextIndex < runto and nextIndex < this.numElem loop
                  this.elements(nextIndex).region_number := regionNum;
                  nextIndex := nextIndex + 1;
               end loop;
               -- 	 lastReg = regionNum;
               lastReg := regionNum;
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
      for element in this.Elements'Range loop
         declare
            --// region index == regnum-1
            region: constant Region_Index :=
              this.elements(element).region_number-1;
            size : Element_Index renames
              this.regions(region).size;
         begin
            size := size + 1;
         end;
      end loop;
      ---    // Second, allocate each region index set
      --x    for (Index_t i=0 ; i<numReg() ; ++i) {
      --x       m_regElemlist[i] = new Index_t[regElemSize(i)];
      --x       regElemSize(i) = 0;
      --x    }
      for region_index in this.Regions'Range loop
         this.regions(region_index).elements := new
           Element_Element_Index_Array (0..this.regions(region_index).size-1);
         this.regions(region_index).size := 0;
      end loop;
      ---    // Third, fill index sets
      --x    for (Index_t i=0 ; i<numElem() ; ++i) {
      --x       Index_t r = regNumList(i)-1;       // region index == regnum-1
      --x       Index_t regndx = regElemSize(r)++; // Note increment
      --x       regElemlist(r,regndx) = i;
      --x    }
      for element in this.Elements'Range loop
         declare
            --       // region index == regnum-1
            region : Region_Index :=
              this.elements(element).region_number-1;
            element : Element_Index renames
              this.regions(region).size;
         begin
            region := region + 1;
            this.regions(region).elements(element) := element;
         end;
      end loop;
      --x }
   end CreateRegionIndexSets;

   --- /////////////////////////////////////////////////////////////
   --x void
   --x Domain::SetupSymmetryPlanes(Int_t edgeNodes)
   --x {
   procedure SetupSymmetryPlanes
     (this      : in out Domain_Record;
      edgeNodes : in Node_Index)
   is
      --x   Index_t nidx = 0 ;
      nidx : Node_Index := 0;
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
      for i in 0..edgeNodes-1 loop
         declare
            planeInc : constant Node_Index := i*edgeNodes**2;
            rowInc   : constant Node_Index := i*edgeNodes;
         begin
            for j in 0..edgeNodes-1 loop
               if this.variables.planeLoc = 0 then
                  this.nodes(nidx).symmetry_plane_nodes(Z) := rowInc + j;
               end if;
               if this.variables.rowLoc = 0 then
                  this.nodes(nidx).symmetry_plane_nodes(Y) := planeInc + j;
               end if;
               if this.variables.colLoc = 0 then
                  this.nodes(nidx).symmetry_plane_nodes(X) := planeInc + j*edgeNodes;
               end if;
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
     (this      : in out Domain_Record;
      edgeElems : in Element_Index) is
   begin
      --x    lxim(0) = 0 ;
      --x    for (Index_t i=1; i<numElem(); ++i) {
      --x       lxim(i)   = i-1 ;
      --x       lxip(i-1) = i ;
      --x    }
      --x    lxip(numElem()-1) = numElem()-1 ;
      this.elements(0).connections(Xi, M) := 0;
      for i in 0..this.numElem-1 loop
         this.elements(i).connections(Xi, M) := i-1;
         this.elements(i-1).connections(Xi, P) := i;
      end loop;
      this.elements(this.numElem-1).connections(Xi, P) := this.numElem-1;
      --x    for (Index_t i=0; i<edgeElems; ++i) {
      --x       letam(i) = i ;
      --x       letap(numElem()-edgeElems+i) = numElem()-edgeElems+i ;
      --x    }
      for i in 0..edgeElems-1 loop
         this.elements(i).connections(Eta, M) := i;
         this.elements(this.numElem-edgeElems+i).connections(Eta, P) :=
           this.numElem-edgeElems+i;
      end loop;
      --x    for (Index_t i=edgeElems; i<numElem(); ++i) {
      --x       letam(i) = i-edgeElems ;
      --x       letap(i-edgeElems) = i ;
      --x    }
      for i in edgeElems..this.numElem-1 loop
         this.elements(i).connections(Eta, M) := i-edgeElems;
         this.elements(i-edgeElems).connections(Eta, P) := i;
      end loop;
      --x    for (Index_t i=0; i<edgeElems*edgeElems; ++i) {
      --x       lzetam(i) = i ;
      --x       lzetap(numElem()-edgeElems*edgeElems+i) = numElem()-edgeElems*edgeElems+i ;
      --x    }
      for i in 0..edgeElems*edgeElems-1 loop
         this.elements(i).connections(Zeta, M) := i;
         this.elements(this.numElem-edgeElems*edgeElems+i).connections(Zeta, P) :=
           this.numElem-edgeElems*edgeElems+i;
      end loop;
      --x    for (Index_t i=edgeElems*edgeElems; i<numElem(); ++i) {
      --x       lzetam(i) = i - edgeElems*edgeElems ;
      --x       lzetap(i-edgeElems*edgeElems) = i ;
      --x    }
      for i in edgeElems*edgeElems..this.numElem-1 loop
         this.elements(i).connections(Zeta, M) := i - edgeElems*edgeElems;
         this.elements(i-edgeElems*edgeElems).connections(Zeta, P) := i;
      end loop;
      --x }
   end SetupElementConnectivities;

   --- /////////////////////////////////////////////////////////////
   --x void
   --x Domain::SetupBoundaryConditions(Int_t edgeElems)
   --x {
   --x   Index_t ghostIdx[6] ;  // offsets to ghost locations
   procedure SetupBoundaryConditions
     (this      : in out Domain_Record;
      edgeElems : in Element_Index)
   is
      ghostIdx : FacesPerNode_Element_Index_Array;
      pidx     : Element_Index;
      size     : constant Cartesian_Size_Array := this.parameters.size;
   begin
      ---   // set up boundary condition information
      --x   for (Index_t i=0; i<numElem(); ++i) {
      --x      elemBC(i) = Int_t(0) ;
      --x   }
      for i in 0..this.numElem-1 loop
         this.elements(i).elemBC := (others => (others => (others => False)));
      end loop;
      --x   for (Index_t i=0; i<6; ++i) {
      --x     ghostIdx[i] = INT_MIN ;
      --x   }
      ghostIdx := (others => Element_Index'Last);
      --x   Int_t pidx = numElem() ;
      --x   if (m_planeMin != 0) {
      --x     ghostIdx[0] = pidx ;
      --x     pidx += sizeX()*sizeY() ;
      --x   }
      pidx := this.numElem;
      if this.variables.planeMin /= 0 then
         ghostIdx(0) := pidx;
         pidx := pidx + size(X)*size(Y);
      end if;
      --x   if (m_planeMax != 0) {
      --x     ghostIdx[1] = pidx ;
      --x     pidx += sizeX()*sizeY() ;
      --x   }
      if this.variables.planeMax /= 0 then
         ghostIdx(1) := pidx;
         pidx := pidx + size(X)*size(Y);
      end if;
      --x   if (m_rowMin != 0) {
      --x     ghostIdx[2] = pidx ;
      --x     pidx += sizeX()*sizeZ() ;
      --x   }
      if this.variables.rowMin /= 0 then
         ghostIdx(2) := pidx;
         pidx := pidx + size(Z)*size(X);
      end if;
      --x   if (m_rowMax != 0) {
      --x     ghostIdx[3] = pidx ;
      --x     pidx += sizeX()*sizeZ() ;
      --x   }
      if this.variables.rowMax /= 0 then
         ghostIdx(3) := pidx;
         pidx := pidx + size(Z)*size(X);
      end if;
      --x   if (m_colMin != 0) {
      --x     ghostIdx[4] = pidx ;
      --x     pidx += sizeY()*sizeZ() ;
      --x   }
      if this.variables.colMin /= 0 then
         ghostIdx(4) := pidx;
         pidx := pidx + size(Y)*size(Z);
      end if;
      --x   if (m_colMax != 0) {
      --x     ghostIdx[5] = pidx ;
      --x   }
      if this.variables.colMax /= 0 then
         ghostIdx(5) := pidx;
      end if;

      ---   // symmetry plane or free surface BCs
      --x   for (Index_t i=0; i<edgeElems; ++i) {
      --x     Index_t planeInc = i*edgeElems*edgeElems ;
      --x     Index_t rowInc   = i*edgeElems ;
      --x     for (Index_t j=0; j<edgeElems; ++j) {
      for i in 0..edgeElems-1 loop
         declare
            planeInc : constant Element_Index := i*edgeElems*edgeElems ;
            rowInc   : constant Element_Index := i*edgeElems ;
            function At_Beginning (domain : in Domain_Index) return Boolean is
               (domain = 0);
            function At_End (domain : in Domain_Index) return Boolean is
               (domain = this.variables.tp-1);
         begin
            for j in 0..edgeElems-1 loop
               --x       if (m_planeLoc == 0) {
               --x 	elemBC(rowInc+j) |= ZETA_M_SYMM ;
               --x       }
               --x       else {
               --x 	elemBC(rowInc+j) |= ZETA_M_COMM ;
               --x 	lzetam(rowInc+j) = ghostIdx[0] + rowInc + j ;
               --x       }
               declare
                  element : Element_Record renames
                    this.elements(rowInc+j);
               begin
                  if At_Beginning(this.variables.planeLoc) then
                     element.elemBC(Zeta, M, Symm) := True;
                  else
                     element.elemBC(Zeta, M, Common) := True;
                     element.connections(Zeta, M) := ghostIdx(0) + rowInc + j ;
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
                  element : Element_Record renames
                    this.elements(rowInc+j+this.numElem-edgeElems*edgeElems);
               begin
                  if At_End(this.variables.planeLoc) then
                     element.elemBC(Zeta, P, Free) := True;
                  else
                     element.elemBC(Zeta, P, Common) := True;
                     element.connections(Zeta, P) := ghostIdx(1) + rowInc + j ;
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
                  element : Element_Record renames
                    this.elements(planeInc+j);
               begin
                  if At_Beginning(this.variables.rowLoc) then
                     element.elemBC(Eta, M, Symm) := True;
                  else
                     element.elemBC(Eta, M, Common) := True;
                     element.connections(Eta, M) := ghostIdx(2) + rowInc + j ;
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
                  element : Element_Record renames
                    this.elements(planeInc+j+edgeElems*edgeElems-edgeElems);
               begin
                  if At_End(this.variables.rowLoc) then
                     element.elemBC(Eta, P, Free) := True;
                  else
                     element.elemBC(Eta, P, Common) := True;
                     element.connections(Eta, P) := ghostIdx(3) + rowInc + j ;
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
                  element : Element_Record renames
                    this.elements(planeInc+j*edgeElems);
               begin
                  if At_Beginning(this.variables.colLoc) then
                     element.elemBC(Xi, M, Symm) := True;
                  else
                     element.elemBC(Xi, M, Common) := True;
                     element.connections(Eta, M) := ghostIdx(4) + rowInc + j ;
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
                  element : Element_Record renames
                    this.elements(planeInc+j*edgeElems+edgeElems-1);
               begin
                  if At_End(this.variables.colLoc) then
                     element.elemBC(Xi, P, Free) := True;
                  else
                     element.elemBC(Xi, P, Common) := True;
                     element.connections(Xi, P) := ghostIdx(5) + rowInc + j ;
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
     (numRanks         : in Rank_Count_Range;
      myRank           : in Rank_Type;
      domain_column    : out Domain_Index;
      domain_row       : out Domain_Index;
      domain_plane     : out Domain_Index;
      domains_per_side : out Domain_Index)
   is
      --x    Int_t testProcs;
      --x    Int_t dx, dy, dz;
      --x    Int_t myDom;
      ranks_root : Rank_Count_Range;
      dx         : Domain_Index;
      dy         : Domain_Index;
      dz         : Domain_Index;
      dxyz       : Domain_Index;
      my_domain  : Domain_Index;
   begin
      ---    // Assume cube processor layout for now
      --x    testProcs = Int_t(cbrt(Real_t(numRanks))+0.5) ;
      --x    if (testProcs*testProcs*testProcs != numRanks) {
      --x       printf("Num processors must be a cube of an integer (1, 8, 27, ...)\n") ;
      -- #if USE_MPI
      --       MPI_Abort(MPI_COMM_WORLD, -1) ;
      -- #else
      --x       exit(-1);
      -- #endif
      --x    }
      --!! same as c++?:
      ranks_root := Rank_Type(cbrt(real10(numRanks)));
      if ranks_root**3 /= numRanks then
         raise Usage_Error
           with "Num processors (" & numRanks'Img &
           ") must be a cube of an integer (1, 8, 27, ...)";
      end if;
      --x    if (sizeof(Real_t) != 4 && sizeof(Real_t) != 8) {
      --x       printf("MPI operations only support float and double right now...\n");
      -- #if USE_MPI
      --       MPI_Abort(MPI_COMM_WORLD, -1) ;
      -- #else
      --x       exit(-1);
      -- #endif
      --x    }
      if Real_Type'Size /= 4 and Real_Type'Size /= 8 then
         raise Usage_Error
           with "MPI operations only support float and double right now." &
           "  Size of Real_t is " & Real_Type'Size'Img;
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
      dx := Domain_Index(ranks_root);
      dy := Domain_Index(ranks_root);
      dz := Domain_Index(ranks_root);
      dxyz := dx*dy*dz;
      if Process_Count_Range(dxyz) /= numRanks then
         raise Usage_Error with
           "Domain count:" & dxyz'Img &
           " does nor equal proc count:" & numRanks'Img;
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
         remainder : constant Rank_Count_Range := Rank_Count_Range(dxyz) rem numRanks;
         -- Always 1, since dxyz = numRanks:
         quotient  : constant Rank_Count_Range := Rank_Count_Range(dxyz) / numRanks;
      begin
         -- Always false:
         if myRank < remainder then
            my_domain := Domain_Index(myRank * (1 + quotient)) ;
         else
            -- my_domain := 0 * (1 + 1) + (myRank - 0) * 1:
            -- my_domain := myRank:
            my_domain := Domain_Index(remainder * (1 + quotient) +
              (myRank - remainder) * quotient);
         end if;
      end;
      --x    *col = myDom % dx ;
      --x    *row = (myDom / dx) % dy ;
      --x    *plane = myDom / (dx*dy) ;
      --x    *side = testProcs;
      domain_column    := my_domain rem dx;
      domain_row       := (my_domain / dx) rem dy;
      domain_plane     := my_domain / (dx * dy);
      domains_per_side := Domain_Index (ranks_root);
      --x    return;
      --x }
   end InitMeshDecomp;

end LULESH.Init;
