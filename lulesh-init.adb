-- #include <math.h>
-- #if USE_MPI
-- # include <mpi.h>
-- #endif
-- #if _OPENMP
-- #include <omp.h>
-- #endif
-- #include <stdio.h>
-- #include <stdlib.h>
-- #include <string.h>
-- #include <limits.h>
-- #include <cstdlib>
-- #include "lulesh.h"


package body LULESH.Init is

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
     (numRanks : in Int_t;
      colLoc   : in Index_Type;
      rowLoc   : in Index_Type;
      planeLoc : in Index_Type;
      nx       : in Index_Type;
      tp       : in Int_t;
      nr       : in Int_t;
      balance  : in Int_t;
      cost     : in Int_t)
      return Domain_Record
   is
      this : Domain_Record;
      --x    Index_t edgeElems = nx ;
      --x    Index_t edgeNodes = edgeElems+1 ;
      edgeElems : constant Element_Index_Type := Element_Index_Type(nx);
      edgeNode  : constant Node_Index_Type := Node_Index_Type(nx) + 1;
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
         dynamic_pressure_tolerance => 1.0e-7,
         relative_volume_tolerance  => 1.0e-10,
         velocity_tolerance         => 1.0e-7,
         hgcoef                   => 3.0,
         four_thirds                => 4.0/3.0,
         qstop                    => 1.0e+12,
         monoq_max_slope          => 1.0,
         monoq_limiter_mult       => 2.0,
         qlc_monoq                => 0.5,
         qqc_monoq                => 2.0/3.0,
         qqc                      => 2.0,
         eosvmax                  => 1.0e+9,
         eosvmin                  => 1.0e-9,
         pressure_floor             => 0.0,
         energy_floor               => -1.0e+15,
         volume_delta_max           => 0.1,
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
      this.numElem := edgeElems*edgeElems*edgeElems;
      this.numNode := edgeNodes*edgeNodes*edgeNodes;

      --x    m_regNumList = new Index_t[numElem()] ;  // material indexset

      --x    // Elem-centered
      --x    AllocateElemPersistent(numElem()) ;
      this.elements := new Element_Array (0..this.numElem-1);

      --x    // Node-centered
      --x    AllocateNodePersistent(numNode()) ;
      this.nodes := new Node_Array (0..this.numNode-1);

--    SetupCommBuffers(edgeNodes);

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
         this.elements(index).energy           := 0.0;
         this.elements(index).pressure         := 0.0;
         this.elements(index).dynamic_pressure := 0.0;
         this.elements(index).sound_speed      := 0.0;
         this.elements(index).volume           := 1.0;
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
         this.nodes.d := (others => 0.0);
         this.nodes.dd := (others => 0.0);
         this.nodes.mass := 0.0;
      end loop;

      --    BuildMesh(nx, edgeNodes, edgeElems);

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
      CreateRegionIndexSets(nr, balance);
      ---    // Setup symmetry nodesets
      --x    SetupSymmetryPlanes(edgeNodes);
      SetupSymmetryPlanes(edgeNodes);
      ---    // Setup element connectivities
      --x    SetupElementConnectivities(edgeElems);
      SetupElementConnectivities(edgeElems);
      ---    // Setup symmetry planes and free surface boundary arrays
      --x    SetupBoundaryConditions(edgeElems);
      SetupBoundaryConditions(edgeElems);

      ---    // Setup defaults

      ---    // These can be changed (requires recompile) if you want to run
      ---    // with a fixed timestep, or to a different end time, but it's
      ---    // probably easier/better to just run a fixed number of timesteps
      ---    // using the -i flag in 2.x

      --x    dtfixed() = Real_t(-1.0e-6) ; // Negative means use courant condition
      --x    stoptime()  = Real_t(1.0e-2); // *Real_t(edgeElems*tp/45.0) ;
      this.variables.dtfixed  := -1.0e-6;
      this.variables.stoptime := -1.0e-2;

      ---    // Initial conditions
      --x    deltatimemultlb() = Real_t(1.1) ;
      --x    deltatimemultub() = Real_t(1.2) ;
      --x    dtcourant() = Real_t(1.0e+20) ;
      --x    dthydro()   = Real_t(1.0e+20) ;
      --x    dtmax()     = Real_t(1.0e-2) ;
      --x    time()    = Real_t(0.) ;
      --x    cycle()   = Int_t(0) ;
      this.variables.deltatimemultlb := 1.1;
      this.variables.deltatimemultub := 1.2;
      this.variables.dtcourant       := 1.0e+20;
      this.variables.dthydro         := 1.0e+20;
      this.variables.dtmax           := 1.0e-2;
      this.variables.time            := 0.0;
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
      for elem_index in 0..this.numElem-1 loop
         declare
            local_coords : NodesPerElement_Coordinate_Array;
            elemToNode   : constant NodesPerElement_Index_Array :=
              this.elements(elem_index).node_indexes;
         begin
            for node_index in NodesPerElement_Index_Type loop
               local_coords(node_index) :=
                 this.nodes(elemToNode(node_index)).coordinate;
            end loop;
            declare
               volume     : constant volume_type :=
                 CalcElemVolume (local_coords);
               mass_share : constant Mass_Type :=
                 Mass_Type(volume) / Mass_Type(NODES_PER_ELEMENT);
            begin
               this.elements(elem_index).reference_volume := volume;
               this.elements(elem_index).mass := Mass_Type(volume);
               for node_index in NodesPerElement_Index_Type loop
                  declare
                     node_mass : Mass_Type renames
                       this.nodes(elemToNode(node_index)).mass;
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

---    const Real_t ebase = Real_t(3.948746e+7);
---    Real_t scale = (nx*m_tp)/Real_t(45.0);
---    Real_t einit = ebase*scale*scale*scale;
--    if (m_rowLoc + m_colLoc + m_planeLoc == 0) {
--       // Dump into the first zone (which we know is in the corner)
--       // of the domain that sits at the origin
--       e(0) = einit;
--    }
      declare
         ebase : constant Energy_Type := 3.948746e+7;
         scale : constant Real_Type   := (nx*m_tp)/45.0;
         einit : constant Energy_Type :=
           ebase*Energy_Type(scale*scale*scale);
      begin
         if (this.variables.rowLoc
             + this.variables.colLoc
             + this.variables.planeLoc = 0) then
            this.nodes(0) := einit;
         end if;
--  }
      end;
      --    //set initial deltatime base on analytic CFL calculation
--    deltatime() = (Real_t(.5)*cbrt(volo(0)))/sqrt(Real_t(2.0)*einit);

      return this;
      -- } // End constructor
   end Create;


-- ////////////////////////////////////////////////////////////////////////////////
-- void
-- Domain::BuildMesh(Int_t nx, Int_t edgeNodes, Int_t edgeElems)
-- {
--   Index_t meshEdgeElems = m_tp*nx ;

--   // initialize nodal coordinates
--   Index_t nidx = 0 ;
--   Real_t tz = Real_t(1.125)*Real_t(m_planeLoc*nx)/Real_t(meshEdgeElems) ;
--   for (Index_t plane=0; plane<edgeNodes; ++plane) {
--     Real_t ty = Real_t(1.125)*Real_t(m_rowLoc*nx)/Real_t(meshEdgeElems) ;
--     for (Index_t row=0; row<edgeNodes; ++row) {
--       Real_t tx = Real_t(1.125)*Real_t(m_colLoc*nx)/Real_t(meshEdgeElems) ;
--       for (Index_t col=0; col<edgeNodes; ++col) {
-- 	x(nidx) = tx ;
-- 	y(nidx) = ty ;
-- 	z(nidx) = tz ;
-- 	++nidx ;
-- 	// tx += ds ; // may accumulate roundoff...
-- 	tx = Real_t(1.125)*Real_t(m_colLoc*nx+col+1)/Real_t(meshEdgeElems) ;
--       }
--       // ty += ds ;  // may accumulate roundoff...
--       ty = Real_t(1.125)*Real_t(m_rowLoc*nx+row+1)/Real_t(meshEdgeElems) ;
--     }
--     // tz += ds ;  // may accumulate roundoff...
--     tz = Real_t(1.125)*Real_t(m_planeLoc*nx+plane+1)/Real_t(meshEdgeElems) ;
--   }


--   // embed hexehedral elements in nodal point lattice
--   Index_t zidx = 0 ;
--   nidx = 0 ;
--   for (Index_t plane=0; plane<edgeElems; ++plane) {
--     for (Index_t row=0; row<edgeElems; ++row) {
--       for (Index_t col=0; col<edgeElems; ++col) {
-- 	Index_t *localNode = nodelist(zidx) ;
-- 	localNode[0] = nidx                                       ;
-- 	localNode[1] = nidx                                   + 1 ;
-- 	localNode[2] = nidx                       + edgeNodes + 1 ;
-- 	localNode[3] = nidx                       + edgeNodes     ;
-- 	localNode[4] = nidx + edgeNodes*edgeNodes                 ;
-- 	localNode[5] = nidx + edgeNodes*edgeNodes             + 1 ;
-- 	localNode[6] = nidx + edgeNodes*edgeNodes + edgeNodes + 1 ;
-- 	localNode[7] = nidx + edgeNodes*edgeNodes + edgeNodes     ;
-- 	++zidx ;
-- 	++nidx ;
--       }
--       ++nidx ;
--     }
--     nidx += edgeNodes ;
--   }
-- }


-- ////////////////////////////////////////////////////////////////////////////////
-- void
-- Domain::SetupThreadSupportStructures()
-- {
-- #if _OPENMP
--    Index_t numthreads = omp_get_max_threads();
-- #else
--    Index_t numthreads = 1;
-- #endif

--   if (numthreads > 1) {
--     // set up node-centered indexing of elements
--     Index_t *nodeElemCount = new Index_t[numNode()] ;

--     for (Index_t i=0; i<numNode(); ++i) {
--       nodeElemCount[i] = 0 ;
--     }

--     for (Index_t i=0; i<numElem(); ++i) {
--       Index_t *nl = nodelist(i) ;
--       for (Index_t j=0; j < 8; ++j) {
-- 	++(nodeElemCount[nl[j]] );
--       }
--     }

--     m_nodeElemStart = new Index_t[numNode()+1] ;

--     m_nodeElemStart[0] = 0;

--     for (Index_t i=1; i <= numNode(); ++i) {
--       m_nodeElemStart[i] =
-- 	m_nodeElemStart[i-1] + nodeElemCount[i-1] ;
--     }
--
--     m_nodeElemCornerList = new Index_t[m_nodeElemStart[numNode()]];

--     for (Index_t i=0; i < numNode(); ++i) {
--       nodeElemCount[i] = 0;
--     }

--     for (Index_t i=0; i < numElem(); ++i) {
--       Index_t *nl = nodelist(i) ;
--       for (Index_t j=0; j < 8; ++j) {
-- 	Index_t m = nl[j];
-- 	Index_t k = i*8 + j ;
-- 	Index_t offset = m_nodeElemStart[m] + nodeElemCount[m] ;
-- 	m_nodeElemCornerList[offset] = k;
-- 	++(nodeElemCount[m]) ;
--       }
--     }

--     Index_t clSize = m_nodeElemStart[numNode()] ;
--     for (Index_t i=0; i < clSize; ++i) {
--       Index_t clv = m_nodeElemCornerList[i] ;
--       if ((clv < 0) || (clv > numElem()*8)) {
-- 	fprintf(stderr,
-- 		"AllocateNodeElemIndexes(): nodeElemCornerList entry out of range!\n");
-- #if USE_MPI
-- 	MPI_Abort(MPI_COMM_WORLD, -1);
-- #else
-- 	exit(-1);
-- #endif
--       }
--     }

--     delete [] nodeElemCount ;
--   }
--   else {
--     // These arrays are not used if we're not threaded
--     m_nodeElemStart = NULL;
--     m_nodeElemCornerList = NULL;
--   }
-- }


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
-- #if USE_MPI
--    Index_t myRank;
--    MPI_Comm_rank(MPI_COMM_WORLD, &myRank) ;
--    srand(myRank);
-- #else
--    srand(0);
--    Index_t myRank = 0;
-- #endif
--    this->numReg() = nr;
--    m_regElemSize = new Index_t[numReg()];
--    m_regElemlist = new Index_t*[numReg()];
--    Index_t nextIndex = 0;
--    //if we only have one region just fill it
--    // Fill out the regNumList with material numbers, which are always
--    // the region index plus one
--    if(numReg() == 1) {
--       while (nextIndex < numElem()) {
-- 	 this->regNumList(nextIndex) = 1;
--          nextIndex++;
--       }
--       regElemSize(0) = 0;
--    }
--    //If we have more than one region distribute the elements.
--    else {
--       Int_t regionNum;
--       Int_t regionVar;
--       Int_t lastReg = -1;
--       Int_t binSize;
--       Index_t elements;
--       Index_t runto = 0;
--       Int_t costDenominator = 0;
--       Int_t* regBinEnd = new Int_t[numReg()];
--       //Determine the relative weights of all the regions.  This is based off the -b flag.  Balance is the value passed into b.
--       for (Index_t i=0 ; i<numReg() ; ++i) {
--          regElemSize(i) = 0;
-- 	 costDenominator += pow((i+1), balance);  //Total sum of all regions weights
-- 	 regBinEnd[i] = costDenominator;  //Chance of hitting a given region is (regBinEnd[i] - regBinEdn[i-1])/costDenominator
--       }
--       //Until all elements are assigned
--       while (nextIndex < numElem()) {
-- 	 //pick the region
-- 	 regionVar = rand() % costDenominator;
-- 	 Index_t i = 0;
--          while(regionVar >= regBinEnd[i])
-- 	    i++;
--          //rotate the regions based on MPI rank.  Rotation is Rank % NumRegions this makes each domain have a different region with
--          //the highest representation
-- 	 regionNum = ((i + myRank) % numReg()) + 1;
-- 	 // make sure we don't pick the same region twice in a row
--          while(regionNum == lastReg) {
-- 	    regionVar = rand() % costDenominator;
-- 	    i = 0;
--             while(regionVar >= regBinEnd[i])
-- 	       i++;
-- 	    regionNum = ((i + myRank) % numReg()) + 1;
--          }
-- 	 //Pick the bin size of the region and determine the number of elements.
--          binSize = rand() % 1000;
-- 	 if(binSize < 773) {
-- 	   elements = rand() % 15 + 1;
-- 	 }
-- 	 else if(binSize < 937) {
-- 	   elements = rand() % 16 + 16;
-- 	 }
-- 	 else if(binSize < 970) {
-- 	   elements = rand() % 32 + 32;
-- 	 }
-- 	 else if(binSize < 974) {
-- 	   elements = rand() % 64 + 64;
-- 	 }
-- 	 else if(binSize < 978) {
-- 	   elements = rand() % 128 + 128;
-- 	 }
-- 	 else if(binSize < 981) {
-- 	   elements = rand() % 256 + 256;
-- 	 }
-- 	 else
-- 	    elements = rand() % 1537 + 512;
-- 	 runto = elements + nextIndex;
-- 	 //Store the elements.  If we hit the end before we run out of elements then just stop.
--          while (nextIndex < runto && nextIndex < numElem()) {
-- 	    this->regNumList(nextIndex) = regionNum;
-- 	    nextIndex++;
-- 	 }
-- 	 lastReg = regionNum;
--       }
--    }
--    // Convert regNumList to region index sets
--    // First, count size of each region
--    for (Index_t i=0 ; i<numElem() ; ++i) {
--       int r = this->regNumList(i)-1; // region index == regnum-1
--       regElemSize(r)++;
--    }
--    // Second, allocate each region index set
--    for (Index_t i=0 ; i<numReg() ; ++i) {
--       m_regElemlist[i] = new Index_t[regElemSize(i)];
--       regElemSize(i) = 0;
--    }
--    // Third, fill index sets
--    for (Index_t i=0 ; i<numElem() ; ++i) {
--       Index_t r = regNumList(i)-1;       // region index == regnum-1
--       Index_t regndx = regElemSize(r)++; // Note increment
--       regElemlist(r,regndx) = i;
--    }
--
-- }

-- /////////////////////////////////////////////////////////////
-- void
-- Domain::SetupSymmetryPlanes(Int_t edgeNodes)
-- {
--   Index_t nidx = 0 ;
--   for (Index_t i=0; i<edgeNodes; ++i) {
--     Index_t planeInc = i*edgeNodes*edgeNodes ;
--     Index_t rowInc   = i*edgeNodes ;
--     for (Index_t j=0; j<edgeNodes; ++j) {
--       if (m_planeLoc == 0) {
-- 	m_symmZ[nidx] = rowInc   + j ;
--       }
--       if (m_rowLoc == 0) {
-- 	m_symmY[nidx] = planeInc + j ;
--       }
--       if (m_colLoc == 0) {
-- 	m_symmX[nidx] = planeInc + j*edgeNodes ;
--       }
--       ++nidx ;
--     }
--   }
-- }



-- /////////////////////////////////////////////////////////////
-- void
-- Domain::SetupElementConnectivities(Int_t edgeElems)
-- {
--    lxim(0) = 0 ;
--    for (Index_t i=1; i<numElem(); ++i) {
--       lxim(i)   = i-1 ;
--       lxip(i-1) = i ;
--    }
--    lxip(numElem()-1) = numElem()-1 ;

--    for (Index_t i=0; i<edgeElems; ++i) {
--       letam(i) = i ;
--       letap(numElem()-edgeElems+i) = numElem()-edgeElems+i ;
--    }
--    for (Index_t i=edgeElems; i<numElem(); ++i) {
--       letam(i) = i-edgeElems ;
--       letap(i-edgeElems) = i ;
--    }

--    for (Index_t i=0; i<edgeElems*edgeElems; ++i) {
--       lzetam(i) = i ;
--       lzetap(numElem()-edgeElems*edgeElems+i) = numElem()-edgeElems*edgeElems+i ;
--    }
--    for (Index_t i=edgeElems*edgeElems; i<numElem(); ++i) {
--       lzetam(i) = i - edgeElems*edgeElems ;
--       lzetap(i-edgeElems*edgeElems) = i ;
--    }
-- }

-- /////////////////////////////////////////////////////////////
-- void
-- Domain::SetupBoundaryConditions(Int_t edgeElems)
-- {
--   Index_t ghostIdx[6] ;  // offsets to ghost locations

--   // set up boundary condition information
--   for (Index_t i=0; i<numElem(); ++i) {
--      elemBC(i) = Int_t(0) ;
--   }

--   for (Index_t i=0; i<6; ++i) {
--     ghostIdx[i] = INT_MIN ;
--   }

--   Int_t pidx = numElem() ;
--   if (m_planeMin != 0) {
--     ghostIdx[0] = pidx ;
--     pidx += sizeX()*sizeY() ;
--   }

--   if (m_planeMax != 0) {
--     ghostIdx[1] = pidx ;
--     pidx += sizeX()*sizeY() ;
--   }

--   if (m_rowMin != 0) {
--     ghostIdx[2] = pidx ;
--     pidx += sizeX()*sizeZ() ;
--   }

--   if (m_rowMax != 0) {
--     ghostIdx[3] = pidx ;
--     pidx += sizeX()*sizeZ() ;
--   }

--   if (m_colMin != 0) {
--     ghostIdx[4] = pidx ;
--     pidx += sizeY()*sizeZ() ;
--   }

--   if (m_colMax != 0) {
--     ghostIdx[5] = pidx ;
--   }

--   // symmetry plane or free surface BCs
--   for (Index_t i=0; i<edgeElems; ++i) {
--     Index_t planeInc = i*edgeElems*edgeElems ;
--     Index_t rowInc   = i*edgeElems ;
--     for (Index_t j=0; j<edgeElems; ++j) {
--       if (m_planeLoc == 0) {
-- 	elemBC(rowInc+j) |= ZETA_M_SYMM ;
--       }
--       else {
-- 	elemBC(rowInc+j) |= ZETA_M_COMM ;
-- 	lzetam(rowInc+j) = ghostIdx[0] + rowInc + j ;
--       }

--       if (m_planeLoc == m_tp-1) {
-- 	elemBC(rowInc+j+numElem()-edgeElems*edgeElems) |=
-- 	  ZETA_P_FREE;
--       }
--       else {
-- 	elemBC(rowInc+j+numElem()-edgeElems*edgeElems) |=
-- 	  ZETA_P_COMM ;
-- 	lzetap(rowInc+j+numElem()-edgeElems*edgeElems) =
-- 	  ghostIdx[1] + rowInc + j ;
--       }

--       if (m_rowLoc == 0) {
-- 	elemBC(planeInc+j) |= ETA_M_SYMM ;
--       }
--       else {
-- 	elemBC(planeInc+j) |= ETA_M_COMM ;
-- 	letam(planeInc+j) = ghostIdx[2] + rowInc + j ;
--       }

--       if (m_rowLoc == m_tp-1) {
-- 	elemBC(planeInc+j+edgeElems*edgeElems-edgeElems) |=
-- 	  ETA_P_FREE ;
--       }
--       else {
-- 	elemBC(planeInc+j+edgeElems*edgeElems-edgeElems) |=
-- 	  ETA_P_COMM ;
-- 	letap(planeInc+j+edgeElems*edgeElems-edgeElems) =
-- 	  ghostIdx[3] +  rowInc + j ;
--       }

--       if (m_colLoc == 0) {
-- 	elemBC(planeInc+j*edgeElems) |= XI_M_SYMM ;
--       }
--       else {
-- 	elemBC(planeInc+j*edgeElems) |= XI_M_COMM ;
-- 	lxim(planeInc+j*edgeElems) = ghostIdx[4] + rowInc + j ;
--       }

--       if (m_colLoc == m_tp-1) {
-- 	elemBC(planeInc+j*edgeElems+edgeElems-1) |= XI_P_FREE ;
--       }
--       else {
-- 	elemBC(planeInc+j*edgeElems+edgeElems-1) |= XI_P_COMM ;
-- 	lxip(planeInc+j*edgeElems+edgeElems-1) =
-- 	  ghostIdx[5] + rowInc + j ;
--       }
--     }
--   }
-- }

--- ///////////////////////////////////////////////////////////////////////////
--x void InitMeshDecomp(Int_t numRanks, Int_t myRank,
--x                     Int_t *col, Int_t *row, Int_t *plane, Int_t *side)
--x {
   procedure InitMeshDecomp
     (numRanks : in Int_t;
      myRank   : in Int_t;
      col      : out Int_t;
      row      : out Int_t;
      plane    : out Int_t;
      side     : out Int_t)
   is
      --x    Int_t testProcs;
      --x    Int_t dx, dy, dz;
      --x    Int_t myDom;
      testProcs : Int_t;
      dx        : Int_t;
      dy        : Int_t;
      dz        : Int_t;
      myDom     : Int_t;
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
      testProcs := Int_t(cbrt(Real_t(numRanks))+0.5) ;
      if (testProcs*testProcs*testProcs /= numRanks) then
         raise Usage_Error
           with "Num processors (" & numRanks'Img &
           ") must be a cube of an integer (1, 8, 27, ...)"
      end if;
      --x    if (sizeof(Real_t) != 4 && sizeof(Real_t) != 8) {
      --x       printf("MPI operations only support float and double right now...\n");
      -- #if USE_MPI
      --       MPI_Abort(MPI_COMM_WORLD, -1) ;
      -- #else
      --x       exit(-1);
      -- #endif
      --x    }
      if (sizeof(Real_t) /= 4 and sizeof(Real_t) /= 8) then
          raise Usage_Error
           with "MPI operations only support float and double right now." &
           "  Size of Real_t is " & sizeof(Real_t);
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
      dx := testProcs;
      dy := testProcs;
      dz := testProcs;

--    // temporary test
--    if (dx*dy*dz != numRanks) {
--       printf("error -- must have as many domains as procs\n") ;
-- #if USE_MPI
--       MPI_Abort(MPI_COMM_WORLD, -1) ;
-- #else
--       exit(-1);
-- #endif
--    }
--    Int_t remainder = dx*dy*dz % numRanks ;
--    if (myRank < remainder) {
--       myDom = myRank*( 1+ (dx*dy*dz / numRanks)) ;
--    }
--    else {
--       myDom = remainder*( 1+ (dx*dy*dz / numRanks)) +
--          (myRank - remainder)*(dx*dy*dz/numRanks) ;
--    }

--    *col = myDom % dx ;
--    *row = (myDom / dx) % dy ;
--    *plane = myDom / (dx*dy) ;
--    *side = testProcs;

--    return;
-- }
end InitMeshDecomp;

end LULESH.Init;
