with Ada.Calendar;
with Ada.Containers.Vectors;
with Ada.Numerics.Generic_Elementary_Functions;
with Ada.Numerics.Generic_Real_Arrays;
with Ada.Real_Time;

with Interfaces;

package LULESH is

   package AC renames Ada.Calendar;
   package ART renames Ada.Real_Time;

   Usage_Error  : Exception;
   Coding_Error : Exception;

   --- #if !defined(USE_MPI)
   --- # error "You should specify USE_MPI=0 or USE_MPI=1 on the compile line"
   --- #endif
   USE_MPI : constant Boolean := False;


   --- // OpenMP will be compiled in if this flag is set to 1 AND the compiler beging
   --- // used supports it (i.e. the _OPENMP symbol is defined)
   --x #define USE_OMP 1
   USE_OMP : constant Boolean := False;

   -- #if USE_MPI
   -- #include <mpi.h>

   -- /*
   --    define one of these three symbols:

   --    SEDOV_SYNC_POS_VEL_NONE
   --    SEDOV_SYNC_POS_VEL_EARLY
   --    SEDOV_SYNC_POS_VEL_LATE
   -- */

   -- #define SEDOV_SYNC_POS_VEL_EARLY 1
   -- #endif

   --- #include <math.h>
   --- #include <vector>

   --- //**************************************************
   --- // Allow flexibility for arithmetic representations
   --- //**************************************************

   -- #define MAX(a, b) ( ((a) > (b)) ? (a) : (b))


   --- // Precision specification
   --x typedef float        real4 ;
   --x typedef double       real8 ;
   --x typedef long double  real10 ;  // 10 bytes on x86
   --x typedef int    Index_t ; // array subscript and loop index
   --x typedef real8  Real_t ;  // floating point representation
   --x typedef int    Int_t ;   // integer representation
   type real4        is new Interfaces.IEEE_Float_32;
   type real8        is new Interfaces.IEEE_Float_64;
   type real10       is new Interfaces.IEEE_Extended_Float;
   type Index_Type   is new Natural;
   type Real_Type    is new real8;
   type Int_t        is new Integer;
   type Balance_Type is new Natural;
   type Cost_Type    is new Natural;
   ---------------------------------------

   type Element_Index is new Index_Type;
   type Region_Index  is new Index_Type;
   type Bisect_Range  is new Index_Type range 0..1;
   type Domain_Index  is new Index_Type;

   type Rank_Type              is new Natural;
   subtype Process_ID_Type     is Rank_Type;
   subtype Rank_Count_Range    is Rank_Type range 1..Rank_Type'Last;
   subtype Process_Count_Range is Rank_Count_Range;

   --- Looking down on element:
   --- 1-2-3-4 are nodes going CCW on bottom.
   --- 5-6-7-8 are nodes directly above each.
   NODES_PER_ELEMENT : constant := 8;
   NODES_PER_FACE    : constant := 4;
   type Node_Index               is new Index_Type;
   subtype NodesPerElement_Range is Node_Index range 0..NODES_PER_ELEMENT-1;
   subtype NodesPerFace_Range    is Node_Index range 0..NODES_PER_FACE-1;

   FACES_PER_ELEMENT : constant := 6;
   type Face_Range is new Index_Type range 0..FACES_PER_ELEMENT-1;

   type Real_Array is array (Index_Type range <>) of Real_Type;
   --     type Real_Array_Access is access Real_Array;
   type Index_Array is array (Index_Type range <>) of Index_Type;
   type Element_Index_Array is array (Element_Index range <>)
     of Element_Index;
   type Element_Index_Array_Access is access Element_Index_Array;

   type Node_Node_Index_Array is array (Node_Index range <>)
     of Node_Index;
   type Node_Node_Index_Array_Access is access Node_Node_Index_Array;
   type Node_Element_Index_Array is array (Node_Index range <>)
     of Element_Index;
   type Node_Element_Index_Array_Access is access Node_Element_Index_Array;

   type Region_Bin_End_Array is array (Region_Index range <>)
     of Cost_Type;
   type Region_Bin_End_Array_Access is access Region_Bin_End_Array;

   --- Cartesian coordinates:
   --- X is positive in direction node 4 -> node 1
   --- Y is positive in direction node 4 -> node 3
   --- X is positive in direction node 4 -> node 8
   subtype Cartesian_Axes is integer range 1..3;
   X : constant Cartesian_Axes := 1;
   Y : constant Cartesian_Axes := 2;
   Z : constant Cartesian_Axes := 3;
   -- type Cartesian_Axes is (X, Y, Z);

   --- Natural Coordinates (isoparametric representation):
   --- Xi   is positive in direction +Y: from face 1584 to 2673.
   --- Eta  is positive in direction -X: from face 1562 to 4873.
   --- Zeta is positive in direction +Z: from face 1234 to 5678.
   type Natural_Axes is (Xi, Eta, Zeta);
   Et : constant Natural_Axes := Eta;
   Ze : constant Natural_Axes := Zeta;
   subtype Natural_Magnitude_Range is Real_Type range -1.0 .. 1.0;

   type FacesPerNode_Element_Index_Array is
     array (Face_Range) of Element_Index;
   type NodesPerElement_Index_Array is
     array (NodesPerElement_Range) of Node_Index;
   type NodesPerElement_Index_Array_Array is array
     (Element_Index range <>) of NodesPerElement_Index_Array;
   type Face_NodesPerFace_NodesPerElement_Array is array
     (Face_Range, NodesPerFace_Range) of NodesPerElement_Range;

   type Domain_Record is private;
   type Domain_Access is access Domain_Record;

   ---------------------------------------

   --x enum { VolumeError = -1, QStopError = -2 } ;
   VolumeError : exception;
   QStop_Error : exception;

   --x inline real4  SQRT(real4  arg) { return sqrtf(arg) ; }
   --x inline real8  SQRT(real8  arg) { return sqrt(arg) ; }
   --x inline real10 SQRT(real10 arg) { return sqrtl(arg) ; }
   package real4_Elementary_functions is new
     Ada.Numerics.Generic_Elementary_Functions (real4);
   package real8_Elementary_functions is new
     Ada.Numerics.Generic_Elementary_Functions (real8);
   package real10_Elementary_functions is new
     Ada.Numerics.Generic_Elementary_Functions (real10);

   function SQRT (This : in real4) return real4 renames
     real4_Elementary_functions.Sqrt;
   function SQRT (This : in real8) return real8 renames
     real8_Elementary_functions.Sqrt;
   function SQRT (This : in real10) return real10 renames
     real10_Elementary_functions.Sqrt;

   --x inline real4  CBRT(real4  arg) { return cbrtf(arg) ; }
   --x inline real8  CBRT(real8  arg) { return cbrt(arg) ; }
   --x inline real10 CBRT(real10 arg) { return cbrtl(arg) ; }
   function CBRT (This : in real4) return real4 is
     (real4_Elementary_functions."**"(This, 1.0/3.0));
   function CBRT (This : in real8) return real8 is
     (real8_Elementary_functions."**"(This, 1.0/3.0));
   function CBRT (This : in real10) return real10 is
     (real10_Elementary_functions."**"(This, 1.0/3.0));
   --x inline real4  FABS(real4  arg) { return fabsf(arg) ; }
   --x inline real8  FABS(real8  arg) { return fabs(arg) ; }
   --x inline real10 FABS(real10 arg) { return fabsl(arg) ; }

   --- // Stuff needed for boundary conditions
   --- // 2 BCs on each of 6 hexahedral faces (12 bits)
   --x #define XI_M        0x00007
   --x #define XI_M_SYMM   0x00001
   --x #define XI_M_FREE   0x00002
   --x #define XI_M_COMM   0x00004
   --x #define XI_P        0x00038
   --x #define XI_P_SYMM   0x00008
   --x #define XI_P_FREE   0x00010
   --x #define XI_P_COMM   0x00020
   --x #define ETA_M       0x001c0
   --x #define ETA_M_SYMM  0x00040
   --x #define ETA_M_FREE  0x00080
   --x #define ETA_M_COMM  0x00100
   --x #define ETA_P       0x00e00
   --x #define ETA_P_SYMM  0x00200
   --x #define ETA_P_FREE  0x00400
   --x #define ETA_P_COMM  0x00800
   --x #define ZETA_M      0x07000
   --x #define ZETA_M_SYMM 0x01000
   --x #define ZETA_M_FREE 0x02000
   --x #define ZETA_M_COMM 0x04000
   --x #define ZETA_P      0x38000
   --x #define ZETA_P_SYMM 0x08000
   --x #define ZETA_P_FREE 0x10000
   --x #define ZETA_P_COMM 0x20000

   --- // MPI Message Tags
   --x #define MSG_COMM_SBN      1024
   --x #define MSG_SYNC_POS_VEL  2048
   --x #define MSG_MONOQ         3072
   MSG_COMM_SBN     : constant := 2#0010_0000_0000#;
   MSG_SYNC_POS_VEL : constant := 2#0100_0000_0000#;
   MSG_MONOQ        : constant := 2#0110_0000_0000#;

   --x #define MAX_FIELDS_PER_MPI_COMM 6
   MAX_FIELDS_PER_MPI_COMM : constant := 6;

   --- // Assume 128 byte coherence
   --- // Assume Real_t is an "integral power of 2" bytes wide
   --x #define CACHE_COHERENCE_PAD_REAL (128 / sizeof(Real_t))
   Byte_Size: constant := 8;
   CACHE_COHERENCE_PAD_REAL : constant := 128 / (Real_Type'Size / Byte_Size);

   -- #define CACHE_ALIGN_REAL(n) \
   --    (((n) + (CACHE_COHERENCE_PAD_REAL - 1)) & ~(CACHE_COHERENCE_PAD_REAL-1))

   --- //////////////////////////////////////////////////////
   --- // Primary data structure
   --- //////////////////////////////////////////////////////

   --- /*
   ---  * The implementation of the data abstraction used for lulesh
   ---  * resides entirely in the Domain class below.  You can change
   ---  * grouping and interleaving of fields here to maximize data layout
   ---  * efficiency for your underlying architecture or compiler.
   ---  *
   ---  * For example, fields can be implemented as STL objects or
   ---  * raw array pointers.  As another example, individual fields
   ---  * m_x, m_y, m_z could be budled into
   ---  *
   ---  *    struct { Real_t x, y, z ; } *m_coord ;
   ---  *
   ---  * allowing accessor functions such as
   ---  *
   ---  *  "Real_t &x(Index_t idx) { return m_coord[idx].x ; }"
   ---  *  "Real_t &y(Index_t idx) { return m_coord[idx].y ; }"
   ---  *  "Real_t &z(Index_t idx) { return m_coord[idx].z ; }"
   ---  */


   -- class Domain {

   --    // Nodes on symmertry planes
   --    bool symmXempty()          { return m_symmX.empty(); }
   --    bool symmYempty()          { return m_symmY.empty(); }
   --    bool symmZempty()          { return m_symmZ.empty(); }

   ---    //
   ---    // Element-centered
   ---    //
   --x    Index_t*  nodelist(Index_t idx)    { return &m_nodelist[Index_t(8)*idx] ; }

   --x    Index_t nodeElemCount(Index_t idx)
   --x    { return m_nodeElemStart[idx+1] - m_nodeElemStart[idx] ; }
   --x    Index_t *nodeElemCornerList(Index_t idx)
   --x    { return &m_nodeElemCornerList[m_nodeElemStart[idx]] ; }
--     function nodeElemCount 
--       (this : in Domain_Record;
--        idx  : in Index_Type) 
--        return Index_Type is 
--       (this.nodeElemStart(idx+1) - this.nodeElemStart(idx))
--     with inline;
--     function nodeElemCornerList 
--       (this : in Domain_Record;
--        idx  : in Index_Type) 
--        return Corner_List is 
--       (this.nodeElemCornerList(this.nodeElemStart(idx)..this.nodeElemStart(idx+1)))
--     with inline;
   
   --    //
   --    // MPI-Related additional data
   --    //

   -- #if USE_MPI
   --    // Communication Work space
   --    Real_t *commDataSend ;
   --    Real_t *commDataRecv ;

   --    // Maximum number of block neighbors
   --    MPI_Request recvRequest[26] ; // 6 faces + 12 edges + 8 corners
   --    MPI_Request sendRequest[26] ; // 6 faces + 12 edges + 8 corners
   -- #endif


   --   private:
   --- International System of Units:
   type SI_Type is new Real_Type;
   
   type Cubic_Meters                 is new SI_Type; -- = m**3
   type Square_Meters                is new SI_Type; -- = m**2
   type Joules                       is new SI_Type; -- = N*m
   type Kilograms                    is new SI_Type; -- "kg"
   type Kilograms_Per_Cubic_Meter    is new SI_Type; -- kg/m**3
   type Meters                       is new SI_Type; -- "m"
   type Meters_Per_Second            is new SI_Type; -- = m/s
   type Meters_Per_Second_Per_Second is new SI_Type; -- = m/s**2 
   type Newtons                      is new SI_Type; -- "N" = kg*m/s**2
   type Pascals                      is new SI_Type; -- "Pa" = N/m**2
   
   subtype Acceleration is Meters_Per_Second_Per_Second;
   subtype Density      is Kilograms_Per_Cubic_Meter;  
   subtype Energy       is Joules;  
   subtype Force        is Newtons;   
   subtype Length       is Meters;
   subtype Mass    is Kilograms;
   subtype Pressure     is Pascals;
   subtype Velocity     is Meters_Per_Second;
   subtype Time         is ART.Time;
   subtype Time_Span    is ART.Time_Span;
   subtype Volume       is Cubic_Meters;
   
   type Gradient_Type is new Real_Type;

   -- Provides matrix and vector math:
   package Length_Matrices   is new Ada.Numerics.Generic_Real_Arrays (Length);

   type Acceleration_Vector  is array (Cartesian_Axes) of Acceleration;
   type Coordinate_Vector    is new Length_Matrices.Real_Vector (Cartesian_Axes);
   type Force_Vector         is array (Cartesian_Axes) of Force;
   type Gradient_Vector      is array (Natural_Axes)   of Gradient_Type;
   type Strain_Vector        is array (Cartesian_Axes) of Length;
   type Velocity_Vector      is array (Cartesian_Axes) of Velocity;

   subtype Size_Type         is Element_Index;
   type Cartesian_Size_Array is array (Cartesian_Axes) of Size_Type;
   

   type Force_Vector_Array is array (Node_Index range <>) of Force_Vector;
   type Force_Vector_Array_Access is access Force_Vector_Array;
--     type Gradient_Array     is array (Node_Index range <>) of Gradient_Vector;
   type Coordinate_Array   is array (Node_Index range <>) of Coordinate_Vector;
--
   type NodesPerElement_Coordinate_Array is new Coordinate_Array
     (NodesPerElement_Range);
   type NodesPerFace_Coordinate_Array is new Coordinate_Array
     (NodesPerFace_Range);
--     type NodesPerElement_Force_Array is new Force_Array
--       (NodesPerElement_Range);
   type NodesPerElement_Force_Vector_Array is new Force_Vector_Array
     (NodesPerElement_Range);
--     type NodesPerElement_Force_Vector_Array_Array
--       is array (Element_Index) of NodesPerElement_Force_Vector_Array;
--     type NodesPerFace_Coordinate_Array is new Coordinate_Array
--       (NodesPerFace_Range);
--
--     type Force_Vector_Array_Access is access Force_Vector_Array;
--     type NodesPerElement_Force_Vector_Array_Array_Access is access
--       NodesPerElement_Force_Vector_Array_Array;
--
--     package Element_Index_Vectors is
--       new Ada.Containers.Vectors
--         (Index_Type   => Element_Index,
--          Element_Type => Element_Index);
--
--     type Region_Element_Vector_Array is array (Region_Index range <>)
--       of Element_Index_Vectors.Vector;
--
--     type Element_To_Region_Array is array (Element_Index range <>)
--       of Region_Index;
--
   type Symmetry_Node_Array is array (Cartesian_Axes) of Node_Index;

   type Node_Record is record
      --x    std::vector<Real_t> m_x ;  /* coordinates */
      --x    std::vector<Real_t> m_y ;
      --x    std::vector<Real_t> m_z ;
      coordinate          : Coordinate_Vector;
      --x    std::vector<Index_t> m_symmX ;  /* symmetry plane nodesets */
      --x    std::vector<Index_t> m_symmY ;
      --x    std::vector<Index_t> m_symmZ ;
      symmetry_plane_nodes : Symmetry_Node_Array;
      --x    std::vector<Real_t> m_xd ; /* velocities */
      --x    std::vector<Real_t> m_yd ;
      --x    std::vector<Real_t> m_zd ;
      velocity            : Velocity_Vector;
      --x    std::vector<Real_t> m_xdd ; /* accelerations */
      --x    std::vector<Real_t> m_ydd ;
      --x    std::vector<Real_t> m_zdd ;
      acceleration        : Acceleration_Vector;
      --x    std::vector<Real_t> m_fx ;  /* forces */
      --x    std::vector<Real_t> m_fy ;
      --x    std::vector<Real_t> m_fz ;
      force               : Force_Vector;
      --x    std::vector<Real_t> m_nodalMass ;  /* mass */
      nmass                : Mass;
   end record;
   type Node_Array is array (Node_Index range <>) of Node_Record;
   type Node_Array_Access is access Node_Array;

   type MP_Type is (M, P);
   type Connectivity_Array is
     array (Natural_Axes, MP_Type) of Element_Index;

   type Boundary_Condition_Type is (Symm, Free, Common);
   type Boundary_Condition_Array is array
     (Natural_Axes, MP_Type, Boundary_Condition_Type) of Boolean
     with pack;

   type Element_Record is record
      --x    std::vector<Index_t>  m_nodelist ;     /* elemToNode connectivity */
      node_indexes : NodesPerElement_Index_Array;
      --x    std::vector<Index_t>  m_lxim ;  /* element connectivity across each face */
      --x    std::vector<Index_t>  m_lxip ;
      --x    std::vector<Index_t>  m_letam ;
      --x    std::vector<Index_t>  m_letap ;
      --x    std::vector<Index_t>  m_lzetam ;
      --x    std::vector<Index_t>  m_lzetap ;
      connections : Connectivity_Array;
      ---    // elem face symm/free-surface flag
      --x    std::vector<Int_t>    m_elemBC ;  /* symmetry/free-surface flags for each elem face */
      elemBC : Boundary_Condition_Array;
      --x    std::vector<Real_t> m_dxx ;  /* principal strains -- temporary */
      --x    std::vector<Real_t> m_dyy ;
      --x    std::vector<Real_t> m_dzz ;
      principal_strain : Strain_Vector;
      --x    std::vector<Real_t> m_delv_xi ;    /* velocity gradient -- temporary */
      --x    std::vector<Real_t> m_delv_eta ;
      --x    std::vector<Real_t> m_delv_zeta ;
      velocity_gradient : Gradient_Vector;
      ---    // Position gradient - temporary
      --x    std::vector<Real_t> m_delx_xi ;    /* coordinate gradient -- temporary */
      --x    std::vector<Real_t> m_delx_eta ;
      --x    std::vector<Real_t> m_delx_zeta ;
      position_gradient : Gradient_Vector;
      --x    std::vector<Real_t> m_e ;   /* energy */
      eenergy : Energy;
      --x    std::vector<Real_t> m_p ;   /* pressure */
      --x    std::vector<Real_t> m_q ;   /* q */
      --x    std::vector<Real_t> m_ql ;  /* linear term for q */
      --x    std::vector<Real_t> m_qq ;  /* quadratic term for q */
      static_pressure            : Pressure;
      dynamic_pressure           : Pressure;
      dynamic_pressure_linear    : Pressure;
      dynamic_pressure_quadratic : Pressure;
      sig                        : Force_Vector;
      --x    std::vector<Real_t> m_v ;     /* relative volume */
      --x    std::vector<Real_t> m_volo ;  /* reference volume */
      --x    std::vector<Real_t> m_vnew ;  /* new relative volume -- temporary */
      --x    std::vector<Real_t> m_delv ;  /* m_vnew - m_v */
      --x    std::vector<Real_t> m_vdov ;  /* volume derivative over volume */
      relative_volume               : Volume;
      reference_volume              : Volume;
      new_relative_volume           : Volume;
      relative_volume_delta         : Volume;
      volume_derivative_over_volume : Volume;
      --x    std::vector<Real_t> m_arealg ;  /* characteristic length of an element */
      characteristic_length : Length;
      --x    std::vector<Real_t> m_ss ;      /* "sound speed" */
      sound_speed : Velocity;
      --x    std::vector<Real_t> m_elemMass ;  /* mass */
      --x    // Element mass
      emass : Mass;
      --x    Index_t *m_regNumList ;    // Region number per domain element
      region_number : Region_Index;
   end record;
   type Element_Array is array (Element_Index range <>) of Element_Record;
   type Element_Array_Access is access Element_Array;

   type Region_Record is record
      --x    Index_t *m_regElemSize ;   // Size of region sets
      size     : Element_Index;
      --x    Index_t **m_regElemlist ;  // region indexset
      elements : Element_Index_Array_Access;
   end record;
   type Region_Array is array (Region_Index range <>) of Region_Record;
   type Region_Array_Access is access Region_Array;

   type Variables_Record is record
      ---    // Variables to keep track of timestep, simulation time, and cycle
      --x    Real_t  m_dtcourant ;         // courant constraint
      --x    Real_t  m_dthydro ;           // volume change constraint
      --x    Int_t   m_cycle ;             // iteration count for simulation
      --x    Real_t  m_dtfixed ;           // fixed time increment
      --x    Real_t  m_time ;              // current time
      --x    Real_t  m_deltatime ;         // variable time increment
      --x    Real_t  m_deltatimemultlb ;
      --x    Real_t  m_deltatimemultub ;
      --x    Real_t  m_dtmax ;             // maximum allowable time increment
      --x    Real_t  m_stoptime ;          // end time for simulation
      dtcourant                         : Time_Span;
      dthydro                           : Time_Span;
      cycle                             : Int_t;
      dtfixed                           : Time_Span;
      use_courant_condition             : Boolean;
      current_time                      : Time;
      deltatime                         : Time_Span;
      delta_time_multiplier_lower_bound : Real_Type;
      delta_time_multiplier_upper_bound : Real_Type;
      dtmax                             : Time_Span;
      stoptime                          : Time;

      --x    Int_t   m_numRanks ;
      numRanks : Rank_Count_Range;

      --x    Index_t m_colLoc ;
      --x    Index_t m_rowLoc ;
      --x    Index_t m_planeLoc ;
      --x    Index_t m_tp ;
      colLoc   : Domain_Index;
      rowLoc   : Domain_Index;
      planeLoc : Domain_Index;
      tp       : Domain_Index;

      ---    // These arrays are not used if we're not threaded
      ---    // OMP hack
      --    Index_t *m_nodeElemStart ;
      --    Index_t *m_nodeElemCornerList ;
      nodeElemStart      : Node_Element_Index_Array_Access;
      nodeElemCornerList : Node_Element_Index_Array_Access;

      --x    Index_t m_maxPlaneSize ;
      --x    Index_t m_maxEdgeSize ;
      maxPlaneSize : Element_Index;
      maxEdgeSize  : Element_Index;

      ---    // Used in setup
      --x    Index_t m_rowMin, m_rowMax;
      --x    Index_t m_colMin, m_colMax;
      --x    Index_t m_planeMin, m_planeMax ;
      rowMin   : Element_Index;
      rowMax   : Element_Index;
      colMin   : Element_Index;
      colMax   : Element_Index;
      planeMin : Element_Index;
      planeMax : Element_Index;
   end record;

   type Parameters_Record is record
      ---    // Parameters

      ---    // Cutoffs (treat as constants)

      ---    const Real_t  m_e_cut ;             // energy tolerance
      ---    const Real_t  m_p_cut ;             // pressure tolerance
      ---    const Real_t  m_q_cut ;             // q tolerance
      ---    const Real_t  m_v_cut ;             // relative volume tolerance
      ---    const Real_t  m_u_cut ;             // velocity tolerance
      energy_tolerance           : Energy;
      pressure_tolerance         : Pressure;
      dynamic_pressure_tolerance : Pressure;
      relative_volume_tolerance  : Volume;
      velocity_tolerance         : Velocity;

      ---    // Other constants (usually setable, but hardcoded in this proxy app)

      ---    const Real_t  m_hgcoef ;            // hourglass control
      ---    const Real_t  m_ss4o3 ;
      ---    const Real_t  m_qstop ;             // excessive q indicator
      ---    const Real_t  m_monoq_max_slope ;
      ---    const Real_t  m_monoq_limiter_mult ;
      ---    const Real_t  m_qlc_monoq ;         // linear term coef for q
      ---    const Real_t  m_qqc_monoq ;         // quadratic term coef for q
      ---    const Real_t  m_qqc ;
      ---    const Real_t  m_eosvmax ;
      ---    const Real_t  m_eosvmin ;
      ---    const Real_t  m_pmin ;              // pressure floor
      ---    const Real_t  m_emin ;              // energy floor
      ---    const Real_t  m_dvovmax ;           // maximum allowable volume change
      ---    const Real_t  m_refdens ;           // reference density
      hgcoef             : Real_Type;
      four_thirds        : Real_Type;
      qstop              : Pressure;
      monoq_max_slope    : Real_Type;
      monoq_limiter_mult : Real_Type;
      qlc_monoq          : Real_Type;
      qqc_monoq          : Real_Type;
      qqc                : Real_Type;
      eosvmax            : Real_Type;
      eosvmin            : Real_Type;
      pressure_floor     : Pressure;
      energy_floor       : Energy;
      volume_delta_max   : Volume;
      reference_density  : Density;
      --x    Int_t    m_cost; //imbalance cost
      imbalance_cost     : Cost_Type;

      --x    Index_t m_sizeX ;
      --x    Index_t m_sizeY ;
      --x    Index_t m_sizeZ ;
      size : Cartesian_Size_Array;
   end record;

   -- typedef Real_t &(Domain::* Domain_member )(Index_t) ;
   -- (Address of (pointer to) a Domain class function (accessor):
   --   that takes an Index_t and returns the Real_t Domain component value)
--     type Accessor_access is access function (Index : in Index_Type) return Real_Type;
--     type Domain_member is access all Real_Type;

private

   function "*" (Left : Time_Span; Right : Real_Type) return Time_Span is
     (ART.To_Time_Span (ART.To_Duration(Left) * Duration (Right)))
   with Inline;
   
   function "*" (Left  : Real_Type; Right : Time_Span) return Time_Span is
     (Right * Left)
   with Inline;
 
   function "/" (Left  : Time_Span; Right : Real_Type) return Time_Span is
     (ART.To_Time_Span (ART.To_Duration(Left) / Duration (Right)))
   with Inline;

   type Domain_Record is record
      --x    Index_t m_numElem ;
      --x    Index_t m_numNode ;
      numNode    : Node_Index;
      numElem    : Element_Index;
      numReg     : Region_Index;
      nodes      : Node_Array_Access;
      regions    : Region_Array_Access;
      elements   : Element_Array_Access;
      variables  : Variables_Record;
      parameters : Parameters_Record;
   end record;

   --- /******************************************/

   --- /* Work Routines */

   --x static inline
   --x void TimeIncrement(Domain& domain)
   procedure TimeIncrement (domain : in out Domain_Record)
   with Inline;

   --x static inline
   --x void LagrangeLeapFrog(Domain& domain)
   procedure LagrangeLeapFrog (domain : in out Domain_Record)
   with Inline;

end LULESH;
