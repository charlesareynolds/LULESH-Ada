with Ada.Calendar;
with Ada.Containers.Vectors;
with Ada.Numerics.Generic_Elementary_Functions;
with Interfaces;

package LULESH is

   package AC renames Ada.Calendar;

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
   type real4      is new Interfaces.IEEE_Float_32;
   type real8      is new Interfaces.IEEE_Float_64;
   type real10     is new Interfaces.IEEE_Extended_Float;
   type Index_Type is new Natural;
   type Real_Type  is new real8;
   type Int_t      is new Integer;

   ---------------------------------------

   subtype Bisect_Range               is Index_Type range 0..1;
   subtype Element_Index_Type         is Index_Type;
   subtype Face_Index_Type            is Index_Type range 0..5;
   subtype NodesPerElement_Index_Type is Index_Type range 0..7;
   subtype NodesPerFace_Index_Type    is Index_Type range 0..3;
   subtype Node_Index_Type            is Index_Type;
   subtype Region_Index_Type          is Index_Type;

   type Real_Array is array (Index_Type range <>) of Real_Type;
   --     type Real_Array_Access is access Real_Array;
   type Index_Array is array (Index_Type range <>) of Index_Type;

   type XYZ_Names           is (X, Y, Z);
   subtype Dimensions       is XYZ_Names;
   subtype Coordinate_Names is XYZ_Names;

   type XEZ_Type is (Xi, Eta, Zeta);
   Et : constant XEZ_Type := Eta;
   Ze : constant XEZ_Type := Zeta;

--     type XYZ_Real_Array is array
--       (XYZ_Names) of Real_Type;

--     type Bisect_XYZ_Real_Array is array
--       (Bisect_Range, XYZ_Names) of Real_Type;
   subtype NodesPerElement_Index_Array is
     Index_Array (NodesPerElement_Index_Type);
   type NodesPerElement_Index_Array_Array is array
     (Index_Type range <>) of NodesPerElement_Index_Array;
--     type NodesPerElement_XYZ_Real_Array is array
--       (NodesPerElement_Index_Type, XYZ_Names) of Real_Type;
--     type NodesPerElement_XYZ_Real_Array_Array is array
--       (NodesPerElement_Index_Type) of XYZ_Real_Array;
   type Face_NodesPerFace_NodesPerElement_Array is array
     (Face_Index_Type, NodesPerFace_Index_Type) of NodesPerElement_Index_Type;
--     subtype NodesPerFace_Real_Array is
--       Real_Array (NodesPerFace_Index_Type);
--     type NodesPerFace_XYZ_Real_Array_Array is array
--       (NodesPerFace_Index_Type) of XYZ_Real_Array;
--     type XYZ_NodesPerElement_Real_Array is array
--       (XYZ_Names, NodesPerElement_Index_Type) of Real_Type;
--     type XYZ_XEZ_Real_Array is array
--       (XYZ_Names, XEZ_Type) of Real_Type;



   type Domain_Record
     (numNode : Node_Index_Type;
      numElem : Element_Index_Type;
      numReg  : Region_Index_Type)
   is private;

   type Access_Domain is access Domain_Record;

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

   -- inline real4  CBRT(real4  arg) { return cbrtf(arg) ; }
   -- inline real8  CBRT(real8  arg) { return cbrt(arg) ; }
   -- inline real10 CBRT(real10 arg) { return cbrtl(arg) ; }

   -- inline real4  FABS(real4  arg) { return fabsf(arg) ; }
   -- inline real8  FABS(real8  arg) { return fabs(arg) ; }
   -- inline real10 FABS(real10 arg) { return fabsl(arg) ; }
   function FABS (This : in real4) return real4 is
      (abs (This));

   --- // Stuff needed for boundary conditions
   --- // 2 BCs on each of 6 hexahedral faces (12 bits)
   --x #define XI_M        0x00007
   --x #define XI_M_SYMM   0x00001
   --x #define XI_M_FREE   0x00002
   --x #define XI_M_COMM   0x00004
   XI_M      : constant := 2#0000_0000_0000_0000_0111#;
   XI_M_SYMM : constant := 2#0000_0000_0000_0000_0001#;
   XI_M_FREE : constant := 2#0000_0000_0000_0000_0010#;
   XI_M_COMM : constant := 2#0000_0000_0000_0000_0100#;

   --x #define XI_P        0x00038
   --x #define XI_P_SYMM   0x00008
   --x #define XI_P_FREE   0x00010
   --x #define XI_P_COMM   0x00020
   XI_P      : constant := 2#0000_0000_0000_0011_1000#;
   XI_P_SYMM : constant := 2#0000_0000_0000_0000_1000#;
   XI_P_FREE : constant := 2#0000_0000_0000_0001_0000#;
   XI_P_COMM : constant := 2#0000_0000_0000_0010_0000#;

   --x #define ETA_M       0x001c0
   --x #define ETA_M_SYMM  0x00040
   --x #define ETA_M_FREE  0x00080
   --x #define ETA_M_COMM  0x00100
   ETA_M      : constant := 2#0000_0000_0001_1100_0000#;
   ETA_M_SYMM : constant := 2#0000_0000_0000_0100_0000#;
   ETA_M_FREE : constant := 2#0000_0000_0000_1000_0000#;
   ETA_M_COMM : constant := 2#0000_0000_0001_0000_0000#;

   --x #define ETA_P       0x00e00
   --x #define ETA_P_SYMM  0x00200
   --x #define ETA_P_FREE  0x00400
   --x #define ETA_P_COMM  0x00800
   ETA_P      : constant := 2#0000_0000_1110_0000_0000#;
   ETA_P_SYMM : constant := 2#0000_0000_0010_0000_0000#;
   ETA_P_FREE : constant := 2#0000_0000_0100_0000_0000#;
   ETA_P_COMM : constant := 2#0000_0000_1000_0000_0000#;

   --x #define ZETA_M      0x07000
   --x #define ZETA_M_SYMM 0x01000
   --x #define ZETA_M_FREE 0x02000
   --x #define ZETA_M_COMM 0x04000
   ZETA_M      : constant := 2#0000_0111_0000_0000_0000#;
   ZETA_M_SYMM : constant := 2#0000_0001_0000_0000_0000#;
   ZETA_M_FREE : constant := 2#0000_0010_0000_0000_0000#;
   ZETA_M_COMM : constant := 2#0000_0100_0000_0000_0000#;

   --x #define ZETA_P      0x38000
   --x #define ZETA_P_SYMM 0x08000
   --x #define ZETA_P_FREE 0x10000
   --x #define ZETA_P_COMM 0x20000
   ZETA_P      : constant := 2#0011_1000_0000_0000_0000#;
   ZETA_P_SYMM : constant := 2#0000_1000_0000_0000_0000#;
   ZETA_P_FREE : constant := 2#0001_0000_0000_0000_0000#;
   ZETA_P_COMM : constant := 2#0010_0000_0000_0000_0000#;

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

   --    public:

   --    //
   --    // ALLOCATION
   --    //

   --    void AllocateNodePersistent(Int_t numNode) // Node-centered
   --    {
   --       m_x.resize(numNode);  // coordinates
   --       m_y.resize(numNode);
   --       m_z.resize(numNode);

   --       m_xd.resize(numNode); // velocities
   --       m_yd.resize(numNode);
   --       m_zd.resize(numNode);

   --       m_xdd.resize(numNode); // accelerations
   --       m_ydd.resize(numNode);
   --       m_zdd.resize(numNode);

   --       m_fx.resize(numNode);  // forces
   --       m_fy.resize(numNode);
   --       m_fz.resize(numNode);

   --       m_nodalMass.resize(numNode);  // mass
   --    }

   --    void AllocateElemPersistent(Int_t numElem) // Elem-centered
   --    {
   --       m_nodelist.resize(8*numElem);

   --       // elem connectivities through face
   --       m_lxim.resize(numElem);
   --       m_lxip.resize(numElem);
   --       m_letam.resize(numElem);
   --       m_letap.resize(numElem);
   --       m_lzetam.resize(numElem);
   --       m_lzetap.resize(numElem);

   --       m_elemBC.resize(numElem);

   --       m_e.resize(numElem);
   --       m_p.resize(numElem);

   --       m_q.resize(numElem);
   --       m_ql.resize(numElem);
   --       m_qq.resize(numElem);

   --       m_v.resize(numElem);

   --       m_volo.resize(numElem);
   --       m_delv.resize(numElem);
   --       m_vdov.resize(numElem);

   --       m_arealg.resize(numElem);

   --       m_ss.resize(numElem);

   --       m_elemMass.resize(numElem);
   --    }

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




   --    // Nodes on symmertry planes
   --    bool symmXempty()          { return m_symmX.empty(); }
   --    bool symmYempty()          { return m_symmY.empty(); }
   --    bool symmZempty()          { return m_symmZ.empty(); }

   --    //
   --    // Element-centered
   --    //
   --    Index_t*  nodelist(Index_t idx)    { return &m_nodelist[Index_t(8)*idx] ; }

   --    Index_t nodeElemCount(Index_t idx)
   --    { return m_nodeElemStart[idx+1] - m_nodeElemStart[idx] ; }

   --    Index_t *nodeElemCornerList(Index_t idx)
   --    { return &m_nodeElemCornerList[m_nodeElemStart[idx]] ; }



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

   --x   private:
private

   subtype Size_Type is Index_Type;

   type Gradient_Magnitude           is new Real_Type;
   type Grams                        is new Real_Type;
   type Joules                       is new Real_Type;
   type Meters                       is new Real_Type;
   type Meters_Per_Second            is new Real_Type;
   type Meters_Per_Second_Per_Second is new Real_Type;
   type Newtons                      is new Real_Type;

   type Acceleration_Vector is array (XYZ_Names) of Meters_Per_Second_Per_Second;
   type Force_Vector        is array (XYZ_Names) of Newtons;
   type Gradient_Type       is array (XEZ_Type)  of Gradient_Magnitude;
   type Location_Type       is array (XYZ_Names) of Meters;
   type XYZ_Size_Array      is array (XYZ_Names) of Size_Type;
   type Strain_Vector       is array (XYZ_Names) of Newtons;
   type Velocity_Vector     is array (XYZ_Names) of Meters_Per_Second;

   type Acceleration_Array  is array (Index_Type range <>) of Acceleration_Vector;
   type Force_Array         is array (Index_Type range <>) of Force_Vector;
   type Gradient_Array      is array (Index_Type range <>) of Gradient_Type;
   type Location_Array      is array (Index_Type range <>) of Location_Type;
   type Mass_Array          is array (Index_Type range <>) of Grams;
   type Strain_Array        is array (Index_Type range <>) of Strain_Vector;
   type Velocity_Array      is array (Index_Type range <>) of Velocity_Vector;

--     type NodesPerElement_Location_Array is array (Index_Type range <>) of Location_Type;

   package Element_Index_Vectors is
     new Ada.Containers.Vectors
       (Index_Type   => Index_Type,
        Element_Type => Element_Index_Type);

   type Region_Element_Vector_Array is array (Region_Index_Type range <>)
     of Element_Index_Vectors.Vector;

   ---    // Node-centered

   type Nodes_Record (numNode : Node_Index_Type) is record
      --x    std::vector<Real_t> m_x ;  /* coordinates */
      --x    std::vector<Real_t> m_y ;
      --x    std::vector<Real_t> m_z ;
      location : Location_Array (0..numNode);

      --x    std::vector<Real_t> m_xd ; /* velocities */
      --x    std::vector<Real_t> m_yd ;
      --x    std::vector<Real_t> m_zd ;
      velocity : Velocity_Array (0..numNode);

      --x    std::vector<Real_t> m_xdd ; /* accelerations */
      --x    std::vector<Real_t> m_ydd ;
      --x    std::vector<Real_t> m_zdd ;
      acceleration : Acceleration_Array (0..numNode);

      --x    std::vector<Real_t> m_fx ;  /* forces */
      --x    std::vector<Real_t> m_fy ;
      --x    std::vector<Real_t> m_fz ;
      force : Force_Array (0..numNode);

      --x    std::vector<Real_t> m_nodalMass ;  /* mass */
      mass : Mass_Array (0..numNode);

      --x    std::vector<Index_t> m_symmX ;  /* symmetry plane nodesets */
      --x    std::vector<Index_t> m_symmY ;
      --x    std::vector<Index_t> m_symmZ ;
--        symmX : Index_Vector;
--        symmY : Index_Vector;
--        symmZ : Index_Vector;

   end record;

   type Elements_Record
     (numElem : Element_Index_Type;
      numReg  : Region_Index_Type)
   is record
      ---    // Element-centered

      ---    // Region information

      --x    Int_t    m_numReg ;
      -- Moved to discriminant.
      --x    Int_t    m_cost; //imbalance cost
      --x    Index_t *m_regElemSize ;   // Size of region sets
      --x    Index_t *m_regNumList ;    // Region number per domain element
      --x    Index_t **m_regElemlist ;  // region indexset
      imbalance_cost : Int_t;
      regElemSize    : Index_Array (0..numElem);
      regNumList     : Index_Array (0..numElem);
      regElemList    : Region_Element_Vector_Array (0..numReg);

      --x    std::vector<Index_t>  m_nodelist ;     /* elemToNode connectivity */
      node_indexes : NodesPerElement_Index_Array_Array (0..numElem);

      --x    std::vector<Index_t>  m_lxim ;  /* element connectivity across each face */
      --x    std::vector<Index_t>  m_lxip ;
      --x    std::vector<Index_t>  m_letam ;
      --x    std::vector<Index_t>  m_letap ;
      --x    std::vector<Index_t>  m_lzetam ;
      --x    std::vector<Index_t>  m_lzetap ;
      --! Or should this be XEZ_Face_Index_Array_Array?
--        connectivity_m : XEZ_Real_Array_Array (1..numNode);
--        connectivity_p : XEZ_Real_Array_Array (1..numNode);

   --    // elem face symm/free-surface flag
      --x    std::vector<Int_t>    m_elemBC ;  /* symmetry/free-surface flags for each elem face */
--      elemBC : Face_Int_t_Array_Array (0..numElem-1);

      --x    std::vector<Real_t> m_dxx ;  /* principal strains -- temporary */
      --x    std::vector<Real_t> m_dyy ;
      --x    std::vector<Real_t> m_dzz ;
      principal_strain : Strain_Array (0..numElem);

      --x    std::vector<Real_t> m_delv_xi ;    /* velocity gradient -- temporary */
      --x    std::vector<Real_t> m_delv_eta ;
      --x    std::vector<Real_t> m_delv_zeta ;
      velocity_gradient : Gradient_Array (0..numElem);

      --    // Position gradient - temporary
      --x    std::vector<Real_t> m_delx_xi ;    /* coordinate gradient -- temporary */
      --x    std::vector<Real_t> m_delx_eta ;
      --x    std::vector<Real_t> m_delx_zeta ;
      position_gradient : Gradient_Array (0..numElem);

      --x    std::vector<Real_t> m_e ;   /* energy */
      energy : Joules;

      --x    std::vector<Real_t> m_p ;   /* pressure */
      --x    std::vector<Real_t> m_q ;   /* q */
      --x    std::vector<Real_t> m_ql ;  /* linear term for q */
      --x    std::vector<Real_t> m_qq ;  /* quadratic term for q */
      pressure    : Real_Array (0..numElem);
      --    // Artificial viscosity
      q           : Real_Array (0..numElem);
      q_linear    : Real_Array (0..numElem);
      q_quadratic : Real_Array (0..numElem);

      --x    std::vector<Real_t> m_v ;     /* relative volume */
      --x    std::vector<Real_t> m_volo ;  /* reference volume */
      --x    std::vector<Real_t> m_vnew ;  /* new relative volume -- temporary */
      --x    std::vector<Real_t> m_delv ;  /* m_vnew - m_v */
      --x    std::vector<Real_t> m_vdov ;  /* volume derivative over volume */
      relative_volume               : Real_Array (0..numElem);
      reference_volume              : Real_Array (0..numElem);
      new_relative_volume           : Real_Array (0..numElem);
      relative_volume_delta         : Real_Array (0..numElem);
      volume_derivative_over_volume : Real_Array (0..numElem);

      --x    std::vector<Real_t> m_arealg ;  /* characteristic length of an element */
      characteristic_length : Real_Array (0..numElem);

      --x    std::vector<Real_t> m_ss ;      /* "sound speed" */
      sound_speed : Velocity_Array (0..numElem);

      --x    std::vector<Real_t> m_elemMass ;  /* mass */
      --    // Element mass
      elemMass : Mass_Array (0..numElem);
   end record;

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
      dtcourant       : Real_Type;
      dthydro         : Real_Type;
      cycle           : Int_t;
      dtfixed         : Real_Type;
      time            : Real_Type;
      deltatime       : Real_Type;
      deltatimemultlb : Real_Type;
      deltatimemultub : Real_Type;
      dtmax           : Real_Type;
      stoptime        : Real_Type;

      --x    Int_t   m_numRanks ;
      numRanks : Int_t;

      --x    Index_t m_colLoc ;
      --x    Index_t m_rowLoc ;
      --x    Index_t m_planeLoc ;
      --x    Index_t m_tp ;
      colLoc   : Index_Type;
      rowLoc   : Index_Type;
      planeLoc : Index_Type;
      tp       : Index_Type;
   end record;

   type Parameters_Record is record
      ---    // Parameters

      ---    // Cutoffs (treat as constants)

      ---    const Real_t  m_e_cut ;             // energy tolerance
      ---    const Real_t  m_p_cut ;             // pressure tolerance
      ---    const Real_t  m_q_cut ;             // q tolerance
      ---    const Real_t  m_v_cut ;             // relative volume tolerance
      ---    const Real_t  m_u_cut ;             // velocity tolerance
      energy_tolerance          : Joules;
      pressure_tolerance        : Real_Type;
      q_tolerance               : Real_Type;
      relative_volume_tolerance : Real_Type;
      velocity_tolerance        : Meters_Per_Second;

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
      ss4o3              : Real_Type;
      qstop              : Real_Type;
      monoq_max_slope    : Real_Type;
      monoq_limiter_mult : Real_Type;
      qlc_monoq          : Real_Type;
      qqc_monoq          : Real_Type;
      qqc                : Real_Type;
      eosvmax            : Real_Type;
      eosvmin            : Real_Type;
      pmin               : Real_Type;
      dvovmax            : Real_Type;
      refdens            : Real_Type;


      --x    Index_t m_sizeX ;
      --x    Index_t m_sizeY ;
      --x    Index_t m_sizeZ ;
      size : XYZ_Size_Array;
      --x    Index_t m_numElem ;
      --x    Index_t m_numNode ;
      -- Moved to discriminant.

      --x    Index_t m_maxPlaneSize ;
      --x    Index_t m_maxEdgeSize ;
      maxPlaneSize : Index_Type;
      maxEdgeSize  : Index_Type;

   --    // OMP hack
   --    Index_t *m_nodeElemStart ;
   --    Index_t *m_nodeElemCornerList ;

      ---    // Used in setup
      --x    Index_t m_rowMin, m_rowMax;
      --x    Index_t m_colMin, m_colMax;
      --x    Index_t m_planeMin, m_planeMax ;
      rowMin   : Index_Type;
      rowMax   : Index_Type;
      colMin   : Index_Type;
      colMax   : Index_Type;
      planeMin : Index_Type;
      planeMax : Index_Type;
   end record;

   type Domain_Record
     (numNode : Node_Index_Type;
      numElem : Element_Index_Type;
      numReg  : Region_Index_Type)
   is record
      nodes      : Nodes_Record (numNode);
      elements   : Elements_Record (numElem, numReg);
      variables  : Variables_Record;
      parameters : Parameters_Record;
   end record;

   -- typedef Real_t &(Domain::* Domain_member )(Index_t) ;
   -- (Address of (pointer to) a Domain class function (accessor):
   --   that takes an Index_t and returns the Real_t Domain component value)
   type Accessor_access is access function (Index : in Index_Type) return Real_Type;
   type Domain_member is access all Real_Type;

   --- /******************************************/

   --- /* Work Routines */

   --x static inline
   --x void TimeIncrement(Domain& domain)
   procedure TimeIncrement (domain : not null access Domain_Record)
   with Inline;

   --x static inline
   --x void LagrangeLeapFrog(Domain& domain)
   procedure LagrangeLeapFrog (domain : not null access Domain_Record)
   with Inline;

end LULESH;
