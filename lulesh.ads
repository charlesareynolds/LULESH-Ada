with Ada.Calendar;
with Ada.Containers.Vectors;
with Ada.Numerics.Generic_Elementary_Functions;
with Interfaces.C;

package LULESH is

   package AC renames Ada.Calendar;
   package IC renames Interfaces.C;

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
   type real4   is new Interfaces.C.C_float;
   type real8   is new Interfaces.C.double;
   type real10  is new Interfaces.C.long_double;
   type Index_t is new Interfaces.C.int range 0..Interfaces.C.int'Last;
   type Real_t  is new real8;
   type Int_t   is new Interfaces.C.int;

   ---------------------------------------

   type Real_t_Access is access all Real_t;
   type Real_t_Array is array (Index_t range <>) of Real_t;
   type Real_t_Array_Access is access Real_t_Array;
   type Real_t_Array_Access_nn is not null access Real_t_Array;

   type Index_t_Array is array (Index_t range <>) of Index_t;

   type Dimensions is (X, Y, Z);
   subtype XYZ_Type is Dimensions;
   subtype Coord_Names is Dimensions;

   type XEZ_Type is (Xi, Eta, Zeta);
   Et : constant XEZ_Type := Eta;
   Ze : constant XEZ_Type := Zeta;

   subtype Elem_Node_Range is Index_t range 0..7;
   subtype Face_Node_Range is Index_t range 0..3;
   subtype Bisect_Range is Index_t range 0..1;
   subtype Face_Range is Index_t range 1..6;

   subtype Elem_Node_Index_Array is Index_t_Array (Elem_Node_Range);
   subtype Elem_Node_Value_Array is Real_t_Array (Elem_Node_Range);
   subtype Face_Node_Value_Array is Real_t_Array (Face_Node_Range);
   type Coord_Value_Array is array (Coord_Names) of Real_t;

   type Coord_Elem_Node_Value_Array is array
     (Coord_Names, Elem_Node_Range) of Real_t;
   type Elem_Node_Coord_Value_Array is array
     (Elem_Node_Range, Coord_Names) of Real_t;
   type Elem_Node_Coord_Value_Array_Array is array
     (Elem_Node_Range) of Coord_Value_Array;

   type Coord_Face_Node_Value_Array is array
     (Coord_Names, Face_Node_Range) of Real_t;
   type Coord_Face_Node_Value_Array_Array is array
     (Coord_Names) of Face_Node_Value_Array;
   type Face_Node_Coord_Value_Array is array
     (Face_Node_Range, Coord_Names) of Real_t;
   type Face_Node_Coord_Value_Array_Array is array
     (Face_Node_Range) of Coord_Value_Array;

   type Coord_XEZ_Value_Array is array
     (Coord_Names, XEZ_Type) of Real_t;

   type Bisect_Coord_Value_Array is array
     (Bisect_Range, Coord_Names) of Real_t;
   subtype Area_Array is Coord_Value_Array;

   type Face_Face_Node_Elem_Node_Array is array
     (Face_Range, Face_Node_Range) of Elem_Node_Range;

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
   CACHE_COHERENCE_PAD_REAL : constant := 128 / (Real_t'Size / Byte_Size);

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
   type Domain_Record
     (numNode : Index_t;
      numElem : Index_t) is tagged private;
   type Access_Domain is access Domain_Record;

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

   --    //
   --    // ACCESSORS
   --    //

   --    // Node-centered

   --    // Nodal coordinates
   --    Real_t& x(Index_t idx)    { return m_x[idx] ; }
--     function x
--       (Self : not null access Domain_Record;
--        idx  : in Index_t)
--        return Real_t_Access is
--       (Self.x (idx)'Access);
   --    Real_t& y(Index_t idx)    { return m_y[idx] ; }
   --    Real_t& z(Index_t idx)    { return m_z[idx] ; }

   --    // Nodal velocities
   --    Real_t& xd(Index_t idx)   { return m_xd[idx] ; }
   --    Real_t& yd(Index_t idx)   { return m_yd[idx] ; }
   --    Real_t& zd(Index_t idx)   { return m_zd[idx] ; }

   --    // Nodal accelerations
   --    Real_t& xdd(Index_t idx)  { return m_xdd[idx] ; }
   --    Real_t& ydd(Index_t idx)  { return m_ydd[idx] ; }
   --    Real_t& zdd(Index_t idx)  { return m_zdd[idx] ; }

   --    // Nodal forces
   --    Real_t& fx(Index_t idx)   { return m_fx[idx] ; }
   --    Real_t& fy(Index_t idx)   { return m_fy[idx] ; }
   --    Real_t& fz(Index_t idx)   { return m_fz[idx] ; }

   --    // Nodal mass
   --    Real_t& nodalMass(Index_t idx) { return m_nodalMass[idx] ; }

   --    // Nodes on symmertry planes
   --    Index_t symmX(Index_t idx) { return m_symmX[idx] ; }
   --    Index_t symmY(Index_t idx) { return m_symmY[idx] ; }
   --    Index_t symmZ(Index_t idx) { return m_symmZ[idx] ; }
   --    bool symmXempty()          { return m_symmX.empty(); }
   --    bool symmYempty()          { return m_symmY.empty(); }
   --    bool symmZempty()          { return m_symmZ.empty(); }

   --    //
   --    // Element-centered
   --    //
   --    Index_t&  regElemSize(Index_t idx) { return m_regElemSize[idx] ; }
   --    Index_t&  regNumList(Index_t idx) { return m_regNumList[idx] ; }
   --    Index_t*  regNumList()            { return &m_regNumList[0] ; }
   --    Index_t*  regElemlist(Int_t r)    { return m_regElemlist[r] ; }
   --    Index_t&  regElemlist(Int_t r, Index_t idx) { return m_regElemlist[r][idx] ; }

   --    Index_t*  nodelist(Index_t idx)    { return &m_nodelist[Index_t(8)*idx] ; }

   --    // elem connectivities through face
   --    Index_t&  lxim(Index_t idx) { return m_lxim[idx] ; }
   --    Index_t&  lxip(Index_t idx) { return m_lxip[idx] ; }
   --    Index_t&  letam(Index_t idx) { return m_letam[idx] ; }
   --    Index_t&  letap(Index_t idx) { return m_letap[idx] ; }
   --    Index_t&  lzetam(Index_t idx) { return m_lzetam[idx] ; }
   --    Index_t&  lzetap(Index_t idx) { return m_lzetap[idx] ; }

   --    // elem face symm/free-surface flag
   --    Int_t&  elemBC(Index_t idx) { return m_elemBC[idx] ; }

   --    // Principal strains - temporary
   --    Real_t& dxx(Index_t idx)  { return m_dxx[idx] ; }
   --    Real_t& dyy(Index_t idx)  { return m_dyy[idx] ; }
   --    Real_t& dzz(Index_t idx)  { return m_dzz[idx] ; }

   --    // Velocity gradient - temporary
   --    Real_t& delv_xi(Index_t idx)    { return m_delv_xi[idx] ; }
   --    Real_t& delv_eta(Index_t idx)   { return m_delv_eta[idx] ; }
   --    Real_t& delv_zeta(Index_t idx)  { return m_delv_zeta[idx] ; }

   --    // Position gradient - temporary
   --    Real_t& delx_xi(Index_t idx)    { return m_delx_xi[idx] ; }
   --    Real_t& delx_eta(Index_t idx)   { return m_delx_eta[idx] ; }
   --    Real_t& delx_zeta(Index_t idx)  { return m_delx_zeta[idx] ; }

   --    // Energy
   --    Real_t& e(Index_t idx)          { return m_e[idx] ; }

   --    // Pressure
   --    Real_t& p(Index_t idx)          { return m_p[idx] ; }

   --    // Artificial viscosity
   --    Real_t& q(Index_t idx)          { return m_q[idx] ; }

   --    // Linear term for q
   --    Real_t& ql(Index_t idx)         { return m_ql[idx] ; }
   --    // Quadratic term for q
   --    Real_t& qq(Index_t idx)         { return m_qq[idx] ; }

   --    // Relative volume
   --    Real_t& v(Index_t idx)          { return m_v[idx] ; }
   --    Real_t& delv(Index_t idx)       { return m_delv[idx] ; }

   --    // Reference volume
   --    Real_t& volo(Index_t idx)       { return m_volo[idx] ; }

   --    // volume derivative over volume
   --    Real_t& vdov(Index_t idx)       { return m_vdov[idx] ; }

   --    // Element characteristic length
   --    Real_t& arealg(Index_t idx)     { return m_arealg[idx] ; }

   --    // Sound speed
   --    Real_t& ss(Index_t idx)         { return m_ss[idx] ; }

   --    // Element mass
   --    Real_t& elemMass(Index_t idx)  { return m_elemMass[idx] ; }

   --    Index_t nodeElemCount(Index_t idx)
   --    { return m_nodeElemStart[idx+1] - m_nodeElemStart[idx] ; }

   --    Index_t *nodeElemCornerList(Index_t idx)
   --    { return &m_nodeElemCornerList[m_nodeElemStart[idx]] ; }

   --    // Parameters

   --    // Cutoffs
   --    Real_t u_cut() const               { return m_u_cut ; }
   --    Real_t e_cut() const               { return m_e_cut ; }
   --    Real_t p_cut() const               { return m_p_cut ; }
   --    Real_t q_cut() const               { return m_q_cut ; }
   --    Real_t v_cut() const               { return m_v_cut ; }

   --    // Other constants (usually are settable via input file in real codes)
   --    Real_t hgcoef() const              { return m_hgcoef ; }
   --    Real_t qstop() const               { return m_qstop ; }
   --    Real_t monoq_max_slope() const     { return m_monoq_max_slope ; }
   --    Real_t monoq_limiter_mult() const  { return m_monoq_limiter_mult ; }
   --    Real_t ss4o3() const               { return m_ss4o3 ; }
   --    Real_t qlc_monoq() const           { return m_qlc_monoq ; }
   --    Real_t qqc_monoq() const           { return m_qqc_monoq ; }
   --    Real_t qqc() const                 { return m_qqc ; }

   --    Real_t eosvmax() const             { return m_eosvmax ; }
   --    Real_t eosvmin() const             { return m_eosvmin ; }
   --    Real_t pmin() const                { return m_pmin ; }
   --    Real_t emin() const                { return m_emin ; }
   --    Real_t dvovmax() const             { return m_dvovmax ; }
   --    Real_t refdens() const             { return m_refdens ; }

   --    // Timestep controls, etc...
   --    Real_t& time()                 { return m_time ; }
   --    Real_t& deltatime()            { return m_deltatime ; }
   --    Real_t& deltatimemultlb()      { return m_deltatimemultlb ; }
   --    Real_t& deltatimemultub()      { return m_deltatimemultub ; }
   --    Real_t& stoptime()             { return m_stoptime ; }
   --    Real_t& dtcourant()            { return m_dtcourant ; }
   --    Real_t& dthydro()              { return m_dthydro ; }
   --    Real_t& dtmax()                { return m_dtmax ; }
   --    Real_t& dtfixed()              { return m_dtfixed ; }

   --    Int_t&  cycle()                { return m_cycle ; }
   --    Index_t&  numRanks()           { return m_numRanks ; }

   --    Index_t&  colLoc()             { return m_colLoc ; }
   --    Index_t&  rowLoc()             { return m_rowLoc ; }
   --    Index_t&  planeLoc()           { return m_planeLoc ; }
   --    Index_t&  tp()                 { return m_tp ; }

   --    Index_t&  sizeX()              { return m_sizeX ; }
   --    Index_t&  sizeY()              { return m_sizeY ; }
   --    Index_t&  sizeZ()              { return m_sizeZ ; }
   --    Index_t&  numReg()             { return m_numReg ; }
   --    Int_t&  cost()             { return m_cost ; }
   --    Index_t&  numElem()            { return m_numElem ; }
   --    Index_t&  numNode()            { return m_numNode ; }

   --    Index_t&  maxPlaneSize()       { return m_maxPlaneSize ; }
   --    Index_t&  maxEdgeSize()        { return m_maxEdgeSize ; }

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

   subtype C_STL_Vector_Index_Type is Natural;

--     package Real_t_Vectors is new Ada.Containers.Vectors
--       (Index_Type   => C_STL_Vector_Index_Type,
--        Element_Type => Real_t);
--     -- To make operations directly visible:
--     type Real_t_Vector is new Real_t_Vectors.Vector with null record;
   type Real_t_Vector is
     array (Index_t range <>) of Real_t;

   package Int_t_Vectors is new Ada.Containers.Vectors
     (Index_Type   => C_STL_Vector_Index_Type,
      Element_Type => Int_t);
   subtype Int_t_Vector is Int_t_Vectors.Vector;

   package Index_t_Vectors is new Ada.Containers.Vectors
     (Index_Type   => C_STL_Vector_Index_Type,
      Element_Type => Index_t);
   subtype Index_t_Vector is Index_t_Vectors.Vector;

   type Elem_Node_Index_Array_Vector is
     array (Index_t range <>) of Elem_Node_Index_Array;

   type Coord_Value_Array_Vector is
     array (Index_t range <>) of Coord_Value_Array;

   ---    //
   ---    // IMPLEMENTATION
   ---    //
   type Domain_Record
     (numNode : Index_t;
      numElem : Index_t) is tagged record
      ---    /* Node-centered */
      --x    std::vector<Real_t> m_x ;  /* coordinates */
      --x    std::vector<Real_t> m_y ;
      --x    std::vector<Real_t> m_z ;
      coords : Coord_Value_Array_Vector (1..numNode);

--        x : Real_t_Vector (1..numNode);
--        y : Real_t_Vector (1..numNode);
--        z : Real_t_Vector (1..numNode);

      --x    std::vector<Real_t> m_xd ; /* velocities */
      --x    std::vector<Real_t> m_yd ;
      --x    std::vector<Real_t> m_zd ;
      xd : Real_t_Vector (1..numNode);
      yd : Real_t_Vector (1..numNode);
      zd : Real_t_Vector (1..numNode);

      --x    std::vector<Real_t> m_xdd ; /* accelerations */
      --x    std::vector<Real_t> m_ydd ;
      --x    std::vector<Real_t> m_zdd ;
      xdd : Real_t_Vector (1..numNode);
      ydd : Real_t_Vector (1..numNode);
      zdd : Real_t_Vector (1..numNode);

      --x    std::vector<Real_t> m_fx ;  /* forces */
      --x    std::vector<Real_t> m_fy ;
      --x    std::vector<Real_t> m_fz ;
      fx : Real_t_Vector (1..numNode);
      fy : Real_t_Vector (1..numNode);
      fz : Real_t_Vector (1..numNode);

      --x    std::vector<Real_t> m_nodalMass ;  /* mass */
      nodalMass : Real_t_Vector (1..numNode);

      --x    std::vector<Index_t> m_symmX ;  /* symmetry plane nodesets */
      --x    std::vector<Index_t> m_symmY ;
      --x    std::vector<Index_t> m_symmZ ;
      symmX : Index_t_Vector;
      symmY : Index_t_Vector;
      symmZ : Index_t_Vector;

      ---    // Element-centered

      ---    // Region information
      --x    Int_t    m_numReg ;
      --x    Int_t    m_cost; //imbalance cost
   --    Index_t *m_regElemSize ;   // Size of region sets
   --    Index_t *m_regNumList ;    // Region number per domain element
   --    Index_t **m_regElemlist ;  // region indexset
      numReg : Int_t;
      cost   : Int_t;

      --x    std::vector<Index_t>  m_nodelist ;     /* elemToNode connectivity */
      nodelist : Elem_Node_Index_Array_Vector (0..numElem);

      --x    std::vector<Index_t>  m_lxim ;  /* element connectivity across each face */
      --x    std::vector<Index_t>  m_lxip ;
      --x    std::vector<Index_t>  m_letam ;
      --x    std::vector<Index_t>  m_letap ;
      --x    std::vector<Index_t>  m_lzetam ;
      --x    std::vector<Index_t>  m_lzetap ;
      lxim   : Index_t_Vector;
      lxip   : Index_t_Vector;
      letam  : Index_t_Vector;
      letap  : Index_t_Vector;
      lzetam : Index_t_Vector;
      lzetap : Index_t_Vector;

      --x    std::vector<Int_t>    m_elemBC ;  /* symmetry/free-surface flags for each elem face */
      elemBC : Int_t_Vector;

      --x    std::vector<Real_t> m_dxx ;  /* principal strains -- temporary */
      --x    std::vector<Real_t> m_dyy ;
      --x    std::vector<Real_t> m_dzz ;
      dxx : Real_t_Vector (0..numElem);
      dyy : Real_t_Vector (0..numElem);
      dzz : Real_t_Vector (0..numElem);

      --x    std::vector<Real_t> m_delv_xi ;    /* velocity gradient -- temporary */
      --x    std::vector<Real_t> m_delv_eta ;
      --x    std::vector<Real_t> m_delv_zeta ;
      delv_xi   : Real_t_Vector (0..numElem);
      delv_eta  : Real_t_Vector (0..numElem);
      delv_zeta : Real_t_Vector (0..numElem);

      --x    std::vector<Real_t> m_delx_xi ;    /* coordinate gradient -- temporary */
      --x    std::vector<Real_t> m_delx_eta ;
      --x    std::vector<Real_t> m_delx_zeta ;
      delx_xi   : Real_t_Vector (0..numElem);
      delx_eta  : Real_t_Vector (0..numElem);
      delx_zeta : Real_t_Vector (0..numElem);

      --x    std::vector<Real_t> m_e ;   /* energy */
      e : Real_t;

      --x    std::vector<Real_t> m_p ;   /* pressure */
      --x    std::vector<Real_t> m_q ;   /* q */
      --x    std::vector<Real_t> m_ql ;  /* linear term for q */
      --x    std::vector<Real_t> m_qq ;  /* quadratic term for q */
      p  : Real_t_Vector (0..numElem);
      q  : Real_t_Vector (0..numElem);
      ql : Real_t_Vector (0..numElem);
      qq : Real_t_Vector (0..numElem);

      --x    std::vector<Real_t> m_v ;     /* relative volume */
      --x    std::vector<Real_t> m_volo ;  /* reference volume */
      --x    std::vector<Real_t> m_vnew ;  /* new relative volume -- temporary */
      --x    std::vector<Real_t> m_delv ;  /* m_vnew - m_v */
      --x    std::vector<Real_t> m_vdov ;  /* volume derivative over volume */
      v    : Real_t_Vector (0..numElem);
      volo : Real_t_Vector (0..numElem);
      vnew : Real_t_Vector (0..numElem);
      delv : Real_t_Vector (0..numElem);
      vdov : Real_t_Vector (0..numElem);

      --x    std::vector<Real_t> m_arealg ;  /* characteristic length of an element */
      arealg : Real_t_Vector (0..numElem);

      --x    std::vector<Real_t> m_ss ;      /* "sound speed" */
      ss : Real_t_Vector (0..numElem);

      --x    std::vector<Real_t> m_elemMass ;  /* mass */
      elemMass : Real_t_Vector (0..numElem);

      ---    // Cutoffs (treat as constants)
      ---    const Real_t  m_e_cut ;             // energy tolerance
      ---    const Real_t  m_p_cut ;             // pressure tolerance
      ---    const Real_t  m_q_cut ;             // q tolerance
      ---    const Real_t  m_v_cut ;             // relative volume tolerance
      ---    const Real_t  m_u_cut ;             // velocity tolerance
      e_cut : Real_t;
      p_cut : Real_t;
      q_cut : Real_t;
      v_cut : Real_t;
      u_cut : Real_t;

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
      hgcoef             : Real_t;
      ss4o3              : Real_t;
      qstop              : Real_t;
      monoq_max_slope    : Real_t;
      monoq_limiter_mult : Real_t;
      qlc_monoq          : Real_t;
      qqc_monoq          : Real_t;
      qqc                : Real_t;
      eosvmax            : Real_t;
      eosvmin            : Real_t;
      pmin               : Real_t;
      dvovmax            : Real_t;
      refdens            : Real_t;

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
      dtcourant       : Real_t;
      dthydro         : Real_t;
      cycle           : Int_t;
      dtfixed         : Real_t;
      time            : Real_t;
      deltatime       : Real_t;
      deltatimemultlb : Real_t;
      deltatimemultub : Real_t;
      dtmax           : Real_t;
      stoptime        : Real_t;

      --x    Int_t   m_numRanks ;
      numRanks : Int_t;

      --x    Index_t m_colLoc ;
      --x    Index_t m_rowLoc ;
      --x    Index_t m_planeLoc ;
      --x    Index_t m_tp ;
      colLoc   : Index_t;
      rowLoc   : Index_t;
      planeLoc : Index_t;
      tp       : Index_t;

      --x    Index_t m_sizeX ;
      --x    Index_t m_sizeY ;
      --x    Index_t m_sizeZ ;
      --x    Index_t m_numElem ;
      --x    Index_t m_numNode ;
      sizeX   : Index_t;
      sizeY   : Index_t;
      sizeZ   : Index_t;
      -- Moved to discriminant:
--        numElem : Index_t;
--        numNode : Index_t;

      --x    Index_t m_maxPlaneSize ;
      --x    Index_t m_maxEdgeSize ;
      maxPlaneSize : Index_t;
      maxEdgeSize  : Index_t;

   --    // OMP hack
   --    Index_t *m_nodeElemStart ;
   --    Index_t *m_nodeElemCornerList ;

      ---    // Used in setup
      --x    Index_t m_rowMin, m_rowMax;
      --x    Index_t m_colMin, m_colMax;
      --x    Index_t m_planeMin, m_planeMax ;
      rowMin   : Index_t;
      rowMax   : Index_t;
      colMin   : Index_t;
      colMax   : Index_t;
      planeMin : Index_t;
      planeMax : Index_t;

      --x } ;
   end record;

   -- typedef Real_t &(Domain::* Domain_member )(Index_t) ;
   -- (Address of (pointer to) a Domain class function (accessor):
   --   that takes an Index_t and returns the Real_t Domain component value)
   type Accessor_access is access function (Index : in Index_t) return Real_t;
   type Domain_member is access all Real_t;

   --- /******************************************/

   --- /* Work Routines */

   --x static inline
   --x void TimeIncrement(Domain& domain)
   procedure TimeIncrement (domain : not null access Domain_Record);

   --x static inline
   --x void LagrangeLeapFrog(Domain& domain)
   procedure LagrangeLeapFrog (domain : not null access Domain_Record);

end LULESH;
