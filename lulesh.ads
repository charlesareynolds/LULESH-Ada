with Ada.Calendar;
with Ada.Containers.Bounded_Vectors;
with Ada.Exceptions;
with Ada.Numerics.Generic_Elementary_Functions;
with Ada.Numerics.Generic_Real_Arrays;
with Ada.Real_Time;
with Ada.Unchecked_Deallocation;
with Interfaces;
with MPI;
with OMP;
with System.Dim.Mks;
with System.Dim.Mks_IO;

-- For types and operations:
use System.Dim.Mks;

package LULESH is

   -----------------------------------------------------------------------------
   --- General declarations:
   -----------------------------------------------------------------------------
   
   package AC    renames Ada.Calendar;
   package ART   renames Ada.Real_Time;
   package MKS   renames System.Dim.Mks;
   package MKSIO renames System.Dim.Mks_IO;

   Usage_Error  : exception;
   Coding_Error : exception;

   --- #if !defined(USE_MPI)
   --- # error "You should specify USE_MPI=0 or USE_MPI=1 on the compile line"
   --- #endif
   USE_MPI : constant Boolean := False;
   --- // OpenMP will be compiled in if this flag is set to 1 AND the compiler beging
   --- // used supports it (i.e. the _OPENMP symbol is defined)
   --x #define USE_OMP 1
   USE_OMP                  : constant Boolean := False;
   COMPILER_SUPPORTS_OPENMP : constant Boolean := False;
   --x #if USE_MPI
   --x #include <mpi.h>
   --- /*
   ---    define one of these three symbols:
   ---    SEDOV_SYNC_POS_VEL_NONE
   ---    SEDOV_SYNC_POS_VEL_EARLY
   ---    SEDOV_SYNC_POS_VEL_LATE
   --- */
   --x #define SEDOV_SYNC_POS_VEL_EARLY 1
   --x #endif
   SEDOV_SYNC_POS_VEL_NONE  : constant Boolean := False;
   SEDOV_SYNC_POS_VEL_EARLY : constant Boolean := USE_MPI;
   SEDOV_SYNC_POS_VEL_LATE  : constant Boolean := False;
   --- // MPI Message Tags
   --x #define MSG_COMM_SBN      1024
   --x #define MSG_SYNC_POS_VEL  2048
   --x #define MSG_MONOQ         3072
   MSG_COMM_SBN     : constant := 2#0010_0000_0000#;
   MSG_SYNC_POS_VEL : constant := 2#0100_0000_0000#;
   MSG_MONOQ        : constant := 2#0110_0000_0000#;
   --x #define MAX_FIELDS_PER_MPI_COMM 6
   MAX_FIELDS_PER_MPI_COMM : constant := 6;

   --x enum { VolumeError = -1, QStopError = -2 } ;
   VolumeError : constant MPI.Errorcode_Type := MPI.Errorcode_Type (-1);
   QStop_Error : constant MPI.Errorcode_Type := MPI.Errorcode_Type (-2);

   ---    // Maximum number of block neighbors
   --- // 6 faces + 12 edges + 8 corners
   MAXIMUM_BLOCK_NEIGHBORS : constant := 26;
   type MPI_Request_Array is array (0 .. MAXIMUM_BLOCK_NEIGHBORS) of MPI.Request;
   
   -----------------------------------------------------------------------------
   --- Unconstrained base types:
   -----------------------------------------------------------------------------

   --- #include <math.h>
   --- #include <vector>
   --- //**************************************************
   --- // Allow flexibility for arithmetic representations
   --- //**************************************************
   --- // Precision specification
   --x typedef float        real4 ;
   --x typedef double       real8 ;
   --x typedef long double  real10 ;  // 10 bytes on x86
   --x typedef int    Index_t ; // array subscript and loop index
   --x typedef real8  Real_t ;  // floating point representation
   --x typedef int    Int_t ;   // integer representation
   type Real4        is new Interfaces.IEEE_Float_32;
   type Real8        is new Interfaces.IEEE_Float_64;
   subtype Real_Type is Real8;
   type Real10       is new Interfaces.IEEE_Extended_Float;
   type Index_Type   is new Natural;
   type Balance_Type is new Natural;
   type Cost_Type    is new Natural;

   --x #define MAX(a, b) ( ((a) > (b)) ? (a) : (b))
   function MAX (L, R : in Index_Type) return Index_Type is
     (if L > R then L else R)
   with Inline;

   --- // Assume 128 byte coherence
   --- // Assume Real_t is an "integral power of 2" bytes wide
   --x #define CACHE_COHERENCE_PAD_REAL (128 / sizeof(Real_t))
   Byte_Size                : constant := 8;
   Real_Byte_Size           : constant := Real_Type'Size / Byte_Size;
   CACHE_COHERENCE_PAD_REAL : constant := 128 / Real_Byte_Size;

   --x #define CACHE_ALIGN_REAL(n) \
   --x    (((n) + (CACHE_COHERENCE_PAD_REAL - 1)) & ~(CACHE_COHERENCE_PAD_REAL-1))
   --- Declared above derived Index_Type declarations so they inherit this:
   --- Round up to a whole multiple of CACHE_COHERENCE_PAD_REAL:
   function CACHE_ALIGN_REAL (This : in Index_Type) return Index_Type is
     (((This + CACHE_COHERENCE_PAD_REAL - 1) / CACHE_COHERENCE_PAD_REAL) * 
        CACHE_COHERENCE_PAD_REAL);

   --x inline real4  SQRT(real4  arg) { return sqrtf(arg) ; }
   --x inline real8  SQRT(real8  arg) { return sqrt(arg) ; }
   --x inline real10 SQRT(real10 arg) { return sqrtl(arg) ; }
   package Real4_Elementary_Functions is new Ada.Numerics
     .Generic_Elementary_Functions (Real4);
   package Real8_Elementary_Functions is new Ada.Numerics
     .Generic_Elementary_Functions (Real8);
   package Real10_Elementary_Functions is new Ada.Numerics
     .Generic_Elementary_Functions (Real10);
   package Mks_Type_Elementary_Functions is new Ada.Numerics
     .Generic_Elementary_Functions (Mks_Type);
   function SQRT
     (This : in Real4) return Real4 renames
     Real4_Elementary_Functions.Sqrt;
   function SQRT
     (This : in Real8) return Real8 renames
     Real8_Elementary_Functions.Sqrt;
   function SQRT
     (This : in Real10) return Real10 renames
     Real10_Elementary_Functions.Sqrt;
   function SQRT
     (This : in Mks_Type) return Mks_Type renames
     Mks_Type_Elementary_Functions.Sqrt;
   --x inline real4  CBRT(real4  arg) { return cbrtf(arg) ; }
   --x inline real8  CBRT(real8  arg) { return cbrt(arg) ; }
   --x inline real10 CBRT(real10 arg) { return cbrtl(arg) ; }
   function CBRT
     (This : in Real4) return Real4 is
     (Real4_Elementary_Functions."**" (This, 1.0 / 3.0));
   function CBRT
     (This : in Real8) return Real8 is
     (Real8_Elementary_Functions."**" (This, 1.0 / 3.0));
   function CBRT
     (This : in Real10) return Real10 is
     (Real10_Elementary_Functions."**" (This, 1.0 / 3.0));
   function CBRT
     (This : in Mks_Type) return Mks_Type is
     (Mks_Type_Elementary_Functions."**" (This, 1.0 / 3.0));
   --x inline real4  FABS(real4  arg) { return fabsf(arg) ; }
   --x inline real8  FABS(real8  arg) { return fabs(arg) ; }
   --x inline real10 FABS(real10 arg) { return fabsl(arg) ; }

   -----------------------------------------------------------------------------
   --- Indexes, ranges, and counts:
   -----------------------------------------------------------------------------

   type Domain_Index     is new Index_Type;
   type Element_Index    is new Index_Type;
   type Node_Index       is new Index_Type;
   type Region_Index     is new Index_Type;
   subtype Element_Count is Element_Index;
   subtype Node_Count    is Node_Index;
   subtype Region_Count  is Region_Index;

   subtype Process_ID_Type is MPI.Rank_Type;
   subtype Rank_Type       is MPI.Rank_Type;
   subtype Rank_Count      is MPI.Rank_Type range 1 .. MPI.Rank_Type'Last;
   subtype Process_Count   is Rank_Count;

   --- Looking down on element:
   --- 0-1-2-3 are nodes on bottom, going CCW.
   --- 4-5-6-7 are nodes directly above each.
   NODES_PER_ELEMENT     : constant := 8;
   NODES_PER_FACE        : constant := 4;
   FACES_PER_ELEMENT     : constant := 6;
   MAX_ELEMENTS_PER_NODE : constant := NODES_PER_ELEMENT;
   subtype Element_Node_Numbers is Node_Index range 0 .. NODES_PER_ELEMENT - 1;
   subtype Face_Node_Numbers    is Node_Index range 0 .. NODES_PER_FACE - 1;
   type Face_Numbers            is new Index_Type range 0 .. FACES_PER_ELEMENT - 1;
   
   -----------------------------------------------------------------------------
   --- Axes:
   -----------------------------------------------------------------------------

   --- Using this as a base type instead of having enum types so that
   --- Ada.Numerics.Generic_Real_Arrays can be used, which requres Integer
   --- array indices:
   subtype Axes_3D is Integer range 0 .. 2;

   --- Cartesian coordinates:
   --- X is positive in direction 3762 -> 0451
   --- Y is positive in direction 3047 -> 2156
   --- Z is positive in direction 3210 -> 7654
   subtype Cartesian_Axes is Axes_3D;
   X : constant Cartesian_Axes := 0;
   Y : constant Cartesian_Axes := 1;
   Z : constant Cartesian_Axes := 2;

   subtype Size_Type is Element_Index;
   type Cartesian_Size_Array is array (Cartesian_Axes) of Size_Type;

   --- Natural Coordinates (isoparametric representation):
   --- Directions are like Cartesian rotated 90 degrees CW horizontally:
   --- Xi   is positive in direction +Y: from face 0473 to 1562.
   --- Eta  is positive in direction -X: from face 0451 to 3762.
   --- Zeta is positive in direction +Z: from face 0123 to 4567.
   subtype Natural_Axes is Axes_3D;
   Xi   : constant Natural_Axes := 0;
   Eta  : constant Natural_Axes := 1;
   Zeta : constant Natural_Axes := 2;

   subtype Natural_Magnitude_Range is Real_Type range -1.0 .. 1.0;
   
   -----------------------------------------------------------------------------
   --- Arrays of indexes/ids:
   -----------------------------------------------------------------------------

   type Index_Array is array (Index_Type range <>) of Index_Type;

   type Element_Element_Index_Array is
     array (Element_Index range <>) of Element_Index;
   type Element_Element_Index_Array_Access is
     access Element_Element_Index_Array;
   procedure Release is new Ada.Unchecked_Deallocation
     (Element_Element_Index_Array,
      Element_Element_Index_Array_Access);

   type Node_Node_Index_Array is array (Node_Index range <>) of Node_Index;
   type Node_Node_Index_Array_Access is access Node_Node_Index_Array;
   procedure Release is new Ada.Unchecked_Deallocation
     (Node_Node_Index_Array,
      Node_Node_Index_Array_Access);

   type Node_Element_Count_Array is
     array (Node_Index range <>) of Element_Count;
   type Node_Element_Count_Array_Access is access Node_Element_Count_Array;
   procedure Release is new Ada.Unchecked_Deallocation
     (Node_Element_Count_Array,
      Node_Element_Count_Array_Access);

   type Faces_Per_Element_Element_Index_Array is
     array (Face_Numbers) of Element_Index;
   type Nodes_Per_Element_Node_Index_Array is
     array (Element_Node_Numbers) of Node_Index;
   
   type Nodes_Per_Face_Nodes_Per_Element_Array is
     array (Face_Node_Numbers) of Element_Node_Numbers;
   type Face_Nodes_Per_Face_Nodes_Per_Element_Array_Array is
     array (Face_Numbers) of Nodes_Per_Face_Nodes_Per_Element_Array;
   Corner_Nodes_Of_Face : constant Face_Nodes_Per_Face_Nodes_Per_Element_Array_Array :=
     (0 => (0, 1, 2, 3),
      1 => (0, 4, 5, 1),
      2 => (1, 5, 6, 2),
      3 => (2, 6, 7, 3),
      4 => (3, 7, 4, 0),
      5 => (4, 7, 6, 5));
   
   -- Not declaring a new local derived type from the vector type in this package
   -- because that introduces too many ambiguous overloaded entities:
   package Element_Count_Element_Index_Vectors is new 
     Ada.Containers.Bounded_Vectors 
       (Index_Type   => Element_Count,
        Element_Type => Element_Index);
   package ECEIV renames Element_Count_Element_Index_Vectors;

   -----------------------------------------------------------------------------
   -- Boundary conditions:
   -----------------------------------------------------------------------------

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
   type MP_Type is (M, P);
   type Boundary_Condition_Type is (Symm, Free, Common);
   type Element_Connection_Array is array (Natural_Axes, MP_Type) of Element_Index;
   type Boundary_Condition_Array is array
      (Natural_Axes, MP_Type, Boundary_Condition_Type) of Boolean with
     Pack;

   -----------------------------------------------------------------------------

   type Symmetry_Plane_Nodes_Array is array (Cartesian_Axes) of Node_Node_Index_Array_Access;
   
   -----------------------------------------------------------------------------
   --- Arrays of values/quantities:
   -----------------------------------------------------------------------------

   type Real_Array is array (Index_Type range <>) of Real_Type;
   type Real_Array_Access is access Real_Array;

   type Element_Real_Array is array (Element_Index range <>) of Real_Type;
   type Element_Real_Array_Access is access Element_Real_Array;
   procedure Release is new Ada.Unchecked_Deallocation
     (Element_Real_Array,
      Element_Real_Array_Access);

   type Region_Bin_End_Array is array (Region_Index range <>) of Cost_Type;
   type Region_Bin_End_Array_Access is access Region_Bin_End_Array;

   type Thread_Element_Count_Array is
     array (OMP.Thread_Index range <>) of Element_Count;

   -----------------------------------------------------------------------------
   -- Physical quantity types:
   -----------------------------------------------------------------------------

   --- <Name>_Type types are to avoid collisions with components, parameters, and
   --- variables named <Name>.
   --- <Name>_Relative types' objects are expected to have an initial value of 1.0.

   subtype Dimensionless    is Mks_Type;

   subtype Acceleration         is Mks_Type 
   with Dimension => (Meter => 1, Second => -2, others => 0);
   --- 1/Volume_Relative:
   subtype Compression_Relative is  Mks_Type 
   with Dimension => (Meter => -3, others => 0);
   subtype Density              is Mks_Type 
   with Dimension => (Meter => -3, Kilogram => 1, others => 0);
   subtype Energy_Type          is Energy;
   subtype Mass_Type            is MKS.Mass;
   subtype Pressure_Relative    is MKS.Pressure;
   subtype Pressure_Type        is MKS.Pressure;
   subtype Velocity             is MKS.Speed;
   -- Viscosity is like Pressure:
   subtype Viscosity            is MKS_Type
   with Dimension => (Meter => -1, Kilogram => 1, Second => -2, others => 0);
   subtype Volume_Relative      is MKS.Volume;
   subtype Volume_Type          is MKS.Volume;

   subtype Derivative_Type  is Dimensionless;
   subtype Determinant_Type is Dimensionless;
   subtype Gradient_Type    is Dimensionless;
   subtype Relative_Change  is Dimensionless;

   --- Multiply a dimensionless quantity by one of these to get a dimensioned quantity:
   M2    : constant Area         := Area (1.0);
   --     Kgpm3 : constant Density      := Density (1.0);
   Mps   : constant Velocity     := Velocity (1.0);
   --     m3    : constant Volume       := Volume (1.0);
   Mps2  : constant Acceleration := Acceleration (1.0);
   
   
   -----------------------------------------------------------------------------

   subtype Time_Span     is MKS.Time;

   Time_First      : constant Time;
   Time_Last       : constant Time;
   Time_Span_First : constant Time_Span;
   Time_Span_Last  : constant Time_Span;
   Time_Span_Zero  : constant Time_Span;
   
   function To_Time_Span
     (This : in Real_Type) return Time_Span is
     (Time_Span (This));
   
   -----------------------------------------------------------------------------
   -- Physical quantity arrays:
   -----------------------------------------------------------------------------

   --- Provides matrix and vector math:
   package Pressure_Arrays is new Ada.Numerics.Generic_Real_Arrays (Pressure);
   
   --- Dimension aspects seem to be lost when using a dimensioned type as a
   --- generic actual parameter. This means that an assignment like:
   ---      Test1 : constant Area_Vector :=
   ---          (X => Area(1.0)...
   --- gets the semantic error:
   ---         635:15 expected dimension [], found [L**2]   

   type Acceleration_Vector is array (Cartesian_Axes) of Acceleration;
   type Area_Vector         is array (Cartesian_Axes) of Area;
   type Coordinate_Vector   is array (Cartesian_Axes) of Length;
   type Force_Vector        is array (Cartesian_Axes) of Force;
   type Gradient_Vector     is array (Natural_Axes) of Gradient_Type;
   type Pressure_Vector     is new Pressure_Arrays.Real_Vector (Cartesian_Axes);
   type Strain_Vector       is array (Cartesian_Axes) of Length;
   type Velocity_Vector     is array (Cartesian_Axes) of Velocity;
   --     type Derivative_Vector is array (Cartesian_Axes) of Derivative_Type;
   
   --- Node arrays:
   
   type Node_Area_Vector_Array         is array 
     (Node_Index range <>) of Area_Vector;
   type Node_Coordinate_Array          is array 
     (Node_Index range <>) of Coordinate_Vector;
   type Node_Force_Vector_Array        is array 
     (Node_Index range <>) of Force_Vector;
   type Node_Force_Vector_Array_Access is access Node_Force_Vector_Array;
   type Node_Volume_Array              is array 
     (Node_Index range <>) of Volume;
   --     type Node_Derivative_Vector_Array       is array (Node_Index range <>) of Derivative_Vector;
   --     type Gradient_Array              is array (Node_Index range <>) of Gradient_Vector;

   subtype Nodes_Per_Element_Area_Vector_Array is
     Node_Area_Vector_Array (Element_Node_Numbers);
   subtype Nodes_Per_Element_Coordinate_Array is
     Node_Coordinate_Array (Element_Node_Numbers);
   subtype Nodes_Per_Element_Force_Vector_Array is
     Node_Force_Vector_Array (Element_Node_Numbers);
   subtype Nodes_Per_Face_Area_Vector_Array is
     Node_Area_Vector_Array (Face_Node_Numbers);
   subtype Nodes_Per_Face_Coordinate_Array is
     Node_Coordinate_Array (Face_Node_Numbers);
   --     subtype NodesPerElement_Derivative_Vector_Array is
   --       Node_Derivative_Vector_Array (NodesPerElement_Range);

   --- Element arrays:

   type Element_Compression_Array is array
     (Element_Index range <>) of Compression_Relative;
   type Element_Compression_Array_Access is access Element_Compression_Array;
   procedure Release is new Ada.Unchecked_Deallocation
     (Element_Compression_Array,
      Element_Compression_Array_Access);

   type Element_Determinant_Array is array 
     (Element_Index range <>) of Determinant_Type;
   type Element_Determinant_Array_Access is access Element_Determinant_Array;
   procedure Release is new Ada.Unchecked_Deallocation
     (Element_Determinant_Array,
      Element_Determinant_Array_Access);

   type Element_Energy_Array is array 
     (Element_Index range <>) of Energy;
   type Element_Energy_Array_Access is access Element_Energy_Array;
   procedure Release is new Ada.Unchecked_Deallocation
     (Element_Energy_Array,
      Element_Energy_Array_Access);

   type Element_Force_Vector_Array is array
     (Element_Index range <>) of Force_Vector;
   type Element_Force_Vector_Array_Access Is access Element_Force_Vector_Array;
   procedure Release is new Ada.Unchecked_Deallocation
     (Element_Force_Vector_Array,
      Element_Force_Vector_Array_Access);
   
   type Element_Pressure_Array is array 
     (Element_Index range <>) of Pressure;
   type Element_Pressure_Array_Access is access Element_Pressure_Array;
   procedure Release is new Ada.Unchecked_Deallocation
     (Element_Pressure_Array,
      Element_Pressure_Array_Access);

   type Element_Pressure_Vector_Array is array 
     (Element_Index range <>) of Pressure_Vector;
   type Element_Pressure_Vector_Array_Access is access Element_Pressure_Vector_Array;
   procedure Release is new Ada.Unchecked_Deallocation
     (Element_Pressure_Vector_Array,
      Element_Pressure_Vector_Array_Access);

   type Element_Volume_Relative_Array is array 
     (Element_Index range <>) of Volume_Relative;
   type Element_Volume_Relative_Array_Access is access Element_Volume_Relative_Array;
   procedure Release is new Ada.Unchecked_Deallocation
     (Element_Volume_Relative_Array,
      Element_Volume_Relative_Array_Access);

   ---

--     type Cartesian_Natural_Real_Array is
--       array (Cartesian_Axes, Natural_Axes) of Real_Type;
   type Cartesian_Natural_Length_Array is
     array (Cartesian_Axes, Natural_Axes) of Length;
   subtype Cofactor is MKS_Type
   with Dimension => (Meter => 2, others => 0);
   type Cartesian_Natural_Cofactor_Array is
     array (Cartesian_Axes, Natural_Axes) of Cofactor;
   --- 0..2**Cartesian_Axes'Length-1?
   --- 0..NODES_PER_ELEMENT-1?
   type Derivatives_Range is new Index_Type range 0 .. 7;
   type Derivative_Cartesian_Cofactor_Array is
     array (Derivatives_Range, Cartesian_Axes) of Cofactor;
   --     type Derivative_Array is array (Derivatives_Range) of C_Coordinate_Vector;
   
   type Thread_Time_Span_Array is array (OMP.Thread_Index range <>) of Time_Span;

   ---

   type Node_Record is record
      --x    std::vector<Real_t> m_x ;  /* coordinates */
      --x    std::vector<Real_t> m_y ;
      --x    std::vector<Real_t> m_z ;
      Coordinate : Coordinate_Vector;
      --x    std::vector<Real_t> m_xd ; /* velocities */
      --x    std::vector<Real_t> m_yd ;
      --x    std::vector<Real_t> m_zd ;
      Velocity : Velocity_Vector;
      --x    std::vector<Real_t> m_xdd ; /* accelerations */
      --x    std::vector<Real_t> m_ydd ;
      --x    std::vector<Real_t> m_zdd ;
      Acceleration : Acceleration_Vector;
      --x    std::vector<Real_t> m_fx ;  /* forces */
      --x    std::vector<Real_t> m_fy ;
      --x    std::vector<Real_t> m_fz ;
      Force : Force_Vector;
      --x    std::vector<Real_t> m_nodalMass ;  /* mass */
      Mass  : Mass_Type;
      ---    // These arrays are not used if we're not threaded
      ---    // OMP hack
      --x    Index_t *m_nodeElemStart ;
      --x    Index_t *m_nodeElemCornerList ;
      Elements_Is_Corner_Of : Element_Count_Element_Index_Vectors.Vector
        (MAX_ELEMENTS_PER_NODE);
   end record;
   type Node_Array is array (Node_Index range <>) of Node_Record;
   type Node_Array_Access is access Node_Array;
   procedure Release is new Ada.Unchecked_Deallocation
     (Node_Array,
      Node_Array_Access);

   ---    //
   ---    // Element-centered
   ---    //
   --x    Index_t*  nodelist(Index_t idx)    { return &m_nodelist[Index_t(8)*idx] ; }

   --x    Index_t nodeElemCount(Index_t idx)
   --x    { return m_nodeElemStart[idx+1] - m_nodeElemStart[idx] ; }
   --x    Index_t *nodeElemCornerList(Index_t idx)
   --x    { return &m_nodeElemCornerList[m_nodeElemStart[idx]] ; }

   type Element_Record is record
      --x    std::vector<Index_t>  m_nodelist ;     /* elemToNode connectivity */
      Corner_Nodes : Nodes_Per_Element_Node_Index_Array;
      --x    std::vector<Index_t>  m_lxim ;  /* element connectivity across each face */
      --x    std::vector<Index_t>  m_lxip ;
      --x    std::vector<Index_t>  m_letam ;
      --x    std::vector<Index_t>  m_letap ;
      --x    std::vector<Index_t>  m_lzetam ;
      --x    std::vector<Index_t>  m_lzetap ;
      Connections : Element_Connection_Array;
      ---    // elem face symm/free-surface flag
      --x    std::vector<Int_t>    m_elemBC ;  /* symmetry/free-surface flags for each elem face */
      ElemBC : Boundary_Condition_Array;
      --x    std::vector<Real_t> m_dxx ;  /* principal strains -- temporary */
      --x    std::vector<Real_t> m_dyy ;
      --x    std::vector<Real_t> m_dzz ;
      Principal_Strain : Strain_Vector;
      --x    std::vector<Real_t> m_delv_xi ;    /* velocity gradient -- temporary */
      --x    std::vector<Real_t> m_delv_eta ;
      --x    std::vector<Real_t> m_delv_zeta ;
      Velocity_Gradient : Gradient_Vector;
      ---    // Position gradient - temporary
      --x    std::vector<Real_t> m_delx_xi ;    /* coordinate gradient -- temporary */
      --x    std::vector<Real_t> m_delx_eta ;
      --x    std::vector<Real_t> m_delx_zeta ;
      Coordinate_Gradient : Gradient_Vector;
      --x    std::vector<Real_t> m_e ;   /* energy */
      Energy : Energy_Type;
      --x    std::vector<Real_t> m_p ;   /* pressure */
      --x    std::vector<Real_t> m_q ;   /* q */
      --x    std::vector<Real_t> m_ql ;  /* linear term for q */
      --x    std::vector<Real_t> m_qq ;  /* quadratic term for q */
      Pressure                       : Pressure_Type;
      Artificial_Viscosity           : Viscosity;
      Artificial_Viscosity_Linear    : Viscosity;
      Artificial_Viscosity_Quadratic : Viscosity;
      --        stress_integrated          : Force_Vector;
      --x    std::vector<Real_t> m_v ;     /* relative volume */
      --x    std::vector<Real_t> m_volo ;  /* reference volume */
      --x    std::vector<Real_t> m_vnew ;  /* new relative volume -- temporary */
      --x    std::vector<Real_t> m_delv ;  /* m_vnew - m_v */
      --x    std::vector<Real_t> m_vdov ;  /* volume derivative over volume */
      Volume           : Volume_Relative;
      Volume_Reference : Volume_Type;
      New_Volume       : Volume_Relative;
      New_Volume_Delta : Volume_Relative;
      --!! delta relative volume over (not reference) volume?  delv/volo?:
      Volume_Derivative_Over_Volume : Compression_Relative;
      --x    std::vector<Real_t> m_arealg ;  /* characteristic length of an element */
      Characteristic_Length         : Length;
      --x    std::vector<Real_t> m_ss ;      /* "sound speed" */
      Sound_Speed                   : Velocity;
      --x    std::vector<Real_t> m_elemMass ;  /* mass */
      --x    // Element mass
      Mass                          : Mass_Type;
      --x    Index_t *m_regNumList ;    // Region number per domain element
      Region : Region_Index;
   end record;
   type Element_Array is array (Element_Index range <>) of Element_Record;
   type Element_Array_Access is access Element_Array;
   procedure Release is new Ada.Unchecked_Deallocation
     (Element_Array,
      Element_Array_Access);

   type Region_Record is record
      --x    Index_t *m_regElemSize ;   // Size of region sets
      Size     : Element_Index;
      --x    Index_t **m_regElemlist ;  // region indexset
      Elements : Element_Element_Index_Array_Access;
   end record;
   type Region_Array is array (Region_Index range <>) of Region_Record;
   type Region_Array_Access is access Region_Array;
   procedure Release is new Ada.Unchecked_Deallocation
     (Region_Array,
      Region_Array_Access);

   type Row_Col_Plane is (Row, Col, Plane);
   function Succ_Wrap (This : in Row_Col_Plane) return Row_Col_Plane is
     (if This = Row_Col_Plane'Last then Row_Col_Plane'First else Row_Col_Plane'Succ (This));
   type Min_Max is (Min, Max);
   type At_Limit_Array is array (Row_Col_Plane, Min_Max) of Boolean;
   
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
      Dtcourant                         : Time_Span;
      Dthydro                           : Time_Span;
      Cycle                             : Natural;
      Dtfixed                           : Time_Span;
      Use_Courant_Condition             : Boolean;
      Current_Time                      : Time;
      Deltatime                         : Time_Span;
      Delta_Time_Multiplier_Lower_Bound : Dimensionless;
      Delta_Time_Multiplier_Upper_Bound : Dimensionless;
      Dtmax                             : Time_Span;
      Stoptime                          : Time;
      --x    Int_t   m_numRanks ;
      NumRanks                          : Rank_Count;
      --x    Index_t m_colLoc ;
      --x    Index_t m_rowLoc ;
      --x    Index_t m_planeLoc ;
      --x    Index_t m_tp ;
      -- Location : Location_Array
      ColLoc   : Domain_Index;
      RowLoc   : Domain_Index;
      PlaneLoc : Domain_Index;
      Tp       : Domain_Index;
      --x    Index_t m_maxPlaneSize ;
      --x    Index_t m_maxEdgeSize ;
      MaxPlaneSize : Node_Index;
      MaxEdgeSize  : Node_Index;
      ---    // Used in setup
      --x    Index_t m_rowMin, m_rowMax;
      --x    Index_t m_colMin, m_colMax;
      --x    Index_t m_planeMin, m_planeMax ;
      At_Limit : At_Limit_Array;
      --x    std::vector<Index_t> m_symmX ;  /* symmetry plane nodesets */
      --x    std::vector<Index_t> m_symmY ;
      --x    std::vector<Index_t> m_symmZ ;
      Symmetry_Plane_Nodes : Symmetry_Plane_Nodes_Array;
      ---    //
      ---    // MPI-Related additional data
      ---    //
      -- #if USE_MPI
      ---    // Communication Work space
      --x    Real_t *commDataSend ;
      --x    Real_t *commDataRecv ;
      CommDataSend : MPI.Comm_Buffer_Access;
      CommDataRecv : MPI.Comm_Buffer_Access;
      ---    // Maximum number of block neighbors
      --x    MPI_Request recvRequest[26] ; // 6 faces + 12 edges + 8 corners
      --x    MPI_Request sendRequest[26] ; // 6 faces + 12 edges + 8 corners
      RecvRequest : MPI_Request_Array;
      SendRequest : MPI_Request_Array;
      -- #endif
   end record;

   type Parameters_Record is record
      ---    // Parameters
      ---    // Cutoffs (treat as constants)
      ---    const Real_t  m_e_cut ;             // energy tolerance
      ---    const Real_t  m_p_cut ;             // pressure tolerance
      ---    const Real_t  m_q_cut ;             // q tolerance
      ---    const Real_t  m_v_cut ;             // relative volume tolerance
      ---    const Real_t  m_u_cut ;             // velocity tolerance
      Energy_Tolerance               : Energy;
      Pressure_Tolerance             : Pressure;
      Artificial_Viscosity_Tolerance : Viscosity;
      Volume_Relative_Tolerance      : Volume_Relative;
      Velocity_Tolerance             : Velocity;
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
      Hgcoef             : Dimensionless;
      Four_Thirds        : Dimensionless;
      Qstop              : Viscosity;
      Monoq_Max_Slope    : Dimensionless;
      Monoq_Limiter_Mult : Dimensionless;
      Qlc_Monoq          : Viscosity;
      Qqc_Monoq          : Viscosity;
      Qqc                : Viscosity;
      Eosvmax            : Volume_Relative;
      Eosvmin            : Volume_Relative;
      Pressure_Floor     : Pressure;
      Energy_Floor       : Energy;
      Dvovmax            : Compression_Relative;
      Reference_Density  : Density;
      --x    Int_t    m_cost; //imbalance cost
      Imbalance_Cost     : Cost_Type;
      --x    Index_t m_sizeX ;
      --x    Index_t m_sizeY ;
      --x    Index_t m_sizeZ ;
      Size               : Cartesian_Size_Array;
   end record;

   --x typedef Real_t &(Domain::* Domain_member )(Index_t) ;
   --x (Address of (pointer to) a Domain class function (accessor):
   ---   that takes an Index_t and returns the Real_t Domain component value)
   --- used by lulesh-comm

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

   --    // Nodes on symmertry planes
   --    bool symmXempty()          { return m_symmX.empty(); }
   --    bool symmYempty()          { return m_symmY.empty(); }
   --    bool symmZempty()          { return m_symmZ.empty(); }

   -- class Domain {
   --- Making Domain_Record tagged just to get Object.Operation notation:
   type Domain_Record is tagged private;
   type Domain_Access is access Domain_Record;

private
   package AEX renames Ada.Exceptions;
   
   Time_Span_First : constant Time_Span := Time_Span'First;
   Time_Span_Last  : constant Time_Span := Time_Span'Last;
   Time_Span_Zero  : constant Time_Span := 0.0 * S;
   Time_First      : constant Time := Time_Span_Zero;
   Time_Last       : constant Time := Time_Span_Last;
   
   --- Making Domain_Record tagged just to get Object.Operation notation:
   type Domain_Record is tagged record
      --x    Index_t m_numElem ;
      --x    Index_t m_numNode ;
      Nodes      : Node_Array_Access;
      Elements   : Element_Array_Access;
      Regions    : Region_Array_Access;
      Variables  : Variables_Record;
      Parameters : Parameters_Record;
   end record;
   
   function NumNode (This : in Domain_Record) return Node_Count is
     (if This.Nodes = null then 0 else This.Nodes.all'Length);

   function NumElem (This : in Domain_Record) return Element_Count is
     (if This.Elements = null then 0 else This.Elements.all'Length);

   function NumReg (This : in Domain_Record) return Region_Count is
     (if This.Regions = null then 0 else This.Regions.all'Length);

   --- /******************************************/

   --- /* Work Routines */

   --x static inline
   --x void TimeIncrement(Domain& domain)
   procedure TimeIncrement (Domain : in out Domain_Record) 
     with Inline;

   --x static inline
   --x void LagrangeLeapFrog(Domain& domain)
   procedure LagrangeLeapFrog (Domain : in out Domain_Record) 
     with Inline;

   --- If using MPI, does an MPI abort, otherwise raises exception X:
   procedure Abort_Or_Raise
     (X       : in AEX.Exception_Id;
      Message : in String);

end LULESH;
