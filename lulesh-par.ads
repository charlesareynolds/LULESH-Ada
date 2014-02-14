package LULESH.Par is
--- // lulesh-par

   --x Real_t CalcElemVolume( const Real_t x[8],
   --x                        const Real_t y[8],
   --x                        const Real_t z[8]);
   type Eight_Reals is array (0..7) of Real_t;
   function CalcElemVolume
     (x : in Eight_Reals;
      y : in Eight_Reals;
      z : in Eight_Reals)
     return Real_t;

   --- // lulesh-par

   --x Real_t CalcElemVolume( const Real_t x[8],
   --x                        const Real_t y[8],
   --x                        const Real_t z[8]);
   type Eight_Reals is array (0..7) of Real_t;
   function CalcElemVolume
     (x : in Eight_Reals;
      y : in Eight_Reals;
      z : in Eight_Reals)
     return Real_t;


end LULESH.Par;
