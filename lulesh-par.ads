package LULESH.Par is

   --x Real_t CalcElemVolume( const Real_t x[8],
   --x                        const Real_t y[8],
   --x                        const Real_t z[8]);
   function CalcElemVolume
     (nodes : in NodesPerElement_Coordinate_Array)
     return Volume;

end LULESH.Par;
