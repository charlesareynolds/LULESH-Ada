package body LULESH.Par is
   --x     static inline
   --x  Real_t CalcElemVolume( const Real_t x0, const Real_t x1,
   --x                 const Real_t x2, const Real_t x3,
   --x                 const Real_t x4, const Real_t x5,
   --x                 const Real_t x6, const Real_t x7,
   --x                 const Real_t y0, const Real_t y1,
   --x                 const Real_t y2, const Real_t y3,
   --x                 const Real_t y4, const Real_t y5,
   --x                 const Real_t y6, const Real_t y7,
   --x                 const Real_t z0, const Real_t z1,
   --x                 const Real_t z2, const Real_t z3,
   --x                 const Real_t z4, const Real_t z5,
   --x                 const Real_t z6, const Real_t z7 )
   --x  {
   --x    Real_t twelveth = Real_t(1.0)/Real_t(12.0);
   --x    Real_t dx61 = x6 - x1;
   --x    Real_t dy61 = y6 - y1;
   --x    Real_t dz61 = z6 - z1;

   --x    Real_t dx70 = x7 - x0;
   --x    Real_t dy70 = y7 - y0;
   --x    Real_t dz70 = z7 - z0;

   --x    Real_t dx63 = x6 - x3;
   --x    Real_t dy63 = y6 - y3;
   --x    Real_t dz63 = z6 - z3;

   --x    Real_t dx20 = x2 - x0;
   --x    Real_t dy20 = y2 - y0;
   --x    Real_t dz20 = z2 - z0;

   --x    Real_t dx50 = x5 - x0;
   --x    Real_t dy50 = y5 - y0;
   --x    Real_t dz50 = z5 - z0;

   --x    Real_t dx64 = x6 - x4;
   --x    Real_t dy64 = y6 - y4;
   --x    Real_t dz64 = z6 - z4;

   --x    Real_t dx31 = x3 - x1;
   --x    Real_t dy31 = y3 - y1;
   --x    Real_t dz31 = z3 - z1;

   --x    Real_t dx72 = x7 - x2;
   --x    Real_t dy72 = y7 - y2;
   --x    Real_t dz72 = z7 - z2;

   --x    Real_t dx43 = x4 - x3;
   --x    Real_t dy43 = y4 - y3;
   --x    Real_t dz43 = z4 - z3;

   --x    Real_t dx57 = x5 - x7;
   --x    Real_t dy57 = y5 - y7;
   --x    Real_t dz57 = z5 - z7;

   --x    Real_t dx14 = x1 - x4;
   --x    Real_t dy14 = y1 - y4;
   --x    Real_t dz14 = z1 - z4;

   --x    Real_t dx25 = x2 - x5;
   --x    Real_t dy25 = y2 - y5;
   --x    Real_t dz25 = z2 - z5;

   twelveth : constant := 1.0/12.0;
   
   function CalcElemVolume
     (nodes : in NodesPerElement_C_Coordinate_Array)
      return Volume
   is
      --- n is just for conciseness below:
      n : constant NodesPerElement_C_Coordinate_Array := nodes;
      result : Volume := 0.0;
      --x  #define TRIPLE_PRODUCT(x1, y1, z1, x2, y2, z2, x3, y3, z3) \
      --x     ((x1)*((y2)*(z3) - (z2)*(y3)) + (x2)*((z1)*(y3) - (y1)*(z3)) + (x3)*((y1)*(z2) - (z1)*(y2)))

      function Triple_Product 
        (C1, C2, C3 : in C_Coordinate_Vector) 
         return Volume 
        with Inline is        
      begin
         return Volume
           (C1(x) * (C2(y)*C3(z) - C2(z)*C3(y)) + 
                C2(x) * (C3(y)*C1(z) - C3(z)*C1(y)) + 
                C3(x) * (C1(y)*C2(z) - C1(z)*C2(y)));
      end Triple_Product;
      
      
   begin
      --x    Real_t volume =
      --x      TRIPLE_PRODUCT(dx31 + dx72, dx63, dx20,
      --x         dy31 + dy72, dy63, dy20,
      --x         dz31 + dz72, dz63, dz20) +
      --x      TRIPLE_PRODUCT(dx43 + dx57, dx64, dx70,
      --x         dy43 + dy57, dy64, dy70,
      --x         dz43 + dz57, dz64, dz70) +
      --x      TRIPLE_PRODUCT(dx14 + dx25, dx61, dx50,
      --x         dy14 + dy25, dy61, dy50,
      --x         dz14 + dz25, dz61, dz50);

      --x  #undef TRIPLE_PRODUCT
      result :=
        Triple_Product((n(3)-n(1)) + (n(7)-n(2)), (n(6)-n(3)), (n(2)-n(0))) +
        Triple_Product((n(4)-n(3)) + (n(5)-n(7)), (n(6)-n(4)), (n(7)-n(0))) +
        Triple_Product((n(1)-n(4)) + (n(2)-n(5)), (n(6)-n(1)), (n(5)-n(0)));
      --x    volume *= twelveth;
      --x    return volume ;
      --x  }
      result := result * twelveth;
      return result ;
   end CalcElemVolume;
   
   ---  /******************************************/

   --x  //inline
   --x  Real_t CalcElemVolume( const Real_t x[8], const Real_t y[8], const Real_t z[8] )
   --x  {
   --x  return CalcElemVolume( x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7],
   --x                         y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7],
   --x                         z[0], z[1], z[2], z[3], z[4], z[5], z[6], z[7]);
   --x  }

end LULESH.Par;
