with LULESH.Comm;
--- /*
---   This is a Version 2.0 MPI + OpenMP Beta implementation of LULESH

---                  Copyright (c) 2010-2013.
---       Lawrence Livermore National Security, LLC.
--- Produced at the Lawrence Livermore National Laboratory.
---                   LLNL-CODE-461231
---                 All rights reserved.

--- This file is part of LULESH, Version 2.0.
--- Please also read this link --- http://www.opensource.org/licenses/index.php

--- //////////////
--- DIFFERENCES BETWEEN THIS VERSION (2.x) AND EARLIER VERSIONS:
--- * Addition of regions to make work more representative of multi-material codes
--- * Default size of each domain is 30^3 (27000 elem) instead of 45^3. This is
---   more representative of our actual working set sizes
--- * Single source distribution supports pure serial, pure OpenMP, MPI-only,
---   and MPI+OpenMP
--- * Addition of ability to visualize the mesh using VisIt
---   https://wci.llnl.gov/codes/visit/download.html
--- * Various command line options (see ./lulesh2.0 -h)
---  -q              : quiet mode - suppress stdout
---  -i <iterations> : number of cycles to run
---  -s <size>       : length of cube mesh along side
---  -r <numregions> : Number of distinct regions (def: 11)
---  -b <balance>    : Load balance between regions of a domain (def: 1)
---  -c <cost>       : Extra cost of more expensive regions (def: 1)
---  -f <filepieces> : Number of file parts for viz output (def: np/9)
---  -p              : Print out progress
---  -v              : Output viz file (requires compiling with -DVIZ_MESH
---  -h              : This message

---  printf("Usage: %s [opts]\n", execname);
---       printf(" where [opts] is one or more of:\n");
---       printf(" -q              : quiet mode - suppress all stdout\n");
---       printf(" -i <iterations> : number of cycles to run\n");
---       printf(" -s <size>       : length of cube mesh along side\n");
---       printf(" -r <numregions> : Number of distinct regions (def: 11)\n");
---       printf(" -b <balance>    : Load balance between regions of a domain (def: 1)\n");
---       printf(" -c <cost>       : Extra cost of more expensive regions (def: 1)\n");
---       printf(" -f <numfiles>   : Number of files to split viz dump into (def: (np+10)/9)\n");
---       printf(" -p              : Print out progress\n");
---       printf(" -v              : Output viz file (requires compiling with -DVIZ_MESH\n");
---       printf(" -h              : This message\n");
---       printf("\n\n");

--- *Notable changes in LULESH 2.0

--- * Split functionality into different files
--- lulesh.cc - where most (all?) of the timed functionality lies
--- lulesh-comm.cc - MPI functionality
--- lulesh-init.cc - Setup code
--- lulesh-viz.cc  - Support for visualization option
--- lulesh-util.cc - Non-timed functions
--- *
--- * The concept of "regions" was added, although every region is the same ideal
--- *    gas material, and the same sedov blast wave problem is still the only
--- *    problem its hardcoded to solve.
--- * Regions allow two things important to making this proxy app more representative:
--- *   Four of the LULESH routines are now performed on a region-by-region basis,
--- *     making the memory access patterns non-unit stride
--- *   Artificial load imbalances can be easily introduced that could impact
--- *     parallelization strategies.
--- * The load balance flag changes region assignment.  Region number is raised to
--- *   the power entered for assignment probability.  Most likely regions changes
--- *   with MPI process id.
--- * The cost flag raises the cost of ~45% of the regions to evaluate EOS by the
--- *   entered multiple. The cost of 5% is 10x the entered multiple.
--- * MPI and OpenMP were added, and coalesced into a single version of the source
--- *   that can support serial builds, MPI-only, OpenMP-only, and MPI+OpenMP
--- * Added support to write plot files using "poor mans parallel I/O" when linked
--- *   with the silo library, which in turn can be read by VisIt.
--- * Enabled variable timestep calculation by default (courant condition), which
--- *   results in an additional reduction.
--- * Default domain (mesh) size reduced from 45^3 to 30^3
--- * Command line options to allow numerous test cases without needing to recompile
--- * Performance optimizations and code cleanup beyond LULESH 1.0
--- * Added a "Figure of Merit" calculation (elements solved per microsecond) and
--- *   output in support of using LULESH 2.0 for the 2017 CORAL procurement
--- *
--- * Possible Differences in Final Release (other changes possible)
--- *
--- * High Level mesh structure to allow data structure transformations
--- * Different default parameters
--- * Minor code performance changes and cleanup

--- TODO in future versions
--- * Add reader for (truly) unstructured meshes, probably serial only
--- * CMake based build system

--- //////////////

--- Redistribution and use in source and binary forms, with or without
--- modification, are permitted provided that the following conditions
--- are met:

---    * Redistributions of source code must retain the above copyright
---      notice, this list of conditions and the disclaimer below.

---    * Redistributions in binary form must reproduce the above copyright
---      notice, this list of conditions and the disclaimer (as noted below)
---      in the documentation and/or other materials provided with the
---      distribution.

---    * Neither the name of the LLNS/LLNL nor the names of its contributors
---      may be used to endorse or promote products derived from this software
---      without specific prior written permission.

--- THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
--- AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
--- IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
--- ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
--- THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
--- INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
--- BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
--- DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
--- OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
--- NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
--- EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


--- Additional BSD Notice

--- 1. This notice is required to be provided under our contract with the U.S.
---    Department of Energy (DOE). This work was produced at Lawrence Livermore
---    National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.

--- 2. Neither the United States Government nor Lawrence Livermore National
---    Security, LLC nor any of their employees, makes any warranty, express
---    or implied, or assumes any liability or responsibility for the accuracy,
---    completeness, or usefulness of any information, apparatus, product, or
---    process disclosed, or represents that its use would not infringe
---    privately-owned rights.

--- 3. Also, reference herein to any specific commercial products, process, or
---    services by trade name, trademark, manufacturer or otherwise does not
---    necessarily constitute or imply its endorsement, recommendation, or
---    favoring by the United States Government or Lawrence Livermore National
---    Security, LLC. The views and opinions of authors expressed herein do not
---    necessarily state or reflect those of the United States Government or
---    Lawrence Livermore National Security, LLC, and shall not be used for
---    advertising or product endorsement purposes.

--- */

-- #include <climits>
-- #include <vector>
-- #include <math.h>
-- #include <stdio.h>
-- #include <stdlib.h>
-- #include <string.h>
-- #include <ctype.h>
-- #include <time.h>
-- #include <sys/time.h>
-- #include <iostream>
-- #include <unistd.h>

-- #if _OPENMP
-- # include <omp.h>
-- #endif


package body LULESH is

   --- /*********************************/
   --- /* Data structure implementation */
   --- /*********************************/

   --- /* might want to add access methods so that memory can be */
   --- /* better managed, as in luleshFT */

   --x template <typename T>
   --x T *Allocate(size_t size)
   --x {
   --x    return static_cast<T *>(malloc(sizeof(T)*size)) ;
   --x }
--     function Allocate_Real_Array
--       (size : in Index_Type)
--        return Real_Array_Access is
--     begin
--        return new Real_Array (0..size-1);
--     end Allocate_Real_Array;

   --x template <typename T>
   --x void Release(T **ptr)
   --x {
   --x    if (*ptr != NULL) {
   --x       free(*ptr) ;
   --x       *ptr = NULL ;
   --x    }
   --x }
--     procedure Release is new Ada.Unchecked_Deallocation
--       (Object => Node_Force_Vector_Array,
--        Name   => Node_Force_Vector_Array_Access);
--
--     procedure Release is new Ada.Unchecked_Deallocation
--       (Object => Sig_Array,
--        Name   => Sig_Array_Access);

   --- /******************************************/

   --- /* Work Routines */

   --x static inline
   --x void TimeIncrement(Domain& domain)
   --x {
   ----------------------
   -- EXPORTED (private):
   ----------------------
   procedure TimeIncrement (domain : in out Domain_Record) is
      targetdt : Time_Span;
      use ART;
   begin
      --x    Real_t targetdt = domain.stoptime() - domain.time() ;
      targetdt := domain.variables.stoptime - domain.variables.current_time;

      --x    if ((domain.dtfixed() <= Real_t(0.0)) && (domain.cycle() != Int_t(0))) {
      if domain.variables.use_courant_condition and domain.variables.cycle /= 0 then
         --x       Real_t ratio ;
         --x       Real_t olddt = domain.deltatime() ;

         ---       /* This will require a reduction in parallel */
         --x       Real_t gnewdt = Real_t(1.0e+20) ;
         --x       Real_t newdt ;
         declare
            ratio  : Real_Type;
            olddt  : constant Time_Span := domain.variables.deltatime;
            gnewdt : Time_Span := Time_Span_Last;
            newdt  : Time_Span;
         begin
            --x       if (domain.dtcourant() < gnewdt) {
            --x          gnewdt = domain.dtcourant() / Real_t(2.0) ;
            --x       }
            --x       if (domain.dthydro() < gnewdt) {
            --x          gnewdt = domain.dthydro() * Real_t(2.0) / Real_t(3.0) ;
            --x       }
            if domain.variables.dtcourant < gnewdt then
               gnewdt := domain.variables.dtcourant / 2.0;
            end if;
            if domain.variables.dthydro < gnewdt then
               gnewdt := domain.variables.dthydro * 2.0 / 3.0;
            end if;

            -- #if USE_MPI
            --       MPI_Allreduce(&gnewdt, &newdt, 1,
            --                     ((sizeof(Real_t) == 4) ? MPI_FLOAT : MPI_DOUBLE),
            --                     MPI_MIN, MPI_COMM_WORLD) ;
            -- #else
            --x       newdt = gnewdt;
            -- #endif
            newdt := gnewdt;

            --x       ratio = newdt / olddt ;
            --x       if (ratio >= Real_t(1.0)) {
            --x          if (ratio < domain.deltatimemultlb()) {
            --x             newdt = olddt ;
            --x          }
            --x          else if (ratio > domain.deltatimemultub()) {
            --x             newdt = olddt*domain.deltatimemultub() ;
            --x          }
            --x       }
            ratio := Real_Type (newdt / olddt);
            if ratio >= 1.0 then
               if ratio < domain.variables.delta_time_multiplier_lower_bound then
                  newdt := olddt;
               elsif ratio > domain.variables.delta_time_multiplier_upper_bound then
                  newdt := olddt * domain.variables.delta_time_multiplier_upper_bound;
               end if;

               --x       if (newdt > domain.dtmax()) {
               --x          newdt = domain.dtmax() ;
               --x       }
               --x       domain.deltatime() = newdt ;
               if newdt > domain.variables.dtmax then
                  newdt := domain.variables.dtmax;
               end if;
               domain.variables.deltatime := newdt;
            end if;
         end;
         --x    }
      end if;

      ---    /* TRY TO PREVENT VERY SMALL SCALING ON THE NEXT CYCLE */
      --x    if ((targetdt > domain.deltatime()) &&
      --x        (targetdt < (Real_t(4.0) * domain.deltatime() / Real_t(3.0))) ) {
      --x       targetdt = Real_t(2.0) * domain.deltatime() / Real_t(3.0) ;
      --x    }
      if targetdt > domain.variables.deltatime and
        targetdt < (4.0 * domain.variables.deltatime / 3.0) then
         targetdt := 2.0 * domain.variables.deltatime / 3.0;
      end if;

      --x    if (targetdt < domain.deltatime()) {
      --x       domain.deltatime() = targetdt ;
      --x    }
      if targetdt < domain.variables.deltatime then
         domain.variables.deltatime := targetdt;
      end if;

      --x    domain.time() += domain.deltatime() ;
      domain.variables.current_time := domain.variables.current_time + domain.variables.deltatime;

      --x    ++domain.cycle() ;
      domain.variables.cycle := domain.variables.cycle + 1;
      --x }
   end TimeIncrement;

   -- /******************************************/

   --x static inline
   --x void CollectDomainNodesToElemNodes(Domain &domain,
   --x                                    const Index_t* elemToNode,
   --x                                    Real_t elemX[8],
   --x                                    Real_t elemY[8],
   --x                                    Real_t elemZ[8])
   --x {
   procedure CollectDomainNodesToElemNodes
     (domain               : in out Domain_Record;
      element_node_indexes : in     NodesPerElement_Element_Index_Array;
      element_node_coords  : out    NodesPerElement_Coordinate_Array)
     with Inline is
   begin
      for element_node in element_node_indexes'Range loop
         element_node_coords (element_node) :=
           domain.nodes (element_node_indexes (element_node)).coordinate;
      end loop;
   end CollectDomainNodesToElemNodes;

   --- /******************************************/

   --x static inline
   --x void InitStressTermsForElems(Domain &domain,
   --x                              Real_t *sigxx, Real_t *sigyy, Real_t *sigzz,
   --x                              Index_t numElem)
   --x {
   procedure InitStressTermsForElems
     (domain           : in Domain_Record;
     stress_integrated : access Element_Pressure_Vector_Array)
     with Inline
   is
      ---    //
      ---    // pull in the stresses appropriate to the hydro integration
      ---    //
   begin
      -- #pragma omp parallel for firstprivate(numElem)
      --x    for (Index_t i = 0 ; i < numElem ; ++i){
      --x       sigxx[i] = sigyy[i] = sigzz[i] =  - domain.p(i) - domain.q(i) ;
      --x    }
      for element in domain.elements'Range loop
         stress_integrated(element) :=
           (others =>
              - domain.elements(element).pressure_static
              - domain.elements(element).pressure_dynamic);
      end loop;
      --x }
   end InitStressTermsForElems;

   -- /******************************************/

   --x static inline
   --x void CalcElemShapeFunctionDerivatives( Real_t const x[],
   --x                                        Real_t const y[],
   --x                                        Real_t const z[],
   --x                                        Real_t b[][8],
   --x                                        Real_t* const volume )
   procedure CalcElemShapeFunctionDerivatives
     (enodes         : in  NodesPerElement_Coordinate_Array;
      b              : out Derivative_Cartesian_Real_Array;
      element_volume : out Volume)
     with Inline
   is
      --x {
      --x   const Real_t x0 = x[0] ;   const Real_t x1 = x[1] ;
      --x   const Real_t x2 = x[2] ;   const Real_t x3 = x[3] ;
      --x   const Real_t x4 = x[4] ;   const Real_t x5 = x[5] ;
      --x   const Real_t x6 = x[6] ;   const Real_t x7 = x[7] ;

      --x   const Real_t y0 = y[0] ;   const Real_t y1 = y[1] ;
      --x   const Real_t y2 = y[2] ;   const Real_t y3 = y[3] ;
      --x   const Real_t y4 = y[4] ;   const Real_t y5 = y[5] ;
      --x   const Real_t y6 = y[6] ;   const Real_t y7 = y[7] ;

      --x   const Real_t z0 = z[0] ;   const Real_t z1 = z[1] ;
      --x   const Real_t z2 = z[2] ;   const Real_t z3 = z[3] ;
      --x   const Real_t z4 = z[4] ;   const Real_t z5 = z[5] ;
      --x   const Real_t z6 = z[6] ;   const Real_t z7 = z[7] ;

      --x   Real_t fjxxi, fjxet, fjxze;
      --x   Real_t fjyxi, fjyet, fjyze;
      --x   Real_t fjzxi, fjzet, fjzze;
      --x   Real_t cjxxi, cjxet, cjxze;
      --x   Real_t cjyxi, cjyet, cjyze;
      --x   Real_t cjzxi, cjzet, cjzze;
      fj : Cartesian_Natural_Real_Array;
      cj : Cartesian_Natural_Real_Array;

      --x   fjxxi = Real_t(.125) * ( (x6-x0) + (x5-x3) - (x7-x1) - (x4-x2) );
      --x   fjxet = Real_t(.125) * ( (x6-x0) - (x5-x3) + (x7-x1) - (x4-x2) );
      --x   fjxze = Real_t(.125) * ( (x6-x0) + (x5-x3) + (x7-x1) + (x4-x2) );

      --x   fjyxi = Real_t(.125) * ( (y6-y0) + (y5-y3) - (y7-y1) - (y4-y2) );
      --x   fjyet = Real_t(.125) * ( (y6-y0) - (y5-y3) + (y7-y1) - (y4-y2) );
      --x   fjyze = Real_t(.125) * ( (y6-y0) + (y5-y3) + (y7-y1) + (y4-y2) );

      --x   fjzxi = Real_t(.125) * ( (z6-z0) + (z5-z3) - (z7-z1) - (z4-z2) );
      --x   fjzet = Real_t(.125) * ( (z6-z0) - (z5-z3) + (z7-z1) - (z4-z2) );
      --x   fjzze = Real_t(.125) * ( (z6-z0) + (z5-z3) + (z7-z1) + (z4-z2) );
      procedure Calculate_Naturals_For_Axis
        (axis : in Cartesian_Axes)
        with Inline
      is
         diag_6_0 : constant Real_Type := Real_Type (enodes(6)(axis) - enodes(0)(axis));
         diag_5_3 : constant Real_Type := Real_Type (enodes(5)(axis) - enodes(3)(axis));
         diag_7_1 : constant Real_Type := Real_Type (enodes(7)(axis) - enodes(1)(axis));
         diag_4_2 : constant Real_Type := Real_Type (enodes(4)(axis) - enodes(2)(axis));
         eighth   : constant := 0.125;
      begin
         fj(axis, Xi) := eighth *
           (+ diag_6_0 + diag_5_3 - diag_7_1 - diag_4_2);
         fj(axis, Eta) := eighth *
           (+ diag_6_0 - diag_5_3 + diag_7_1 - diag_4_2);
         fj(axis, Zeta) := eighth *
           (+ diag_6_0 + diag_5_3 + diag_7_1 + diag_4_2);
      end Calculate_Naturals_For_Axis;

      procedure Calculate_Naturals
        with Inline is
      begin
         for axis in Cartesian_Axes loop
            Calculate_Naturals_For_Axis (axis);
         end loop;
      end Calculate_Naturals;

      procedure Compute_Cofactors is
      begin
         ---   /* compute cofactors */
         --x   cjxxi =    (fjyet * fjzze) - (fjzet * fjyze);
         --x   cjxet =  - (fjyxi * fjzze) + (fjzxi * fjyze);
         --x   cjxze =    (fjyxi * fjzet) - (fjzxi * fjyet);
         --- Order of cjxet components reversed for visual consistency:
         cj(X,Xi)   := (fj(Y,Eta)  * fj(Z,Zeta)) - (fj(Z,Eta)  * fj(Y,Zeta));
         cj(X,Eta)  := (fj(Y,Zeta) * fj(Z,Xi))   - (fj(Z,Zeta) * fj(Y,Xi));
         cj(X,Zeta) := (fj(Y,Xi)   * fj(Z,Eta))  - (fj(Z,Xi)   * fj(Y,Eta));

         --x   cjyxi =  - (fjxet * fjzze) + (fjzet * fjxze);
         --x   cjyet =    (fjxxi * fjzze) - (fjzxi * fjxze);
         --x   cjyze =  - (fjxxi * fjzet) + (fjzxi * fjxet);
         --- Order of cjyxi and cjyze components reversed for visual consistency:
         cj(Y,Xi)   := (fj(Z,Eta)  * fj(X,Zeta)) - (fj(X,Eta)  * fj(Z,Zeta));
         cj(Y,Eta)  := (fj(Z,Zeta) * fj(X,Xi))   - (fj(X,Zeta) * fj(Z,Xi));
         cj(Y,Zeta) := (fj(Z,Xi)   * fj(X,Eta))  - (fj(X,Xi)   * fj(Z,Eta));

         --x   cjzxi =    (fjxet * fjyze) - (fjyet * fjxze);
         --x   cjzet =  - (fjxxi * fjyze) + (fjyxi * fjxze);
         --x   cjzze =    (fjxxi * fjyet) - (fjyxi * fjxet);
         --- Order of cjzet components reversed below for visual consistency:
         cj(Z,Xi)   := (fj(X,Eta)  * fj(Y,Zeta)) - (fj(Y,Eta)  * fj(X,Zeta));
         cj(Z,Eta)  := (fj(X,Zeta) * fj(Y,Xi))   - (fj(Y,Zeta) * fj(X,Xi));
         cj(Z,Zeta) := (fj(X,Xi)   * fj(Y,Eta))  - (fj(Y,Xi)   * fj(X,Eta));
      end Compute_Cofactors;

      ---   /* calculate partials :
      ---      this need only be done for l = 0,1,2,3   since , by symmetry ,
      ---      (6,7,4,5) = - (0,1,2,3) .
      ---   */
      --x   b[0][0] =   -  cjxxi  -  cjxet  -  cjxze;
      --x   b[0][1] =      cjxxi  -  cjxet  -  cjxze;
      --x   b[0][2] =      cjxxi  +  cjxet  -  cjxze;
      --x   b[0][3] =   -  cjxxi  +  cjxet  -  cjxze;
      --x   b[0][4] = -b[0][2];
      --x   b[0][5] = -b[0][3];
      --x   b[0][6] = -b[0][0];
      --x   b[0][7] = -b[0][1];

      --x   b[1][0] =   -  cjyxi  -  cjyet  -  cjyze;
      --x   b[1][1] =      cjyxi  -  cjyet  -  cjyze;
      --x   b[1][2] =      cjyxi  +  cjyet  -  cjyze;
      --x   b[1][3] =   -  cjyxi  +  cjyet  -  cjyze;
      --x   b[1][4] = -b[1][2];
      --x   b[1][5] = -b[1][3];
      --x   b[1][6] = -b[1][0];
      --x   b[1][7] = -b[1][1];

      --x   b[2][0] =   -  cjzxi  -  cjzet  -  cjzze;
      --x   b[2][1] =      cjzxi  -  cjzet  -  cjzze;
      --x   b[2][2] =      cjzxi  +  cjzet  -  cjzze;
      --x   b[2][3] =   -  cjzxi  +  cjzet  -  cjzze;
      --x   b[2][4] = -b[2][2];
      --x   b[2][5] = -b[2][3];
      --x   b[2][6] = -b[2][0];
      --x   b[2][7] = -b[2][1];
      procedure Calculate_Partial
        (axis : in Cartesian_Axes)
        with Inline is
      begin
         b(0, axis) := -  cj(axis, Xi) -  cj(axis, Eta) - cj(axis, Zeta);
         b(1, axis) := +  cj(axis, Xi) -  cj(axis, Eta) - cj(axis, Zeta);
         b(2, axis) := +  cj(axis, Xi) +  cj(axis, Eta) - cj(axis, Zeta);
         b(3, axis) := -  cj(axis, Xi) +  cj(axis, Eta) - cj(axis, Zeta);
         b(4, axis) := - b(2, axis);
         b(5, axis) := - b(3, axis);
         b(6, axis) := - b(0, axis);
         b(7, axis) := - b(1, axis);
      end Calculate_Partial;

      procedure Calculate_Partials is
      begin
         for axis in Cartesian_Axes loop
            Calculate_Partial (axis);
         end loop;
      end Calculate_Partials;

      procedure Calculate_Jacobian_Determinant is
      begin
         ---   /* calculate jacobian determinant (volume) */
         --x   *volume = Real_t(8.) * ( fjxet * cjxet + fjyet * cjyet + fjzet * cjzet);
         element_volume := 8.0 * Volume (
           (+ fj(X, Eta) * cj(X, Eta)
            + fj(Y, Eta) * cj(Y, Eta)
            + fj(Z, Eta) * cj(Z, Eta)));
         -- }
      end Calculate_Jacobian_Determinant;

   begin
      Calculate_Naturals;
      Compute_Cofactors;
      Calculate_Partials;
      Calculate_Jacobian_Determinant;
   end CalcElemShapeFunctionDerivatives;

   -- /******************************************/


   --x static inline
   --x void SumElemFaceNormal(Real_t *normalX0, Real_t *normalY0, Real_t *normalZ0,
   --x                        Real_t *normalX1, Real_t *normalY1, Real_t *normalZ1,
   --x                        Real_t *normalX2, Real_t *normalY2, Real_t *normalZ2,
   --x                        Real_t *normalX3, Real_t *normalY3, Real_t *normalZ3,
   --x                        const Real_t x0, const Real_t y0, const Real_t z0,
   --x                        const Real_t x1, const Real_t y1, const Real_t z1,
   --x                        const Real_t x2, const Real_t y2, const Real_t z2,
   --x                        const Real_t x3, const Real_t y3, const Real_t z3)
   --x {
   procedure SumElemFaceNormal
   -- Using four parms here instead of an array because this will be called with
   -- disjoint elements of an array. You can't use an aggregate in the call
   -- because it is not a variable, and these are in out parameters :
     (node0_area  : in out Area_Vector;
      node1_area  : in out Area_Vector;
      node2_area  : in out Area_Vector;
      node3_area  : in out Area_Vector;
      node_coords : in NodesPerFace_Coordinate_Array)
     with Inline
   is
      --x    Real_t bisectX0 = Real_t(0.5) * (x3 + x2 - x1 - x0);
      --x    Real_t bisectY0 = Real_t(0.5) * (y3 + y2 - y1 - y0);
      --x    Real_t bisectZ0 = Real_t(0.5) * (z3 + z2 - z1 - z0);
      --x    Real_t bisectX1 = Real_t(0.5) * (x2 + x1 - x3 - x0);
      --x    Real_t bisectY1 = Real_t(0.5) * (y2 + y1 - y3 - y0);
      --x    Real_t bisectZ1 = Real_t(0.5) * (z2 + z1 - z3 - z0);
      --x    Real_t areaX = Real_t(0.25) * (bisectY0 * bisectZ1 - bisectZ0 * bisectY1);
      --x    Real_t areaY = Real_t(0.25) * (bisectZ0 * bisectX1 - bisectX0 * bisectZ1);
      --x    Real_t areaZ = Real_t(0.25) * (bisectX0 * bisectY1 - bisectY0 * bisectX1);

      type Bisect_Range  is new Index_Type range 0..1;
      type Bisect_Coordinate_Array is array (Bisect_Range) of Coordinate_Vector;
      bisect : constant Bisect_Coordinate_Array :=
        (0 => 0.5 * (
           - node_coords(0) - node_coords(1) + node_coords(2) + node_coords(3)),
         1 => 0.5 * (
           - node_coords(0) + node_coords(1) + node_coords(2) - node_coords(3)));

      function Calc_Area
        (axis_1, axis_2 : in Cartesian_Axes)
         return Area is
        (0.25 * Area
           (bisect (0)(axis_1) * bisect (1)(axis_2) -
            bisect (0)(axis_2) * bisect (1)(axis_1)))
        with Inline;

      face_area : constant Area_Vector :=
        (X => Calc_Area (Y, Z),
         Y => Calc_Area (Z, X),
         Z => Calc_Area (X, Y));
   begin
      --x    *normalX0 += areaX;
      --x    *normalX1 += areaX;
      --x    *normalX2 += areaX;
      --x    *normalX3 += areaX;

      --x    *normalY0 += areaY;
      --x    *normalY1 += areaY;
      --x    *normalY2 += areaY;
      --x    *normalY3 += areaY;

      --x    *normalZ0 += areaZ;
      --x    *normalZ1 += areaZ;
      --x    *normalZ2 += areaZ;
      --x    *normalZ3 += areaZ;
         node0_area := node0_area + face_area;
         node1_area := node1_area + face_area;
         node2_area := node2_area + face_area;
         node3_area := node3_area + face_area;

      --x }
   end SumElemFaceNormal;

   --- /******************************************/

   --x static inline
   --x void CalcElemNodeNormals(Real_t pfx[8],
   --x                          Real_t pfy[8],
   --x                          Real_t pfz[8],
   --x                          const Real_t x[8],
   --x                          const Real_t y[8],
   --x                          const Real_t z[8])
   -- {

   procedure CalcElemNodeNormals
     (pf                  : in out NodesPerElement_Area_Vector_Array;
      element_node_coords : in     NodesPerElement_Coordinate_Array)
     with Inline is
   begin
      --x    for (Index_t i = 0 ; i < 8 ; ++i) {
      --x       pfx[i] = Real_t(0.0);
      --x       pfy[i] = Real_t(0.0);
      --x       pfz[i] = Real_t(0.0);
      --x    }
      pf := (others => (others => 0.0));

      --x    /* evaluate face one: nodes 0, 1, 2, 3 */
      --x    SumElemFaceNormal(&pfx[0], &pfy[0], &pfz[0],
      --x                   &pfx[1], &pfy[1], &pfz[1],
      --x                   &pfx[2], &pfy[2], &pfz[2],
      --x                   &pfx[3], &pfy[3], &pfz[3],
      --x                   x[0], y[0], z[0], x[1], y[1], z[1],
      --x                   x[2], y[2], z[2], x[3], y[3], z[3]);
      --x    /* evaluate face two: nodes 0, 4, 5, 1 */
      --x    SumElemFaceNormal(&pfx[0], &pfy[0], &pfz[0],
      --x                   &pfx[4], &pfy[4], &pfz[4],
      --x                   &pfx[5], &pfy[5], &pfz[5],
      --x                   &pfx[1], &pfy[1], &pfz[1],
      --x                   x[0], y[0], z[0], x[4], y[4], z[4],
      --x                   x[5], y[5], z[5], x[1], y[1], z[1]);
      --x    /* evaluate face three: nodes 1, 5, 6, 2 */
      --x    SumElemFaceNormal(&pfx[1], &pfy[1], &pfz[1],
      --x                   &pfx[5], &pfy[5], &pfz[5],
      --x                   &pfx[6], &pfy[6], &pfz[6],
      --x                   &pfx[2], &pfy[2], &pfz[2],
      --x                   x[1], y[1], z[1], x[5], y[5], z[5],
      --x                   x[6], y[6], z[6], x[2], y[2], z[2]);
      --x    /* evaluate face four: nodes 2, 6, 7, 3 */
      --x    SumElemFaceNormal(&pfx[2], &pfy[2], &pfz[2],
      --x                   &pfx[6], &pfy[6], &pfz[6],
      --x                   &pfx[7], &pfy[7], &pfz[7],
      --x                   &pfx[3], &pfy[3], &pfz[3],
      --x                   x[2], y[2], z[2], x[6], y[6], z[6],
      --x                   x[7], y[7], z[7], x[3], y[3], z[3]);
      --x    /* evaluate face five: nodes 3, 7, 4, 0 */
      --x    SumElemFaceNormal(&pfx[3], &pfy[3], &pfz[3],
      --x                   &pfx[7], &pfy[7], &pfz[7],
      --x                   &pfx[4], &pfy[4], &pfz[4],
      --x                   &pfx[0], &pfy[0], &pfz[0],
      --x                   x[3], y[3], z[3], x[7], y[7], z[7],
      --x                   x[4], y[4], z[4], x[0], y[0], z[0]);
      --x    /* evaluate face six: nodes 4, 7, 6, 5 */
      --x    SumElemFaceNormal(&pfx[4], &pfy[4], &pfz[4],
      --x                   &pfx[7], &pfy[7], &pfz[7],
      --x                   &pfx[6], &pfy[6], &pfz[6],
      --x                   &pfx[5], &pfy[5], &pfz[5],
      --x                   x[4], y[4], z[4], x[7], y[7], z[7],
      --x                   x[6], y[6], z[6], x[5], y[5], z[5]);
      for face in Face_Range loop
         SumElemFaceNormal
           (node0_area => pf (nodes_of(face)(0)),
            node1_area => pf (nodes_of(face)(1)),
            node2_area => pf (nodes_of(face)(2)),
            node3_area => pf (nodes_of(face)(3)),
            node_coords  =>
              (element_node_coords (nodes_of(face)(0)),
               element_node_coords (nodes_of(face)(1)),
               element_node_coords (nodes_of(face)(2)),
               element_node_coords (nodes_of(face)(3))));
      end loop;
   -- }
   end CalcElemNodeNormals;

   -- /******************************************/

   -- static inline
   -- void SumElemStressesToNodeForces( const Real_t B[][8],
   --                                   const Real_t stress_xx,
   --                                   const Real_t stress_yy,
   --                                   const Real_t stress_zz,
   --                                   Real_t fx[], Real_t fy[], Real_t fz[] )
   -- {
   procedure SumElemStressesToNodeForces
     (node_areas     : in     NodesPerElement_Area_Vector_Array;
      element_stress : in     Pressure_Vector;
      node_forces    : in out NodesPerElement_Force_Vector_Array)
     with Inline
   is begin
   --    for(Index_t i = 0; i < 8; i++) {
   --       fx[i] = -( stress_xx * B[0][i] );
   --       fy[i] = -( stress_yy * B[1][i] );
   --       fz[i] = -( stress_zz * B[2][i] );
   --    }
      for node in NodesPerElement_Range loop
         for axis in Cartesian_Axes loop
            node_forces(node)(axis) :=
              -(element_stress(axis) * node_areas(node)(axis));
         end loop;
      end loop;
   -- }
   end SumElemStressesToNodeForces;

   -- /******************************************/

   --x static inline
   --x void IntegrateStressForElems( Domain &domain,
   --x                               Real_t *sigxx, Real_t *sigyy, Real_t *sigzz,
   --x                               Real_t *determ, Index_t numElem, Index_t numNode)
   --x {
   procedure IntegrateStressForElems
     (domain            : in out Domain_Record;
      stress_integrated : access Element_Pressure_Vector_Array;
      determinants      : access Element_Determinant_Array)
     with Inline is

      -- #if _OPENMP
      --    Index_t numthreads = omp_get_max_threads();
      -- #else
      --t    Index_t numthreads = 1;
      numthreads : constant Index_Type := 1;
      -- #endif

      --x    Index_t numElem8 = numElem * 8 ;
      --x    Real_t *fx_elem;
      --x    Real_t *fy_elem;
      --x    Real_t *fz_elem;
      --x    Real_t fx_local[8] ;
      --x    Real_t fy_local[8] ;
      --x    Real_t fz_local[8] ;
      f_elem  : Element_NodesPerElement_Force_Vector_Array_Array_Access;
      local_nodes_forces : NodesPerElement_Force_Vector_Array;
   begin
      --x   if (numthreads > 1) {
      --x      fx_elem = Allocate<Real_t>(numElem8) ;
      --x      fy_elem = Allocate<Real_t>(numElem8) ;
      --x      fz_elem = Allocate<Real_t>(numElem8) ;
      --x   }
      if numthreads > 1 then
         f_elem := new Element_NodesPerElement_Force_Vector_Array_Array
           (0..domain.numElem-1);
      end if;

      ---   // loop over all elements

      -- #pragma omp parallel for firstprivate(numElem)
      --x   for( Index_t k=0 ; k<numElem ; ++k )
      --x   {
      for element in domain.elements'Range loop
         declare
            --x     const Index_t* const elemToNode = domain.nodelist(k);
            --x     Real_t B[3][8] ;// shape function derivatives
            --x     Real_t x_local[8] ;
            --x     Real_t y_local[8] ;
            --x     Real_t z_local[8] ;
            element_node_coords        : NodesPerElement_Coordinate_Array;
            element_node_indexes       : constant NodesPerElement_Element_Index_Array
              := domain.elements(element).node_indexes;
            B                          : NodesPerElement_Area_Vector_Array;
            shape_function_areas : Derivative_Cartesian_Real_Array;
         begin
            ---     // get nodal coordinates from global arrays and copy into local arrays.
            --x     CollectDomainNodesToElemNodes(domain, elemToNode, x_local, y_local, z_local);
            CollectDomainNodesToElemNodes
              (domain               => domain,
               element_node_indexes => element_node_indexes,
               element_node_coords  => element_node_coords);

            ---     // Volume calculation involves extra work for numerical consistency
            --x     CalcElemShapeFunctionDerivatives(x_local, y_local, z_local,
            --x                                          B, &determ[k]);
            CalcElemShapeFunctionDerivatives
              (enodes         => element_node_coords,
               B              => shape_function_areas,
               element_volume => Volume (determinants (element)));

            --- CalcElemShapeFunctionDerivatives above sets B.
            --- CalcElemNodeNormals below resets B before using it.
            --- Is B not used because of "extra work for numerical consistency"?

            --     CalcElemNodeNormals( B[0] , B[1], B[2],
            --                           x_local, y_local, z_local );
            CalcElemNodeNormals
              (pf                  => B,
               element_node_coords => element_node_coords);

            --     if (numthreads > 1) {
            if (numthreads > 1) then
               --        // Eliminate thread writing conflicts at the nodes by giving
               --        // each element its own copy to write to
               --        SumElemStressesToNodeForces( B, sigxx[k], sigyy[k], sigzz[k],
               --                                     &fx_elem[k*8],
               --                                     &fy_elem[k*8],
               --                                     &fz_elem[k*8] ) ;
               SumElemStressesToNodeForces
                 (node_areas     => B,
                  element_stress => stress_integrated (element),
                  node_forces    => f_elem (element));
               --     }
               --     else {
            else
               -- x       SumElemStressesToNodeForces( B, sigxx[k], sigyy[k], sigzz[k],
               -- x                                    fx_local, fy_local, fz_local ) ;
               SumElemStressesToNodeForces
                 (node_areas     => B,
                  element_stress => stress_integrated (element),
                  node_forces    => local_nodes_forces);

               --        // copy nodal force contributions to global force arrray.
               --        for( Index_t lnode=0 ; lnode<8 ; ++lnode ) {
               --           Index_t gnode = elemToNode[lnode];
               --           domain.fx(gnode) += fx_local[lnode];
               --           domain.fy(gnode) += fy_local[lnode];
               --           domain.fz(gnode) += fz_local[lnode];
               --        }
               --     }
               for element_node in NodesPerElement_Range loop
                  declare
                     node : constant Node_Index :=
                       element_node_indexes (element_node);
                  begin
                     domain.nodes(node).force := domain.nodes(node).force +
                       local_nodes_forces(element_node);
                  end;
               end loop;
            end if;
         end;
         --   }
      end loop;

      --   if (numthreads > 1) {
      --      // If threaded, then we need to copy the data out of the temporary
      --      // arrays used above into the final forces field
      -- #pragma omp parallel for firstprivate(numNode)
      --      for( Index_t gnode=0 ; gnode<numNode ; ++gnode )
      --      {
      --         Index_t count = domain.nodeElemCount(gnode) ;
      --         Index_t *cornerList = domain.nodeElemCornerList(gnode) ;
      --         Real_t fx_tmp = Real_t(0.0) ;
      --         Real_t fy_tmp = Real_t(0.0) ;
      --         Real_t fz_tmp = Real_t(0.0) ;
      --         for (Index_t i=0 ; i < count ; ++i) {
      --            Index_t elem = cornerList[i] ;
      --            fx_tmp += fx_elem[elem] ;
      --            fy_tmp += fy_elem[elem] ;
      --            fz_tmp += fz_elem[elem] ;
      --         }
      --         domain.fx(gnode) = fx_tmp ;
      --         domain.fy(gnode) = fy_tmp ;
      --         domain.fz(gnode) = fz_tmp ;
      --      }
      --      Release(&fz_elem) ;
      --      Release(&fy_elem) ;
      --      Release(&fx_elem) ;
      --   }
      -- }
   end IntegrateStressForElems;

   -- /******************************************/

   twelfth : constant := 1.0 / 12.0;

   -- static inline
   -- void VoluDer(const Real_t x0, const Real_t x1, const Real_t x2,
   --              const Real_t x3, const Real_t x4, const Real_t x5,
   --              const Real_t y0, const Real_t y1, const Real_t y2,
   --              const Real_t y3, const Real_t y4, const Real_t y5,
   --              const Real_t z0, const Real_t z1, const Real_t z2,
   --              const Real_t z3, const Real_t z4, const Real_t z5,
   --              Real_t* dvdx, Real_t* dvdy, Real_t* dvdz)
   -- {
   procedure VoluDer
     (fc  : in  NodesPerElement_XYZ_Real_Array_Array;
      dvd : out XYZ_Real_Array)
     with Inline
   is
   begin
   --    const Real_t twelfth = Real_t(1.0) / Real_t(12.0) ;
   --    *dvdx =
   --       (y1 + y2) * (z0 + z1) - (y0 + y1) * (z1 + z2) +
   --       (y0 + y4) * (z3 + z4) - (y3 + y4) * (z0 + z4) -
   --       (y2 + y5) * (z3 + z5) + (y3 + y5) * (z2 + z5);
   --    *dvdy =
   --       - (x1 + x2) * (z0 + z1) + (x0 + x1) * (z1 + z2) -
   --       (x0 + x4) * (z3 + z4) + (x3 + x4) * (z0 + z4) +
   --       (x2 + x5) * (z3 + z5) - (x3 + x5) * (z2 + z5);
   --    *dvdz =
   --       - (y1 + y2) * (x0 + x1) + (y0 + y1) * (x1 + x2) -
   --       (y0 + y4) * (x3 + x4) + (y3 + y4) * (x0 + x4) +
   --       (y2 + y5) * (x3 + x5) - (y3 + y5) * (x2 + x5);

   --    *dvdx *= twelfth;
   --    *dvdy *= twelfth;
   --    *dvdz *= twelfth;
      pragma Warnings (Off, "unary minus expression should be parenthesized here");
      dvd(X) := twelfth *
        (+ (fc(1)(y) + fc(2)(y)) * (fc(0)(z) + fc(1)(z))
         - (fc(0)(y) + fc(1)(y)) * (fc(1)(z) + fc(2)(z))
         + (fc(0)(y) + fc(4)(y)) * (fc(3)(z) + fc(4)(z))
         - (fc(3)(y) + fc(4)(y)) * (fc(0)(z) + fc(4)(z))
         - (fc(2)(y) + fc(5)(y)) * (fc(3)(z) + fc(5)(z))
         + (fc(3)(y) + fc(5)(y)) * (fc(2)(z) + fc(5)(z)));
      dvd(Y) := twelfth *
        (- (fc(1)(x) + fc(2)(x)) * (fc(0)(z) + fc(1)(z))
         + (fc(0)(x) + fc(1)(x)) * (fc(1)(z) + fc(2)(z))
         - (fc(0)(x) + fc(4)(x)) * (fc(3)(z) + fc(4)(z))
         + (fc(3)(x) + fc(4)(x)) * (fc(0)(z) + fc(4)(z))
         + (fc(2)(x) + fc(5)(x)) * (fc(3)(z) + fc(5)(z))
         - (fc(3)(x) + fc(5)(x)) * (fc(2)(z) + fc(5)(z)));
      dvd(Z) := twelfth *
        (- (fc(1)(y) + fc(2)(y)) * (fc(0)(x) + fc(1)(x))
         + (fc(0)(y) + fc(1)(y)) * (fc(1)(x) + fc(2)(x))
         - (fc(0)(y) + fc(4)(y)) * (fc(3)(x) + fc(4)(x))
         + (fc(3)(y) + fc(4)(y)) * (fc(0)(x) + fc(4)(x))
         + (fc(2)(y) + fc(5)(y)) * (fc(3)(x) + fc(5)(x))
         - (fc(3)(y) + fc(5)(y)) * (fc(2)(x) + fc(5)(x)));
      pragma Warnings (On, "unary minus expression should be parenthesized here");
   -- }
   end VoluDer;

   -- /******************************************/

   -- static inline
   -- void CalcElemVolumeDerivative(Real_t dvdx[8],
   --                               Real_t dvdy[8],
   --                               Real_t dvdz[8],
   --                               const Real_t x[8],
   --                               const Real_t y[8],
   --                               const Real_t z[8])
   -- {
   --    VoluDer(x[1], x[2], x[3], x[4], x[5], x[7],
   --            y[1], y[2], y[3], y[4], y[5], y[7],
   --            z[1], z[2], z[3], z[4], z[5], z[7],
   --            &dvdx[0], &dvdy[0], &dvdz[0]);
   --    VoluDer(x[0], x[1], x[2], x[7], x[4], x[6],
   --            y[0], y[1], y[2], y[7], y[4], y[6],
   --            z[0], z[1], z[2], z[7], z[4], z[6],
   --            &dvdx[3], &dvdy[3], &dvdz[3]);
   --    VoluDer(x[3], x[0], x[1], x[6], x[7], x[5],
   --            y[3], y[0], y[1], y[6], y[7], y[5],
   --            z[3], z[0], z[1], z[6], z[7], z[5],
   --            &dvdx[2], &dvdy[2], &dvdz[2]);
   --    VoluDer(x[2], x[3], x[0], x[5], x[6], x[4],
   --            y[2], y[3], y[0], y[5], y[6], y[4],
   --            z[2], z[3], z[0], z[5], z[6], z[4],
   --            &dvdx[1], &dvdy[1], &dvdz[1]);
   --    VoluDer(x[7], x[6], x[5], x[0], x[3], x[1],
   --            y[7], y[6], y[5], y[0], y[3], y[1],
   --            z[7], z[6], z[5], z[0], z[3], z[1],
   --            &dvdx[4], &dvdy[4], &dvdz[4]);
   --    VoluDer(x[4], x[7], x[6], x[1], x[0], x[2],
   --            y[4], y[7], y[6], y[1], y[0], y[2],
   --            z[4], z[7], z[6], z[1], z[0], z[2],
   --            &dvdx[5], &dvdy[5], &dvdz[5]);
   --    VoluDer(x[5], x[4], x[7], x[2], x[1], x[3],
   --            y[5], y[4], y[7], y[2], y[1], y[3],
   --            z[5], z[4], z[7], z[2], z[1], z[3],
   --            &dvdx[6], &dvdy[6], &dvdz[6]);
   --    VoluDer(x[6], x[5], x[4], x[3], x[2], x[0],
   --            y[6], y[5], y[4], y[3], y[2], y[0],
   --            z[6], z[5], z[4], z[3], z[2], z[0],
   --            &dvdx[7], &dvdy[7], &dvdz[7]);
   -- }

   -- /******************************************/

   -- static inline
   -- void CalcElemFBHourglassForce(Real_t *xd, Real_t *yd, Real_t *zd,  Real_t hourgam[][4],
   --                               Real_t coefficient,
   --                               Real_t *hgfx, Real_t *hgfy, Real_t *hgfz )
   -- {
   --    Real_t hxx[4];
   --    for(Index_t i = 0; i < 4; i++) {
   --       hxx[i] = hourgam[0][i] * xd[0] + hourgam[1][i] * xd[1] +
   --                hourgam[2][i] * xd[2] + hourgam[3][i] * xd[3] +
   --                hourgam[4][i] * xd[4] + hourgam[5][i] * xd[5] +
   --                hourgam[6][i] * xd[6] + hourgam[7][i] * xd[7];
   --    }
   --    for(Index_t i = 0; i < 8; i++) {
   --       hgfx[i] = coefficient *
   --                 (hourgam[i][0] * hxx[0] + hourgam[i][1] * hxx[1] +
   --                  hourgam[i][2] * hxx[2] + hourgam[i][3] * hxx[3]);
   --    }
   --    for(Index_t i = 0; i < 4; i++) {
   --       hxx[i] = hourgam[0][i] * yd[0] + hourgam[1][i] * yd[1] +
   --                hourgam[2][i] * yd[2] + hourgam[3][i] * yd[3] +
   --                hourgam[4][i] * yd[4] + hourgam[5][i] * yd[5] +
   --                hourgam[6][i] * yd[6] + hourgam[7][i] * yd[7];
   --    }
   --    for(Index_t i = 0; i < 8; i++) {
   --       hgfy[i] = coefficient *
   --                 (hourgam[i][0] * hxx[0] + hourgam[i][1] * hxx[1] +
   --                  hourgam[i][2] * hxx[2] + hourgam[i][3] * hxx[3]);
   --    }
   --    for(Index_t i = 0; i < 4; i++) {
   --       hxx[i] = hourgam[0][i] * zd[0] + hourgam[1][i] * zd[1] +
   --                hourgam[2][i] * zd[2] + hourgam[3][i] * zd[3] +
   --                hourgam[4][i] * zd[4] + hourgam[5][i] * zd[5] +
   --                hourgam[6][i] * zd[6] + hourgam[7][i] * zd[7];
   --    }
   --    for(Index_t i = 0; i < 8; i++) {
   --       hgfz[i] = coefficient *
   --                 (hourgam[i][0] * hxx[0] + hourgam[i][1] * hxx[1] +
   --                  hourgam[i][2] * hxx[2] + hourgam[i][3] * hxx[3]);
   --    }
   -- }

   -- /******************************************/

   -- static inline
   -- void CalcFBHourglassForceForElems( Domain &domain,
   --                                    Real_t *determ,
   --                                    Real_t *x8n, Real_t *y8n, Real_t *z8n,
   --                                    Real_t *dvdx, Real_t *dvdy, Real_t *dvdz,
   --                                    Real_t hourg, Index_t numElem,
   --                                    Index_t numNode)
   -- {

   -- #if _OPENMP
   --    Index_t numthreads = omp_get_max_threads();
   -- #else
   --    Index_t numthreads = 1;
   -- #endif
   --    /*************************************************
   --     *
   --     *     FUNCTION: Calculates the Flanagan-Belytschko anti-hourglass
   --     *               force.
   --     *
   --     *************************************************/

   --    Index_t numElem8 = numElem * 8 ;

   --    Real_t *fx_elem;
   --    Real_t *fy_elem;
   --    Real_t *fz_elem;

   --    if(numthreads > 1) {
   --       fx_elem = Allocate<Real_t>(numElem8) ;
   --       fy_elem = Allocate<Real_t>(numElem8) ;
   --       fz_elem = Allocate<Real_t>(numElem8) ;
   --    }

   --    Real_t  gamma[4][8];

   --    gamma[0][0] = Real_t( 1.);
   --    gamma[0][1] = Real_t( 1.);
   --    gamma[0][2] = Real_t(-1.);
   --    gamma[0][3] = Real_t(-1.);
   --    gamma[0][4] = Real_t(-1.);
   --    gamma[0][5] = Real_t(-1.);
   --    gamma[0][6] = Real_t( 1.);
   --    gamma[0][7] = Real_t( 1.);

   --    gamma[1][0] = Real_t( 1.);
   --    gamma[1][1] = Real_t(-1.);
   --    gamma[1][2] = Real_t(-1.);
   --    gamma[1][3] = Real_t( 1.);
   --    gamma[1][4] = Real_t(-1.);
   --    gamma[1][5] = Real_t( 1.);
   --    gamma[1][6] = Real_t( 1.);
   --    gamma[1][7] = Real_t(-1.);

   --    gamma[2][0] = Real_t( 1.);
   --    gamma[2][1] = Real_t(-1.);
   --    gamma[2][2] = Real_t( 1.);
   --    gamma[2][3] = Real_t(-1.);
   --    gamma[2][4] = Real_t( 1.);
   --    gamma[2][5] = Real_t(-1.);
   --    gamma[2][6] = Real_t( 1.);
   --    gamma[2][7] = Real_t(-1.);

   --    gamma[3][0] = Real_t(-1.);
   --    gamma[3][1] = Real_t( 1.);
   --    gamma[3][2] = Real_t(-1.);
   --    gamma[3][3] = Real_t( 1.);
   --    gamma[3][4] = Real_t( 1.);
   --    gamma[3][5] = Real_t(-1.);
   --    gamma[3][6] = Real_t( 1.);
   --    gamma[3][7] = Real_t(-1.);

   -- /*************************************************/
   -- /*    compute the hourglass modes */


   -- #pragma omp parallel for firstprivate(numElem, hourg)
   --    for(Index_t i2=0;i2<numElem;++i2){
   --       Real_t *fx_local, *fy_local, *fz_local ;
   --       Real_t hgfx[8], hgfy[8], hgfz[8] ;

   --       Real_t coefficient;

   --       Real_t hourgam[8][4];
   --       Real_t xd1[8], yd1[8], zd1[8] ;

   --       const Index_t *elemToNode = domain.nodelist(i2);
   --       Index_t i3=8*i2;
   --       Real_t volinv=Real_t(1.0)/determ[i2];
   --       Real_t ss1, mass1, volume13 ;
   --       for(Index_t i1=0;i1<4;++i1){

   --          Real_t hourmodx =
   --             x8n[i3] * gamma[i1][0] + x8n[i3+1] * gamma[i1][1] +
   --             x8n[i3+2] * gamma[i1][2] + x8n[i3+3] * gamma[i1][3] +
   --             x8n[i3+4] * gamma[i1][4] + x8n[i3+5] * gamma[i1][5] +
   --             x8n[i3+6] * gamma[i1][6] + x8n[i3+7] * gamma[i1][7];

   --          Real_t hourmody =
   --             y8n[i3] * gamma[i1][0] + y8n[i3+1] * gamma[i1][1] +
   --             y8n[i3+2] * gamma[i1][2] + y8n[i3+3] * gamma[i1][3] +
   --             y8n[i3+4] * gamma[i1][4] + y8n[i3+5] * gamma[i1][5] +
   --             y8n[i3+6] * gamma[i1][6] + y8n[i3+7] * gamma[i1][7];

   --          Real_t hourmodz =
   --             z8n[i3] * gamma[i1][0] + z8n[i3+1] * gamma[i1][1] +
   --             z8n[i3+2] * gamma[i1][2] + z8n[i3+3] * gamma[i1][3] +
   --             z8n[i3+4] * gamma[i1][4] + z8n[i3+5] * gamma[i1][5] +
   --             z8n[i3+6] * gamma[i1][6] + z8n[i3+7] * gamma[i1][7];

   --          hourgam[0][i1] = gamma[i1][0] -  volinv*(dvdx[i3  ] * hourmodx +
   --                                                   dvdy[i3  ] * hourmody +
   --                                                   dvdz[i3  ] * hourmodz );

   --          hourgam[1][i1] = gamma[i1][1] -  volinv*(dvdx[i3+1] * hourmodx +
   --                                                   dvdy[i3+1] * hourmody +
   --                                                   dvdz[i3+1] * hourmodz );

   --          hourgam[2][i1] = gamma[i1][2] -  volinv*(dvdx[i3+2] * hourmodx +
   --                                                   dvdy[i3+2] * hourmody +
   --                                                   dvdz[i3+2] * hourmodz );

   --          hourgam[3][i1] = gamma[i1][3] -  volinv*(dvdx[i3+3] * hourmodx +
   --                                                   dvdy[i3+3] * hourmody +
   --                                                   dvdz[i3+3] * hourmodz );

   --          hourgam[4][i1] = gamma[i1][4] -  volinv*(dvdx[i3+4] * hourmodx +
   --                                                   dvdy[i3+4] * hourmody +
   --                                                   dvdz[i3+4] * hourmodz );

   --          hourgam[5][i1] = gamma[i1][5] -  volinv*(dvdx[i3+5] * hourmodx +
   --                                                   dvdy[i3+5] * hourmody +
   --                                                   dvdz[i3+5] * hourmodz );

   --          hourgam[6][i1] = gamma[i1][6] -  volinv*(dvdx[i3+6] * hourmodx +
   --                                                   dvdy[i3+6] * hourmody +
   --                                                   dvdz[i3+6] * hourmodz );

   --          hourgam[7][i1] = gamma[i1][7] -  volinv*(dvdx[i3+7] * hourmodx +
   --                                                   dvdy[i3+7] * hourmody +
   --                                                   dvdz[i3+7] * hourmodz );

   --       }

   --       /* compute forces */
   --       /* store forces into h arrays (force arrays) */

   --       ss1=domain.ss(i2);
   --       mass1=domain.elemMass(i2);
   --       volume13=CBRT(determ[i2]);

   --       Index_t n0si2 = elemToNode[0];
   --       Index_t n1si2 = elemToNode[1];
   --       Index_t n2si2 = elemToNode[2];
   --       Index_t n3si2 = elemToNode[3];
   --       Index_t n4si2 = elemToNode[4];
   --       Index_t n5si2 = elemToNode[5];
   --       Index_t n6si2 = elemToNode[6];
   --       Index_t n7si2 = elemToNode[7];

   --       xd1[0] = domain.xd(n0si2);
   --       xd1[1] = domain.xd(n1si2);
   --       xd1[2] = domain.xd(n2si2);
   --       xd1[3] = domain.xd(n3si2);
   --       xd1[4] = domain.xd(n4si2);
   --       xd1[5] = domain.xd(n5si2);
   --       xd1[6] = domain.xd(n6si2);
   --       xd1[7] = domain.xd(n7si2);

   --       yd1[0] = domain.yd(n0si2);
   --       yd1[1] = domain.yd(n1si2);
   --       yd1[2] = domain.yd(n2si2);
   --       yd1[3] = domain.yd(n3si2);
   --       yd1[4] = domain.yd(n4si2);
   --       yd1[5] = domain.yd(n5si2);
   --       yd1[6] = domain.yd(n6si2);
   --       yd1[7] = domain.yd(n7si2);

   --       zd1[0] = domain.zd(n0si2);
   --       zd1[1] = domain.zd(n1si2);
   --       zd1[2] = domain.zd(n2si2);
   --       zd1[3] = domain.zd(n3si2);
   --       zd1[4] = domain.zd(n4si2);
   --       zd1[5] = domain.zd(n5si2);
   --       zd1[6] = domain.zd(n6si2);
   --       zd1[7] = domain.zd(n7si2);

   --       coefficient = - hourg * Real_t(0.01) * ss1 * mass1 / volume13;

   --       CalcElemFBHourglassForce(xd1,yd1,zd1,
   --                       hourgam,
   --                       coefficient, hgfx, hgfy, hgfz);

   --       // With the threaded version, we write into local arrays per elem
   --       // so we don't have to worry about race conditions
   --       if (numthreads > 1) {
   --          fx_local = &fx_elem[i3] ;
   --          fx_local[0] = hgfx[0];
   --          fx_local[1] = hgfx[1];
   --          fx_local[2] = hgfx[2];
   --          fx_local[3] = hgfx[3];
   --          fx_local[4] = hgfx[4];
   --          fx_local[5] = hgfx[5];
   --          fx_local[6] = hgfx[6];
   --          fx_local[7] = hgfx[7];

   --          fy_local = &fy_elem[i3] ;
   --          fy_local[0] = hgfy[0];
   --          fy_local[1] = hgfy[1];
   --          fy_local[2] = hgfy[2];
   --          fy_local[3] = hgfy[3];
   --          fy_local[4] = hgfy[4];
   --          fy_local[5] = hgfy[5];
   --          fy_local[6] = hgfy[6];
   --          fy_local[7] = hgfy[7];

   --          fz_local = &fz_elem[i3] ;
   --          fz_local[0] = hgfz[0];
   --          fz_local[1] = hgfz[1];
   --          fz_local[2] = hgfz[2];
   --          fz_local[3] = hgfz[3];
   --          fz_local[4] = hgfz[4];
   --          fz_local[5] = hgfz[5];
   --          fz_local[6] = hgfz[6];
   --          fz_local[7] = hgfz[7];
   --       }
   --       else {
   --          domain.fx(n0si2) += hgfx[0];
   --          domain.fy(n0si2) += hgfy[0];
   --          domain.fz(n0si2) += hgfz[0];

   --          domain.fx(n1si2) += hgfx[1];
   --          domain.fy(n1si2) += hgfy[1];
   --          domain.fz(n1si2) += hgfz[1];

   --          domain.fx(n2si2) += hgfx[2];
   --          domain.fy(n2si2) += hgfy[2];
   --          domain.fz(n2si2) += hgfz[2];

   --          domain.fx(n3si2) += hgfx[3];
   --          domain.fy(n3si2) += hgfy[3];
   --          domain.fz(n3si2) += hgfz[3];

   --          domain.fx(n4si2) += hgfx[4];
   --          domain.fy(n4si2) += hgfy[4];
   --          domain.fz(n4si2) += hgfz[4];

   --          domain.fx(n5si2) += hgfx[5];
   --          domain.fy(n5si2) += hgfy[5];
   --          domain.fz(n5si2) += hgfz[5];

   --          domain.fx(n6si2) += hgfx[6];
   --          domain.fy(n6si2) += hgfy[6];
   --          domain.fz(n6si2) += hgfz[6];

   --          domain.fx(n7si2) += hgfx[7];
   --          domain.fy(n7si2) += hgfy[7];
   --          domain.fz(n7si2) += hgfz[7];
   --       }
   --    }

   --    if (numthreads > 1) {
   --      // Collect the data from the local arrays into the final force arrays
   -- #pragma omp parallel for firstprivate(numNode)
   --       for( Index_t gnode=0 ; gnode<numNode ; ++gnode )
   --       {
   --          Index_t count = domain.nodeElemCount(gnode) ;
   --          Index_t *cornerList = domain.nodeElemCornerList(gnode) ;
   --          Real_t fx_tmp = Real_t(0.0) ;
   --          Real_t fy_tmp = Real_t(0.0) ;
   --          Real_t fz_tmp = Real_t(0.0) ;
   --          for (Index_t i=0 ; i < count ; ++i) {
   --             Index_t elem = cornerList[i] ;
   --             fx_tmp += fx_elem[elem] ;
   --             fy_tmp += fy_elem[elem] ;
   --             fz_tmp += fz_elem[elem] ;
   --          }
   --          domain.fx(gnode) += fx_tmp ;
   --          domain.fy(gnode) += fy_tmp ;
   --          domain.fz(gnode) += fz_tmp ;
   --       }
   --       Release(&fz_elem) ;
   --       Release(&fy_elem) ;
   --       Release(&fx_elem) ;
   --    }
   -- }

   -- /******************************************/

   -- static inline
   -- void CalcHourglassControlForElems(Domain& domain,
   --                                   Real_t determ[], Real_t hgcoef)
   -- {
   --    Index_t numElem = domain.numElem() ;
   --    Index_t numElem8 = numElem * 8 ;
   --    Real_t *dvdx = Allocate<Real_t>(numElem8) ;
   --    Real_t *dvdy = Allocate<Real_t>(numElem8) ;
   --    Real_t *dvdz = Allocate<Real_t>(numElem8) ;
   --    Real_t *x8n  = Allocate<Real_t>(numElem8) ;
   --    Real_t *y8n  = Allocate<Real_t>(numElem8) ;
   --    Real_t *z8n  = Allocate<Real_t>(numElem8) ;

   --    /* start loop over elements */
   -- #pragma omp parallel for firstprivate(numElem)
   --    for (Index_t i=0 ; i<numElem ; ++i){
   --       Real_t  x1[8],  y1[8],  z1[8] ;
   --       Real_t pfx[8], pfy[8], pfz[8] ;

   --       Index_t* elemToNode = domain.nodelist(i);
   --       CollectDomainNodesToElemNodes(domain, elemToNode, x1, y1, z1);

   --       CalcElemVolumeDerivative(pfx, pfy, pfz, x1, y1, z1);

   --       /* load into temporary storage for FB Hour Glass control */
   --       for(Index_t ii=0;ii<8;++ii){
   --          Index_t jj=8*i+ii;

   --          dvdx[jj] = pfx[ii];
   --          dvdy[jj] = pfy[ii];
   --          dvdz[jj] = pfz[ii];
   --          x8n[jj]  = x1[ii];
   --          y8n[jj]  = y1[ii];
   --          z8n[jj]  = z1[ii];
   --       }

   --       determ[i] = domain.volo(i) * domain.v(i);

   --       /* Do a check for negative volumes */
   --       if ( domain.v(i) <= Real_t(0.0) ) {
   -- #if USE_MPI
   --          MPI_Abort(MPI_COMM_WORLD, VolumeError) ;
   -- #else
   --          exit(VolumeError);
   -- #endif
   --       }
   --    }

   --    if ( hgcoef > Real_t(0.) ) {
   --       CalcFBHourglassForceForElems( domain,
   --                                     determ, x8n, y8n, z8n, dvdx, dvdy, dvdz,
   --                                     hgcoef, numElem, domain.numNode()) ;
   --    }

   --    Release(&z8n) ;
   --    Release(&y8n) ;
   --    Release(&x8n) ;
   --    Release(&dvdz) ;
   --    Release(&dvdy) ;
   --    Release(&dvdx) ;

   --    return ;
   -- }

   -- /******************************************/

   --x static inline
   --x void CalcVolumeForceForElems(Domain& domain)
   --x {
   procedure CalcVolumeForceForElems
     (domain : in out Domain_Record)
   is
      --x    Index_t numElem = domain.numElem() ;
      numElem : constant Element_Index := domain.numElem;
   begin
      --    if (numElem != 0) {
      if numElem > 0 then
         declare
            --       Real_t  hgcoef = domain.hgcoef() ;
            --       Real_t *sigxx  = Allocate<Real_t>(numElem) ;
            --       Real_t *sigyy  = Allocate<Real_t>(numElem) ;
            --       Real_t *sigzz  = Allocate<Real_t>(numElem) ;
            --       Real_t *determ = Allocate<Real_t>(numElem) ;
            hgcoef            : Real_Type := domain.parameters.hgcoef;
            stress_integrated : Element_Pressure_Vector_Array_Access  :=
              new Element_Pressure_Vector_Array (0..numElem - 1);
            determinants            : Element_Determinant_Array_Access :=
              new Element_Determinant_Array (0..numElem - 1);
         begin
            ---       /* Sum contributions to total stress tensor */
            --x       InitStressTermsForElems(domain, sigxx, sigyy, sigzz, numElem);
            ---       // call elemlib stress integration loop to produce nodal forces from
            ---       // material stresses.
            --x       IntegrateStressForElems( domain,
            --x                                sigxx, sigyy, sigzz, determ, numElem,
            --x                                domain.numNode()) ;
            InitStressTermsForElems(domain, stress_integrated);
            IntegrateStressForElems(domain, stress_integrated, determinants);

            ---       // check for negative element volume
            -- #pragma omp parallel for firstprivate(numElem)
            --x       for ( Index_t k=0 ; k<numElem ; ++k ) {
            --x          if (determ[k] <= Real_t(0.0)) {
            -- #if USE_MPI
            --             MPI_Abort(MPI_COMM_WORLD, VolumeError) ;
            -- #else
            --x             exit(VolumeError);
            -- #endif
            --          }
            --       }
            for element in determinants'Range loop
               if determinants(element) <= 0.0 then
                  raise VolumeError;
               end if;
            end loop;

            --x       CalcHourglassControlForElems(domain, determ, hgcoef) ;
            CalcHourglassControlForElems(domain, determinants, hgcoef);

            --x       Release(&determ) ;
            --x       Release(&sigzz) ;
            --x       Release(&sigyy) ;
            --x       Release(&sigxx) ;
            Release (determinants);
            Release (stress_integrated);
         end;
         --x    }
      end if;
      --x }
   end CalcVolumeForceForElems;

   -- /******************************************/

   --x static inline void CalcForceForNodes(Domain& domain)
   --x {
   procedure CalcForceForNodes(domain : in out Domain_Record) is
      --   Index_t numNode = domain.numNode() ;
      numNode : Node_Index := domain.numNode;
   begin

      -- #if USE_MPI
      --   CommRecv(domain, MSG_COMM_SBN, 3,
      --            domain.sizeX() + 1, domain.sizeY() + 1, domain.sizeZ() + 1,
      --            true, false) ;
      -- #endif

      -- #pragma omp parallel for firstprivate(numNode)
      --x   for (Index_t i=0; i<numNode; ++i) {
      --x      domain.fx(i) = Real_t(0.0) ;
      --x      domain.fy(i) = Real_t(0.0) ;
      --x      domain.fz(i) = Real_t(0.0) ;
      --x   }
      -- Looping instead of using an aggregate because an aggreagate might be
      -- very large on the stack:
      for node in domain.nodes'Range loop
         domain.nodes(node).force := (others => 0.0);
      end loop;

      ---   /* Calcforce calls partial, force, hourq */
      --   CalcVolumeForceForElems(domain) ;
      CalcVolumeForceForElems(domain);

      -- #if USE_MPI
      --   Domain_member fieldData[3] ;
      --   fieldData[0] = &Domain::fx ;
      --   fieldData[1] = &Domain::fy ;
      --   fieldData[2] = &Domain::fz ;

      --   CommSend(domain, MSG_COMM_SBN, 3, fieldData,
      --            domain.sizeX() + 1, domain.sizeY() + 1, domain.sizeZ() +  1,
      --            true, false) ;
      --   CommSBN(domain, 3, fieldData) ;
      -- #endif
      --x }
   end CalcForceForNodes;

   -- /******************************************/

   -- static inline
   -- void CalcAccelerationForNodes(Domain &domain, Index_t numNode)
   -- {

   -- #pragma omp parallel for firstprivate(numNode)
   --    for (Index_t i = 0; i < numNode; ++i) {
   --       domain.xdd(i) = domain.fx(i) / domain.nodalMass(i);
   --       domain.ydd(i) = domain.fy(i) / domain.nodalMass(i);
   --       domain.zdd(i) = domain.fz(i) / domain.nodalMass(i);
   --    }
   -- }

   -- /******************************************/

   -- static inline
   -- void ApplyAccelerationBoundaryConditionsForNodes(Domain& domain)
   -- {
   --    Index_t size = domain.sizeX();
   --    Index_t numNodeBC = (size+1)*(size+1) ;

   -- #pragma omp parallel
   --    {
   --       if (!domain.symmXempty() != 0) {
   -- #pragma omp for nowait firstprivate(numNodeBC)
   --          for(Index_t i=0 ; i<numNodeBC ; ++i)
   --             domain.xdd(domain.symmX(i)) = Real_t(0.0) ;
   --       }

   --       if (!domain.symmYempty() != 0) {
   -- #pragma omp for nowait firstprivate(numNodeBC)
   --          for(Index_t i=0 ; i<numNodeBC ; ++i)
   --             domain.ydd(domain.symmY(i)) = Real_t(0.0) ;
   --       }

   --       if (!domain.symmZempty() != 0) {
   -- #pragma omp for nowait firstprivate(numNodeBC)
   --          for(Index_t i=0 ; i<numNodeBC ; ++i)
   --             domain.zdd(domain.symmZ(i)) = Real_t(0.0) ;
   --       }
   --    }
   -- }

   -- /******************************************/

   -- static inline
   -- void CalcVelocityForNodes(Domain &domain, const Real_t dt, const Real_t u_cut,
   --                           Index_t numNode)
   -- {

   -- #pragma omp parallel for firstprivate(numNode)
   --    for ( Index_t i = 0 ; i < numNode ; ++i )
   --    {
   --      Real_t xdtmp, ydtmp, zdtmp ;

   --      xdtmp = domain.xd(i) + domain.xdd(i) * dt ;
   --      if( FABS(xdtmp) < u_cut ) xdtmp = Real_t(0.0);
   --      domain.xd(i) = xdtmp ;

   --      ydtmp = domain.yd(i) + domain.ydd(i) * dt ;
   --      if( FABS(ydtmp) < u_cut ) ydtmp = Real_t(0.0);
   --      domain.yd(i) = ydtmp ;

   --      zdtmp = domain.zd(i) + domain.zdd(i) * dt ;
   --      if( FABS(zdtmp) < u_cut ) zdtmp = Real_t(0.0);
   --      domain.zd(i) = zdtmp ;
   --    }
   -- }

   -- /******************************************/

   -- static inline
   -- void CalcPositionForNodes(Domain &domain, const Real_t dt, Index_t numNode)
   -- {
   -- #pragma omp parallel for firstprivate(numNode)
   --    for ( Index_t i = 0 ; i < numNode ; ++i )
   --    {
   --      domain.x(i) += domain.xd(i) * dt ;
   --      domain.y(i) += domain.yd(i) * dt ;
   --      domain.z(i) += domain.zd(i) * dt ;
   --    }
   -- }

   -- /******************************************/

   --x static inline
   --x void LagrangeNodal(Domain& domain)
   procedure LagrangeNodal (domain : in out Domain_Record) is
      --x {
      -- #ifdef SEDOV_SYNC_POS_VEL_EARLY
      --    Domain_member fieldData[6] ;
      -- #endif

      --x    const Real_t delt = domain.deltatime() ;
      --x    Real_t u_cut = domain.u_cut() ;
      delt  : constant Time_Span := domain.variables.deltatime;
      u_cut : Velocity := domain.parameters.velocity_tolerance;
   begin
      ---   /* time of boundary condition evaluation is beginning of step for force and
      ---    * acceleration boundary conditions. */
      --x   CalcForceForNodes(domain);
      CalcForceForNodes(domain);

      -- #if USE_MPI
      -- #ifdef SEDOV_SYNC_POS_VEL_EARLY
      --    CommRecv(domain, MSG_SYNC_POS_VEL, 6,
      --             domain.sizeX() + 1, domain.sizeY() + 1, domain.sizeZ() + 1,
      --             false, false) ;
      -- #endif
      -- #endif

      --x    CalcAccelerationForNodes(domain, domain.numNode());
      --x    ApplyAccelerationBoundaryConditionsForNodes(domain);
      --x    CalcVelocityForNodes( domain, delt, u_cut, domain.numNode()) ;
      --x    CalcPositionForNodes( domain, delt, domain.numNode() );
      CalcAccelerationForNodes (domain, domain.numNode);
      ApplyAccelerationBoundaryConditionsForNodes (domain);
      CalcVelocityForNodes (domain, delt, u_cut, domain.numNode);
      CalcPositionForNodes (domain, delt, domain.numNode);

      -- #if USE_MPI
      -- #ifdef SEDOV_SYNC_POS_VEL_EARLY
      --   fieldData[0] = &Domain::x ;
      --   fieldData[1] = &Domain::y ;
      --   fieldData[2] = &Domain::z ;
      --   fieldData[3] = &Domain::xd ;
      --   fieldData[4] = &Domain::yd ;
      --   fieldData[5] = &Domain::zd ;

      --    CommSend(domain, MSG_SYNC_POS_VEL, 6, fieldData,
      --             domain.sizeX() + 1, domain.sizeY() + 1, domain.sizeZ() + 1,
      --             false, false) ;
      --    CommSyncPosVel(domain) ;
      -- #endif
      -- #endif

      --x   return;
      --x }
   end LagrangeNodal;

   --- /******************************************/

   type EOS_Info_Record (numElem : Element_Count) is record
      bvc          : Element_Real_Array_Access := new Element_Real_Array (0..numElem);
      compHalfStep : Element_Compression_Array_Access := new Element_Compression_Array (0..numElem);
      compression  : Element_Compression_Array_Access := new Element_Compression_Array (0..numElem);
      delvc        : Element_Volume_Array_Access := new Element_Volume_Array (0..numElem);
      e_cut        : Energy;
      e_new        : Element_Energy_Array_Access := new Element_Energy_Array (0..numElem);
      e_old        : Element_Energy_Array_Access := new Element_Energy_Array (0..numElem);
      emin         : Energy;
      eosvmax      : Volume;
      eosvmin      : Volume;
      p_cut        : Pressure;
      p_new        : Element_Pressure_Array_Access := new Element_Pressure_Array (0..numElem);
      p_old        : Element_Pressure_Array_Access := new Element_Pressure_Array (0..numElem);
      pbvc         : Element_Real_Array_Access := new Element_Real_Array (0..numElem);
      pmin         : Pressure;
      q_cut        : Pressure;
      q_new        : Element_Pressure_Array_Access := new Element_Pressure_Array (0..numElem);
      q_old        : Element_Pressure_Array_Access := new Element_Pressure_Array (0..numElem);
      ql_old       : Element_Pressure_Array_Access := new Element_Pressure_Array (0..numElem);
      qq_old       : Element_Pressure_Array_Access := new Element_Pressure_Array (0..numElem);
      rho0         : Density;
      ss4o3        : Real_Type;
      work         : Element_Real_Array_Access := new Element_Real_Array (0..numElem);
   end record;

   -- /******************************************/

   -- static inline
   -- Real_t CalcElemVolume( const Real_t x0, const Real_t x1,
   --                const Real_t x2, const Real_t x3,
   --                const Real_t x4, const Real_t x5,
   --                const Real_t x6, const Real_t x7,
   --                const Real_t y0, const Real_t y1,
   --                const Real_t y2, const Real_t y3,
   --                const Real_t y4, const Real_t y5,
   --                const Real_t y6, const Real_t y7,
   --                const Real_t z0, const Real_t z1,
   --                const Real_t z2, const Real_t z3,
   --                const Real_t z4, const Real_t z5,
   --                const Real_t z6, const Real_t z7 )
   -- {
   --   Real_t twelveth = Real_t(1.0)/Real_t(12.0);

   --   Real_t dx61 = x6 - x1;
   --   Real_t dy61 = y6 - y1;
   --   Real_t dz61 = z6 - z1;

   --   Real_t dx70 = x7 - x0;
   --   Real_t dy70 = y7 - y0;
   --   Real_t dz70 = z7 - z0;

   --   Real_t dx63 = x6 - x3;
   --   Real_t dy63 = y6 - y3;
   --   Real_t dz63 = z6 - z3;

   --   Real_t dx20 = x2 - x0;
   --   Real_t dy20 = y2 - y0;
   --   Real_t dz20 = z2 - z0;

   --   Real_t dx50 = x5 - x0;
   --   Real_t dy50 = y5 - y0;
   --   Real_t dz50 = z5 - z0;

   --   Real_t dx64 = x6 - x4;
   --   Real_t dy64 = y6 - y4;
   --   Real_t dz64 = z6 - z4;

   --   Real_t dx31 = x3 - x1;
   --   Real_t dy31 = y3 - y1;
   --   Real_t dz31 = z3 - z1;

   --   Real_t dx72 = x7 - x2;
   --   Real_t dy72 = y7 - y2;
   --   Real_t dz72 = z7 - z2;

   --   Real_t dx43 = x4 - x3;
   --   Real_t dy43 = y4 - y3;
   --   Real_t dz43 = z4 - z3;

   --   Real_t dx57 = x5 - x7;
   --   Real_t dy57 = y5 - y7;
   --   Real_t dz57 = z5 - z7;

   --   Real_t dx14 = x1 - x4;
   --   Real_t dy14 = y1 - y4;
   --   Real_t dz14 = z1 - z4;

   --   Real_t dx25 = x2 - x5;
   --   Real_t dy25 = y2 - y5;
   --   Real_t dz25 = z2 - z5;

   -- #define TRIPLE_PRODUCT(x1, y1, z1, x2, y2, z2, x3, y3, z3) \
   --    ((x1)*((y2)*(z3) - (z2)*(y3)) + (x2)*((z1)*(y3) - (y1)*(z3)) + (x3)*((y1)*(z2) - (z1)*(y2)))

   --   Real_t volume =
   --     TRIPLE_PRODUCT(dx31 + dx72, dx63, dx20,
   --        dy31 + dy72, dy63, dy20,
   --        dz31 + dz72, dz63, dz20) +
   --     TRIPLE_PRODUCT(dx43 + dx57, dx64, dx70,
   --        dy43 + dy57, dy64, dy70,
   --        dz43 + dz57, dz64, dz70) +
   --     TRIPLE_PRODUCT(dx14 + dx25, dx61, dx50,
   --        dy14 + dy25, dy61, dy50,
   --        dz14 + dz25, dz61, dz50);

   -- #undef TRIPLE_PRODUCT

   --   volume *= twelveth;

   --   return volume ;
   -- }

   -- /******************************************/

   -- //inline
   -- Real_t CalcElemVolume( const Real_t x[8], const Real_t y[8], const Real_t z[8] )
   -- {
   -- return CalcElemVolume( x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7],
   --                        y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7],
   --                        z[0], z[1], z[2], z[3], z[4], z[5], z[6], z[7]);
   -- }

   -- /******************************************/

   -- static inline
   -- Real_t AreaFace( const Real_t x0, const Real_t x1,
   --                  const Real_t x2, const Real_t x3,
   --                  const Real_t y0, const Real_t y1,
   --                  const Real_t y2, const Real_t y3,
   --                  const Real_t z0, const Real_t z1,
   --                  const Real_t z2, const Real_t z3)
   -- {
   --    Real_t fx = (x2 - x0) - (x3 - x1);
   --    Real_t fy = (y2 - y0) - (y3 - y1);
   --    Real_t fz = (z2 - z0) - (z3 - z1);
   --    Real_t gx = (x2 - x0) + (x3 - x1);
   --    Real_t gy = (y2 - y0) + (y3 - y1);
   --    Real_t gz = (z2 - z0) + (z3 - z1);
   --    Real_t area =
   --       (fx * fx + fy * fy + fz * fz) *
   --       (gx * gx + gy * gy + gz * gz) -
   --       (fx * gx + fy * gy + fz * gz) *
   --       (fx * gx + fy * gy + fz * gz);
   --    return area ;
   -- }

   -- /******************************************/

   -- static inline
   -- Real_t CalcElemCharacteristicLength( const Real_t x[8],
   --                                      const Real_t y[8],
   --                                      const Real_t z[8],
   --                                      const Real_t volume)
   -- {
   --    Real_t a, charLength = Real_t(0.0);

   --    a = AreaFace(x[0],x[1],x[2],x[3],
   --                 y[0],y[1],y[2],y[3],
   --                 z[0],z[1],z[2],z[3]) ;
   --    charLength = std::max(a,charLength) ;

   --    a = AreaFace(x[4],x[5],x[6],x[7],
   --                 y[4],y[5],y[6],y[7],
   --                 z[4],z[5],z[6],z[7]) ;
   --    charLength = std::max(a,charLength) ;

   --    a = AreaFace(x[0],x[1],x[5],x[4],
   --                 y[0],y[1],y[5],y[4],
   --                 z[0],z[1],z[5],z[4]) ;
   --    charLength = std::max(a,charLength) ;

   --    a = AreaFace(x[1],x[2],x[6],x[5],
   --                 y[1],y[2],y[6],y[5],
   --                 z[1],z[2],z[6],z[5]) ;
   --    charLength = std::max(a,charLength) ;

   --    a = AreaFace(x[2],x[3],x[7],x[6],
   --                 y[2],y[3],y[7],y[6],
   --                 z[2],z[3],z[7],z[6]) ;
   --    charLength = std::max(a,charLength) ;

   --    a = AreaFace(x[3],x[0],x[4],x[7],
   --                 y[3],y[0],y[4],y[7],
   --                 z[3],z[0],z[4],z[7]) ;
   --    charLength = std::max(a,charLength) ;

   --    charLength = Real_t(4.0) * volume / SQRT(charLength);

   --    return charLength;
   -- }

   -- /******************************************/

   -- static inline
   -- void CalcElemVelocityGradient( const Real_t* const xvel,
   --                                 const Real_t* const yvel,
   --                                 const Real_t* const zvel,
   --                                 const Real_t b[][8],
   --                                 const Real_t detJ,
   --                                 Real_t* const d )
   -- {
   --   const Real_t inv_detJ = Real_t(1.0) / detJ ;
   --   Real_t dyddx, dxddy, dzddx, dxddz, dzddy, dyddz;
   --   const Real_t* const pfx = b[0];
   --   const Real_t* const pfy = b[1];
   --   const Real_t* const pfz = b[2];

   --   d[0] = inv_detJ * ( pfx[0] * (xvel[0]-xvel[6])
   --                      + pfx[1] * (xvel[1]-xvel[7])
   --                      + pfx[2] * (xvel[2]-xvel[4])
   --                      + pfx[3] * (xvel[3]-xvel[5]) );

   --   d[1] = inv_detJ * ( pfy[0] * (yvel[0]-yvel[6])
   --                      + pfy[1] * (yvel[1]-yvel[7])
   --                      + pfy[2] * (yvel[2]-yvel[4])
   --                      + pfy[3] * (yvel[3]-yvel[5]) );

   --   d[2] = inv_detJ * ( pfz[0] * (zvel[0]-zvel[6])
   --                      + pfz[1] * (zvel[1]-zvel[7])
   --                      + pfz[2] * (zvel[2]-zvel[4])
   --                      + pfz[3] * (zvel[3]-zvel[5]) );

   --   dyddx  = inv_detJ * ( pfx[0] * (yvel[0]-yvel[6])
   --                       + pfx[1] * (yvel[1]-yvel[7])
   --                       + pfx[2] * (yvel[2]-yvel[4])
   --                       + pfx[3] * (yvel[3]-yvel[5]) );

   --   dxddy  = inv_detJ * ( pfy[0] * (xvel[0]-xvel[6])
   --                       + pfy[1] * (xvel[1]-xvel[7])
   --                       + pfy[2] * (xvel[2]-xvel[4])
   --                       + pfy[3] * (xvel[3]-xvel[5]) );

   --   dzddx  = inv_detJ * ( pfx[0] * (zvel[0]-zvel[6])
   --                       + pfx[1] * (zvel[1]-zvel[7])
   --                       + pfx[2] * (zvel[2]-zvel[4])
   --                       + pfx[3] * (zvel[3]-zvel[5]) );

   --   dxddz  = inv_detJ * ( pfz[0] * (xvel[0]-xvel[6])
   --                       + pfz[1] * (xvel[1]-xvel[7])
   --                       + pfz[2] * (xvel[2]-xvel[4])
   --                       + pfz[3] * (xvel[3]-xvel[5]) );

   --   dzddy  = inv_detJ * ( pfy[0] * (zvel[0]-zvel[6])
   --                       + pfy[1] * (zvel[1]-zvel[7])
   --                       + pfy[2] * (zvel[2]-zvel[4])
   --                       + pfy[3] * (zvel[3]-zvel[5]) );

   --   dyddz  = inv_detJ * ( pfz[0] * (yvel[0]-yvel[6])
   --                       + pfz[1] * (yvel[1]-yvel[7])
   --                       + pfz[2] * (yvel[2]-yvel[4])
   --                       + pfz[3] * (yvel[3]-yvel[5]) );
   --   d[5]  = Real_t( .5) * ( dxddy + dyddx );
   --   d[4]  = Real_t( .5) * ( dxddz + dzddx );
   --   d[3]  = Real_t( .5) * ( dzddy + dyddz );
   -- }

   -- /******************************************/

   -- //static inline
   -- void CalcKinematicsForElems( Domain &domain, Real_t *vnew,
   --                              Real_t deltaTime, Index_t numElem )
   -- {

   --   // loop over all elements
   -- #pragma omp parallel for firstprivate(numElem, deltaTime)
   --   for( Index_t k=0 ; k<numElem ; ++k )
   --   {
   --     Real_t B[3][8] ; /** shape function derivatives */
   --     Real_t D[6] ;
   --     Real_t x_local[8] ;
   --     Real_t y_local[8] ;
   --     Real_t z_local[8] ;
   --     Real_t xd_local[8] ;
   --     Real_t yd_local[8] ;
   --     Real_t zd_local[8] ;
   --     Real_t detJ = Real_t(0.0) ;

   --     Real_t volume ;
   --     Real_t relativeVolume ;
   --     const Index_t* const elemToNode = domain.nodelist(k) ;

   --     // get nodal coordinates from global arrays and copy into local arrays.
   --     CollectDomainNodesToElemNodes(domain, elemToNode, x_local, y_local, z_local);

   --     // volume calculations
   --     volume = CalcElemVolume(x_local, y_local, z_local );
   --     relativeVolume = volume / domain.volo(k) ;
   --     vnew[k] = relativeVolume ;
   --     domain.delv(k) = relativeVolume - domain.v(k) ;

   --     // set characteristic length
   --     domain.arealg(k) = CalcElemCharacteristicLength(x_local, y_local, z_local,
   --                                              volume);

   --     // get nodal velocities from global array and copy into local arrays.
   --     for( Index_t lnode=0 ; lnode<8 ; ++lnode )
   --     {
   --       Index_t gnode = elemToNode[lnode];
   --       xd_local[lnode] = domain.xd(gnode);
   --       yd_local[lnode] = domain.yd(gnode);
   --       zd_local[lnode] = domain.zd(gnode);
   --     }

   --     Real_t dt2 = Real_t(0.5) * deltaTime;
   --     for ( Index_t j=0 ; j<8 ; ++j )
   --     {
   --        x_local[j] -= dt2 * xd_local[j];
   --        y_local[j] -= dt2 * yd_local[j];
   --        z_local[j] -= dt2 * zd_local[j];
   --     }

   --     CalcElemShapeFunctionDerivatives( x_local, y_local, z_local,
   --                                       B, &detJ );

   --     CalcElemVelocityGradient( xd_local, yd_local, zd_local,
   --                                B, detJ, D );

   --     // put velocity gradient quantities into their global arrays.
   --     domain.dxx(k) = D[0];
   --     domain.dyy(k) = D[1];
   --     domain.dzz(k) = D[2];
   --   }
   -- }

   -- /******************************************/

   -- static inline
   -- void CalcLagrangeElements(Domain& domain, Real_t* vnew)
   -- {
   --    Index_t numElem = domain.numElem() ;
   --    if (numElem > 0) {
   --       const Real_t deltatime = domain.deltatime() ;

   --       domain.AllocateStrains(numElem);

   --       CalcKinematicsForElems(domain, vnew, deltatime, numElem) ;

   --       // element loop to do some stuff not included in the elemlib function.
   -- #pragma omp parallel for firstprivate(numElem)
   --       for ( Index_t k=0 ; k<numElem ; ++k )
   --       {
   --          // calc strain rate and apply as constraint (only done in FB element)
   --          Real_t vdov = domain.dxx(k) + domain.dyy(k) + domain.dzz(k) ;
   --          Real_t vdovthird = vdov/Real_t(3.0) ;

   --          // make the rate of deformation tensor deviatoric
   --          domain.vdov(k) = vdov ;
   --          domain.dxx(k) -= vdovthird ;
   --          domain.dyy(k) -= vdovthird ;
   --          domain.dzz(k) -= vdovthird ;

   --         // See if any volumes are negative, and take appropriate action.
   --          if (vnew[k] <= Real_t(0.0))
   --         {
   -- #if USE_MPI
   --            MPI_Abort(MPI_COMM_WORLD, VolumeError) ;
   -- #else
   --            exit(VolumeError);
   -- #endif
   --         }
   --       }
   --       domain.DeallocateStrains();
   --    }
   -- }

   -- /******************************************/

   -- static inline
   -- void CalcMonotonicQGradientsForElems(Domain& domain, Real_t vnew[])
   -- {
   --    Index_t numElem = domain.numElem();

   -- #pragma omp parallel for firstprivate(numElem)
   --    for (Index_t i = 0 ; i < numElem ; ++i ) {
   --       const Real_t ptiny = Real_t(1.e-36) ;
   --       Real_t ax,ay,az ;
   --       Real_t dxv,dyv,dzv ;

   --       const Index_t *elemToNode = domain.nodelist(i);
   --       Index_t n0 = elemToNode[0] ;
   --       Index_t n1 = elemToNode[1] ;
   --       Index_t n2 = elemToNode[2] ;
   --       Index_t n3 = elemToNode[3] ;
   --       Index_t n4 = elemToNode[4] ;
   --       Index_t n5 = elemToNode[5] ;
   --       Index_t n6 = elemToNode[6] ;
   --       Index_t n7 = elemToNode[7] ;

   --       Real_t x0 = domain.x(n0) ;
   --       Real_t x1 = domain.x(n1) ;
   --       Real_t x2 = domain.x(n2) ;
   --       Real_t x3 = domain.x(n3) ;
   --       Real_t x4 = domain.x(n4) ;
   --       Real_t x5 = domain.x(n5) ;
   --       Real_t x6 = domain.x(n6) ;
   --       Real_t x7 = domain.x(n7) ;

   --       Real_t y0 = domain.y(n0) ;
   --       Real_t y1 = domain.y(n1) ;
   --       Real_t y2 = domain.y(n2) ;
   --       Real_t y3 = domain.y(n3) ;
   --       Real_t y4 = domain.y(n4) ;
   --       Real_t y5 = domain.y(n5) ;
   --       Real_t y6 = domain.y(n6) ;
   --       Real_t y7 = domain.y(n7) ;

   --       Real_t z0 = domain.z(n0) ;
   --       Real_t z1 = domain.z(n1) ;
   --       Real_t z2 = domain.z(n2) ;
   --       Real_t z3 = domain.z(n3) ;
   --       Real_t z4 = domain.z(n4) ;
   --       Real_t z5 = domain.z(n5) ;
   --       Real_t z6 = domain.z(n6) ;
   --       Real_t z7 = domain.z(n7) ;

   --       Real_t xv0 = domain.xd(n0) ;
   --       Real_t xv1 = domain.xd(n1) ;
   --       Real_t xv2 = domain.xd(n2) ;
   --       Real_t xv3 = domain.xd(n3) ;
   --       Real_t xv4 = domain.xd(n4) ;
   --       Real_t xv5 = domain.xd(n5) ;
   --       Real_t xv6 = domain.xd(n6) ;
   --       Real_t xv7 = domain.xd(n7) ;

   --       Real_t yv0 = domain.yd(n0) ;
   --       Real_t yv1 = domain.yd(n1) ;
   --       Real_t yv2 = domain.yd(n2) ;
   --       Real_t yv3 = domain.yd(n3) ;
   --       Real_t yv4 = domain.yd(n4) ;
   --       Real_t yv5 = domain.yd(n5) ;
   --       Real_t yv6 = domain.yd(n6) ;
   --       Real_t yv7 = domain.yd(n7) ;

   --       Real_t zv0 = domain.zd(n0) ;
   --       Real_t zv1 = domain.zd(n1) ;
   --       Real_t zv2 = domain.zd(n2) ;
   --       Real_t zv3 = domain.zd(n3) ;
   --       Real_t zv4 = domain.zd(n4) ;
   --       Real_t zv5 = domain.zd(n5) ;
   --       Real_t zv6 = domain.zd(n6) ;
   --       Real_t zv7 = domain.zd(n7) ;

   --       Real_t vol = domain.volo(i)*vnew[i] ;
   --       Real_t norm = Real_t(1.0) / ( vol + ptiny ) ;

   --       Real_t dxj = Real_t(-0.25)*((x0+x1+x5+x4) - (x3+x2+x6+x7)) ;
   --       Real_t dyj = Real_t(-0.25)*((y0+y1+y5+y4) - (y3+y2+y6+y7)) ;
   --       Real_t dzj = Real_t(-0.25)*((z0+z1+z5+z4) - (z3+z2+z6+z7)) ;

   --       Real_t dxi = Real_t( 0.25)*((x1+x2+x6+x5) - (x0+x3+x7+x4)) ;
   --       Real_t dyi = Real_t( 0.25)*((y1+y2+y6+y5) - (y0+y3+y7+y4)) ;
   --       Real_t dzi = Real_t( 0.25)*((z1+z2+z6+z5) - (z0+z3+z7+z4)) ;

   --       Real_t dxk = Real_t( 0.25)*((x4+x5+x6+x7) - (x0+x1+x2+x3)) ;
   --       Real_t dyk = Real_t( 0.25)*((y4+y5+y6+y7) - (y0+y1+y2+y3)) ;
   --       Real_t dzk = Real_t( 0.25)*((z4+z5+z6+z7) - (z0+z1+z2+z3)) ;

   --       /* find delvk and delxk ( i cross j ) */

   --       ax = dyi*dzj - dzi*dyj ;
   --       ay = dzi*dxj - dxi*dzj ;
   --       az = dxi*dyj - dyi*dxj ;

   --       domain.delx_zeta(i) = vol / SQRT(ax*ax + ay*ay + az*az + ptiny) ;

   --       ax *= norm ;
   --       ay *= norm ;
   --       az *= norm ;

   --       dxv = Real_t(0.25)*((xv4+xv5+xv6+xv7) - (xv0+xv1+xv2+xv3)) ;
   --       dyv = Real_t(0.25)*((yv4+yv5+yv6+yv7) - (yv0+yv1+yv2+yv3)) ;
   --       dzv = Real_t(0.25)*((zv4+zv5+zv6+zv7) - (zv0+zv1+zv2+zv3)) ;

   --       domain.delv_zeta(i) = ax*dxv + ay*dyv + az*dzv ;

   --       /* find delxi and delvi ( j cross k ) */

   --       ax = dyj*dzk - dzj*dyk ;
   --       ay = dzj*dxk - dxj*dzk ;
   --       az = dxj*dyk - dyj*dxk ;

   --       domain.delx_xi(i) = vol / SQRT(ax*ax + ay*ay + az*az + ptiny) ;

   --       ax *= norm ;
   --       ay *= norm ;
   --       az *= norm ;

   --       dxv = Real_t(0.25)*((xv1+xv2+xv6+xv5) - (xv0+xv3+xv7+xv4)) ;
   --       dyv = Real_t(0.25)*((yv1+yv2+yv6+yv5) - (yv0+yv3+yv7+yv4)) ;
   --       dzv = Real_t(0.25)*((zv1+zv2+zv6+zv5) - (zv0+zv3+zv7+zv4)) ;

   --       domain.delv_xi(i) = ax*dxv + ay*dyv + az*dzv ;

   --       /* find delxj and delvj ( k cross i ) */

   --       ax = dyk*dzi - dzk*dyi ;
   --       ay = dzk*dxi - dxk*dzi ;
   --       az = dxk*dyi - dyk*dxi ;

   --       domain.delx_eta(i) = vol / SQRT(ax*ax + ay*ay + az*az + ptiny) ;

   --       ax *= norm ;
   --       ay *= norm ;
   --       az *= norm ;

   --       dxv = Real_t(-0.25)*((xv0+xv1+xv5+xv4) - (xv3+xv2+xv6+xv7)) ;
   --       dyv = Real_t(-0.25)*((yv0+yv1+yv5+yv4) - (yv3+yv2+yv6+yv7)) ;
   --       dzv = Real_t(-0.25)*((zv0+zv1+zv5+zv4) - (zv3+zv2+zv6+zv7)) ;

   --       domain.delv_eta(i) = ax*dxv + ay*dyv + az*dzv ;
   --    }
   -- }

   -- /******************************************/

   -- static inline
   -- void CalcMonotonicQRegionForElems(Domain &domain, Int_t r,
   --                                   Real_t vnew[], Real_t ptiny)
   -- {
   --    Real_t monoq_limiter_mult = domain.monoq_limiter_mult();
   --    Real_t monoq_max_slope = domain.monoq_max_slope();
   --    Real_t qlc_monoq = domain.qlc_monoq();
   --    Real_t qqc_monoq = domain.qqc_monoq();

   -- #pragma omp parallel for firstprivate(qlc_monoq, qqc_monoq, monoq_limiter_mult, monoq_max_slope, ptiny)
   --    for ( Index_t ielem = 0 ; ielem < domain.regElemSize(r); ++ielem ) {
   --       Index_t i = domain.regElemlist(r,ielem);
   --       Real_t qlin, qquad ;
   --       Real_t phixi, phieta, phizeta ;
   --       Int_t bcMask = domain.elemBC(i) ;
   --       Real_t delvm, delvp ;

   --       /*  phixi     */
   --       Real_t norm = Real_t(1.) / (domain.delv_xi(i)+ ptiny ) ;

   --       switch (bcMask & XI_M) {
   --          case XI_M_COMM: /* needs comm data */
   --          case 0:         delvm = domain.delv_xi(domain.lxim(i)); break ;
   --          case XI_M_SYMM: delvm = domain.delv_xi(i) ;       break ;
   --          case XI_M_FREE: delvm = Real_t(0.0) ;      break ;
   --          default:        /* ERROR */ ;              break ;
   --       }
   --       switch (bcMask & XI_P) {
   --          case XI_P_COMM: /* needs comm data */
   --          case 0:         delvp = domain.delv_xi(domain.lxip(i)) ; break ;
   --          case XI_P_SYMM: delvp = domain.delv_xi(i) ;       break ;
   --          case XI_P_FREE: delvp = Real_t(0.0) ;      break ;
   --          default:        /* ERROR */ ;              break ;
   --       }

   --       delvm = delvm * norm ;
   --       delvp = delvp * norm ;

   --       phixi = Real_t(.5) * ( delvm + delvp ) ;

   --       delvm *= monoq_limiter_mult ;
   --       delvp *= monoq_limiter_mult ;

   --       if ( delvm < phixi ) phixi = delvm ;
   --       if ( delvp < phixi ) phixi = delvp ;
   --       if ( phixi < Real_t(0.)) phixi = Real_t(0.) ;
   --       if ( phixi > monoq_max_slope) phixi = monoq_max_slope;


   --       /*  phieta     */
   --       norm = Real_t(1.) / ( domain.delv_eta(i) + ptiny ) ;

   --       switch (bcMask & ETA_M) {
   --          case ETA_M_COMM: /* needs comm data */
   --          case 0:          delvm = domain.delv_eta(domain.letam(i)) ; break ;
   --          case ETA_M_SYMM: delvm = domain.delv_eta(i) ;        break ;
   --          case ETA_M_FREE: delvm = Real_t(0.0) ;        break ;
   --          default:         /* ERROR */ ;                break ;
   --       }
   --       switch (bcMask & ETA_P) {
   --          case ETA_P_COMM: /* needs comm data */
   --          case 0:          delvp = domain.delv_eta(domain.letap(i)) ; break ;
   --          case ETA_P_SYMM: delvp = domain.delv_eta(i) ;        break ;
   --          case ETA_P_FREE: delvp = Real_t(0.0) ;        break ;
   --          default:         /* ERROR */ ;                break ;
   --       }

   --       delvm = delvm * norm ;
   --       delvp = delvp * norm ;

   --       phieta = Real_t(.5) * ( delvm + delvp ) ;

   --       delvm *= monoq_limiter_mult ;
   --       delvp *= monoq_limiter_mult ;

   --       if ( delvm  < phieta ) phieta = delvm ;
   --       if ( delvp  < phieta ) phieta = delvp ;
   --       if ( phieta < Real_t(0.)) phieta = Real_t(0.) ;
   --       if ( phieta > monoq_max_slope)  phieta = monoq_max_slope;

   --       /*  phizeta     */
   --       norm = Real_t(1.) / ( domain.delv_zeta(i) + ptiny ) ;

   --       switch (bcMask & ZETA_M) {
   --          case ZETA_M_COMM: /* needs comm data */
   --          case 0:           delvm = domain.delv_zeta(domain.lzetam(i)) ; break ;
   --          case ZETA_M_SYMM: delvm = domain.delv_zeta(i) ;         break ;
   --          case ZETA_M_FREE: delvm = Real_t(0.0) ;          break ;
   --          default:          /* ERROR */ ;                  break ;
   --       }
   --       switch (bcMask & ZETA_P) {
   --          case ZETA_P_COMM: /* needs comm data */
   --          case 0:           delvp = domain.delv_zeta(domain.lzetap(i)) ; break ;
   --          case ZETA_P_SYMM: delvp = domain.delv_zeta(i) ;         break ;
   --          case ZETA_P_FREE: delvp = Real_t(0.0) ;          break ;
   --          default:          /* ERROR */ ;                  break ;
   --       }

   --       delvm = delvm * norm ;
   --       delvp = delvp * norm ;

   --       phizeta = Real_t(.5) * ( delvm + delvp ) ;

   --       delvm *= monoq_limiter_mult ;
   --       delvp *= monoq_limiter_mult ;

   --       if ( delvm   < phizeta ) phizeta = delvm ;
   --       if ( delvp   < phizeta ) phizeta = delvp ;
   --       if ( phizeta < Real_t(0.)) phizeta = Real_t(0.);
   --       if ( phizeta > monoq_max_slope  ) phizeta = monoq_max_slope;

   --       /* Remove length scale */

   --       if ( domain.vdov(i) > Real_t(0.) )  {
   --          qlin  = Real_t(0.) ;
   --          qquad = Real_t(0.) ;
   --       }
   --       else {
   --          Real_t delvxxi   = domain.delv_xi(i)   * domain.delx_xi(i)   ;
   --          Real_t delvxeta  = domain.delv_eta(i)  * domain.delx_eta(i)  ;
   --          Real_t delvxzeta = domain.delv_zeta(i) * domain.delx_zeta(i) ;

   --          if ( delvxxi   > Real_t(0.) ) delvxxi   = Real_t(0.) ;
   --          if ( delvxeta  > Real_t(0.) ) delvxeta  = Real_t(0.) ;
   --          if ( delvxzeta > Real_t(0.) ) delvxzeta = Real_t(0.) ;

   --          Real_t rho = domain.elemMass(i) / (domain.volo(i) * vnew[i]) ;

   --          qlin = -qlc_monoq * rho *
   --             (  delvxxi   * (Real_t(1.) - phixi) +
   --                delvxeta  * (Real_t(1.) - phieta) +
   --                delvxzeta * (Real_t(1.) - phizeta)  ) ;

   --          qquad = qqc_monoq * rho *
   --             (  delvxxi*delvxxi     * (Real_t(1.) - phixi*phixi) +
   --                delvxeta*delvxeta   * (Real_t(1.) - phieta*phieta) +
   --                delvxzeta*delvxzeta * (Real_t(1.) - phizeta*phizeta)  ) ;
   --       }

   --       domain.qq(i) = qquad ;
   --       domain.ql(i) = qlin  ;
   --    }
   -- }

   -- /******************************************/

   -- static inline
   -- void CalcMonotonicQForElems(Domain& domain, Real_t vnew[])
   -- {
   --    //
   --    // initialize parameters
   --    //
   --    const Real_t ptiny = Real_t(1.e-36) ;

   --    //
   --    // calculate the monotonic q for all regions
   --    //
   --    for (Index_t r=0 ; r<domain.numReg() ; ++r) {

   --       if (domain.regElemSize(r) > 0) {
   --          CalcMonotonicQRegionForElems(domain, r, vnew, ptiny) ;
   --       }
   --    }
   -- }

   -- /******************************************/

   -- static inline
   -- void CalcQForElems(Domain& domain, Real_t vnew[])
   -- {
   --    //
   --    // MONOTONIC Q option
   --    //

   --    Index_t numElem = domain.numElem() ;

   --    if (numElem != 0) {
   --       Int_t allElem = numElem +  /* local elem */
   --             2*domain.sizeX()*domain.sizeY() + /* plane ghosts */
   --             2*domain.sizeX()*domain.sizeZ() + /* row ghosts */
   --             2*domain.sizeY()*domain.sizeZ() ; /* col ghosts */

   --       domain.AllocateGradients(numElem, allElem);

   -- #if USE_MPI
   --       CommRecv(domain, MSG_MONOQ, 3,
   --                domain.sizeX(), domain.sizeY(), domain.sizeZ(),
   --                true, true) ;
   -- #endif

   --       /* Calculate velocity gradients */
   --       CalcMonotonicQGradientsForElems(domain, vnew);

   -- #if USE_MPI
   --       Domain_member fieldData[3] ;

   --       /* Transfer veloctiy gradients in the first order elements */
   --       /* problem->commElements->Transfer(CommElements::monoQ) ; */

   --       fieldData[0] = &Domain::delv_xi ;
   --       fieldData[1] = &Domain::delv_eta ;
   --       fieldData[2] = &Domain::delv_zeta ;

   --       CommSend(domain, MSG_MONOQ, 3, fieldData,
   --                domain.sizeX(), domain.sizeY(), domain.sizeZ(),
   --                true, true) ;

   --       CommMonoQ(domain) ;
   -- #endif

   --       CalcMonotonicQForElems(domain, vnew) ;

   --       // Free up memory
   --       domain.DeallocateGradients();

   --       /* Don't allow excessive artificial viscosity */
   --       Index_t idx = -1;
   --       for (Index_t i=0; i<numElem; ++i) {
   --          if ( domain.q(i) > domain.qstop() ) {
   --             idx = i ;
   --             break ;
   --          }
   --       }

   --       if(idx >= 0) {
   -- #if USE_MPI
   --          MPI_Abort(MPI_COMM_WORLD, QStopError) ;
   -- #else
   --          exit(QStopError);
   -- #endif
   --       }
   --    }
   -- }

   -- /******************************************/

   -- static inline
   -- void CalcPressureForElems(Real_t* p_new, Real_t* bvc,
   --                           Real_t* pbvc, Real_t* e_old,
   --                           Real_t* compression, Real_t *vnewc,
   --                           Real_t pmin,
   --                           Real_t p_cut, Real_t eosvmax,
   --                           Index_t length, Index_t *regElemList)
   -- {
   -- #pragma omp parallel for firstprivate(length)
   --    for (Index_t i = 0; i < length ; ++i) {
   --       Real_t c1s = Real_t(2.0)/Real_t(3.0) ;
   --       bvc[i] = c1s * (compression[i] + Real_t(1.));
   --       pbvc[i] = c1s;
   --    }

   -- #pragma omp parallel for firstprivate(length, pmin, p_cut, eosvmax)
   --    for (Index_t i = 0 ; i < length ; ++i){
   --       Index_t elem = regElemList[i];

   --       p_new[i] = bvc[i] * e_old[i] ;

   --       if    (FABS(p_new[i]) <  p_cut   )
   --          p_new[i] = Real_t(0.0) ;

   --       if    ( vnewc[elem] >= eosvmax ) /* impossible condition here? */
   --          p_new[i] = Real_t(0.0) ;

   --       if    (p_new[i]       <  pmin)
   --          p_new[i]   = pmin ;
   --    }
   -- }

   --- /******************************************/

   -- static inline
   -- void CalcEnergyForElems(Real_t* p_new, Real_t* e_new, Real_t* q_new,
   --                         Real_t* bvc, Real_t* pbvc,
   --                         Real_t* p_old, Real_t* e_old, Real_t* q_old,
   --                         Real_t* compression, Real_t* compHalfStep,
   --                         Real_t* vnewc, Real_t* work, Real_t* delvc, Real_t pmin,
   --                         Real_t p_cut, Real_t  e_cut, Real_t q_cut, Real_t emin,
   --                         Real_t* qq_old, Real_t* ql_old,
   --                         Real_t rho0,
   --                         Real_t eosvmax,
   --                         Index_t length, Index_t *regElemList)
   -- {
   procedure CalcEnergyForElems
     (domain      : in out Domain_Record;
      vnewc       : access Element_Volume_Array;
      info        : in EOS_Info_Record;
      regElemList : access Element_Element_Index_Array)
    --       (domain      : in out Domain_Record;
--        p_new   : access Element_Pressure_Array;
--        e_new : access Element_Energy_Array;
--         q_new : access Element_Pressure_Array;
--        bvc: access Element_Real_Array;
--        pbvc: access Element_Real_Array;
--       p_old : access Element_Pressure_Array;
--         e_old: access Element_Energy_Array;
--         q_old : access Element_Pressure_Array;
--        Real_t* compression : access Element_Compression_Array;
--        Real_t* compHalfStep,
--         vnewc : access Element_Volume_Array;
--        Real_t* work,
--        Real_t* delvc,
--        Real_t* qq_old,
--        Real_t* ql_old,
--        Index_t length,
--        regElemList : access Element_Element_Index_Array)
   with inline
   is
--        Real_t pmin,
--        Real_t p_cut,
--        Real_t  e_cut,
--        Real_t q_cut,
--        Real_t emin,
--        rho0        : in Density;
--        Real_t eosvmax,
   --    Real_t *pHalfStep = Allocate<Real_t>(length) ;
   begin
null;
   -- #pragma omp parallel for firstprivate(length, emin)
   --    for (Index_t i = 0 ; i < length ; ++i) {
   --       e_new[i] = e_old[i] - Real_t(0.5) * delvc[i] * (p_old[i] + q_old[i])
   --          + Real_t(0.5) * work[i];

   --       if (e_new[i]  < emin ) {
   --          e_new[i] = emin ;
   --       }
   --    }

   --    CalcPressureForElems(pHalfStep, bvc, pbvc, e_new, compHalfStep, vnewc,
   --                         pmin, p_cut, eosvmax, length, regElemList);

   -- #pragma omp parallel for firstprivate(length, rho0)
   --    for (Index_t i = 0 ; i < length ; ++i) {
   --       Real_t vhalf = Real_t(1.) / (Real_t(1.) + compHalfStep[i]) ;

   --       if ( delvc[i] > Real_t(0.) ) {
   --          q_new[i] /* = qq_old[i] = ql_old[i] */ = Real_t(0.) ;
   --       }
   --       else {
   --          Real_t ssc = ( pbvc[i] * e_new[i]
   --                  + vhalf * vhalf * bvc[i] * pHalfStep[i] ) / rho0 ;

   --          if ( ssc <= Real_t(.1111111e-36) ) {
   --             ssc = Real_t(.3333333e-18) ;
   --          } else {
   --             ssc = SQRT(ssc) ;
   --          }

   --          q_new[i] = (ssc*ql_old[i] + qq_old[i]) ;
   --       }

   --       e_new[i] = e_new[i] + Real_t(0.5) * delvc[i]
   --          * (  Real_t(3.0)*(p_old[i]     + q_old[i])
   --               - Real_t(4.0)*(pHalfStep[i] + q_new[i])) ;
   --    }

   -- #pragma omp parallel for firstprivate(length, emin, e_cut)
   --    for (Index_t i = 0 ; i < length ; ++i) {

   --       e_new[i] += Real_t(0.5) * work[i];

   --       if (FABS(e_new[i]) < e_cut) {
   --          e_new[i] = Real_t(0.)  ;
   --       }
   --       if (     e_new[i]  < emin ) {
   --          e_new[i] = emin ;
   --       }
   --    }

   --    CalcPressureForElems(p_new, bvc, pbvc, e_new, compression, vnewc,
   --                         pmin, p_cut, eosvmax, length, regElemList);

   -- #pragma omp parallel for firstprivate(length, rho0, emin, e_cut)
   --    for (Index_t i = 0 ; i < length ; ++i){
   --       const Real_t sixth = Real_t(1.0) / Real_t(6.0) ;
   --       Index_t elem = regElemList[i];
   --       Real_t q_tilde ;

   --       if (delvc[i] > Real_t(0.)) {
   --          q_tilde = Real_t(0.) ;
   --       }
   --       else {
   --          Real_t ssc = ( pbvc[i] * e_new[i]
   --                  + vnewc[elem] * vnewc[elem] * bvc[i] * p_new[i] ) / rho0 ;

   --          if ( ssc <= Real_t(.1111111e-36) ) {
   --             ssc = Real_t(.3333333e-18) ;
   --          } else {
   --             ssc = SQRT(ssc) ;
   --          }

   --          q_tilde = (ssc*ql_old[i] + qq_old[i]) ;
   --       }

   --       e_new[i] = e_new[i] - (  Real_t(7.0)*(p_old[i]     + q_old[i])
   --                                - Real_t(8.0)*(pHalfStep[i] + q_new[i])
   --                                + (p_new[i] + q_tilde)) * delvc[i]*sixth ;

   --       if (FABS(e_new[i]) < e_cut) {
   --          e_new[i] = Real_t(0.)  ;
   --       }
   --       if (     e_new[i]  < emin ) {
   --          e_new[i] = emin ;
   --       }
   --    }

   --    CalcPressureForElems(p_new, bvc, pbvc, e_new, compression, vnewc,
   --                         pmin, p_cut, eosvmax, length, regElemList);

   -- #pragma omp parallel for firstprivate(length, rho0, q_cut)
   --    for (Index_t i = 0 ; i < length ; ++i){
   --       Index_t elem = regElemList[i];

   --       if ( delvc[i] <= Real_t(0.) ) {
   --          Real_t ssc = ( pbvc[i] * e_new[i]
   --                  + vnewc[elem] * vnewc[elem] * bvc[i] * p_new[i] ) / rho0 ;

   --          if ( ssc <= Real_t(.1111111e-36) ) {
   --             ssc = Real_t(.3333333e-18) ;
   --          } else {
   --             ssc = SQRT(ssc) ;
   --          }

   --          q_new[i] = (ssc*ql_old[i] + qq_old[i]) ;

   --          if (FABS(q_new[i]) < q_cut) q_new[i] = Real_t(0.) ;
   --       }
   --    }

   --    Release(&pHalfStep) ;

   --x    return ;
   --x }
   end CalcEnergyForElems;

   --- /******************************************/

   --x static inline
   --x void CalcSoundSpeedForElems(Domain &domain,
   --x                             Real_t *vnewc, Real_t rho0, Real_t *enewc,
   --x                             Real_t *pnewc, Real_t *pbvc,
   --x                             Real_t *bvc, Real_t ss4o3,
   --x                             Index_t len, Index_t *regElemList)
   --x {
   procedure CalcSoundSpeedForElems
     (domain      : in out Domain_Record;
      vnewc       : access Element_Volume_Array;
      rho0        : in Density;
      enewc       : access Element_Energy_Array;
      pnewc       : access Element_Pressure_Array;
      pbvc        : access Element_Real_Array;
      bvc         : access Element_Real_Array;
      ss4o3       : in Real_Type;
      regElemList : access Element_Element_Index_Array)
     with inline is
   begin
   -- #pragma omp parallel for firstprivate(rho0, ss4o3)
   --x    for (Index_t i = 0; i < len ; ++i) {
   --x       Index_t elem = regElemList[i];
   --x       Real_t ssTmp = (pbvc[i] * enewc[i] + vnewc[elem] * vnewc[elem] *
   --x                  bvc[i] * pnewc[i]) / rho0;
   --x       if (ssTmp <= Real_t(.1111111e-36)) {
   --x          ssTmp = Real_t(.3333333e-18);
   --x       }
   --x       else {
   --x          ssTmp = SQRT(ssTmp);
   --x       }
   --x       domain.ss(elem) = ssTmp ;
   --x    }
      for region_element in regElemList'Range loop
         declare
            element : constant Element_Index := regElemList(region_element);
            ssTmp : Real_Type :=
              (pbvc(region_element) * Real_Type(enewc(region_element)) +
                   Real_Type(vnewc(element))**2 * bvc(region_element) *
                   Real_Type(pnewc(region_element))) /
                Real_Type(rho0);
         begin
            if ssTmp <= 0.1111111e-36 then
               ssTmp := 0.3333333e-18;
            else
               ssTmp := SQRT(ssTmp);
            end if;
            domain.elements(element).sound_speed := Velocity (ssTmp);
         end;
     end loop;
   --x }
   end CalcSoundSpeedForElems;

   --x static inline
   --x void EvalEOSForElems(Domain& domain, Real_t *vnewc,
   --x                      Int_t numElemReg, Index_t *regElemList, Int_t rep)
   --x {
   procedure EvalEOSForElems
     (domain      : in out Domain_Record;
      vnewc       : access Element_Volume_Array;
      numElemReg  : in     Element_Count;
      regElemList : access Element_Element_Index_Array;
      rep         : in     Cost_Type)
     with Inline
   is
      --x    Real_t  e_cut = domain.e_cut() ;
      --x    Real_t  p_cut = domain.p_cut() ;
      --x    Real_t  ss4o3 = domain.ss4o3() ;
      --x    Real_t  q_cut = domain.q_cut() ;
--        e_cut : constant Energy    := domain.parameters.energy_tolerance;
--        p_cut : constant Pressure  := domain.parameters.pressure_tolerance;
--        ss4o3 : constant Real_Type := domain.parameters.four_thirds;
--        q_cut : constant Pressure  := domain.parameters.pressure_dynamic_tolerance;
      --x    Real_t eosvmax = domain.eosvmax() ;
      --x    Real_t eosvmin = domain.eosvmin() ;
      --x    Real_t pmin    = domain.pmin() ;
      --x    Real_t emin    = domain.emin() ;
      --x    Real_t rho0    = domain.refdens() ;
--        eosvmax : constant Volume   := domain.parameters.eosvmax;
--        eosvmin : constant Volume   := domain.parameters.eosvmin;
--        pmin    : constant Pressure := domain.parameters.pressure_floor;
--        emin    : constant Energy   := domain.parameters.energy_floor;
--        rho0    : constant Density  := domain.parameters.reference_density;
      ---    // These temporaries will be of different size for
      ---    // each call (due to different sized region element
      ---    // lists)
      --x    Real_t *e_old = Allocate<Real_t>(numElemReg) ;
      --x    Real_t *delvc = Allocate<Real_t>(numElemReg) ;
      --x    Real_t *p_old = Allocate<Real_t>(numElemReg) ;
      --x    Real_t *q_old = Allocate<Real_t>(numElemReg) ;
      --x    Real_t *compression = Allocate<Real_t>(numElemReg) ;
      --x    Real_t *compHalfStep = Allocate<Real_t>(numElemReg) ;
      --x    Real_t *qq_old = Allocate<Real_t>(numElemReg) ;
      --x    Real_t *ql_old = Allocate<Real_t>(numElemReg) ;
      --x    Real_t *work = Allocate<Real_t>(numElemReg) ;
      --x    Real_t *p_new = Allocate<Real_t>(numElemReg) ;
      --x    Real_t *e_new = Allocate<Real_t>(numElemReg) ;
      --x    Real_t *q_new = Allocate<Real_t>(numElemReg) ;
      --x    Real_t *bvc = Allocate<Real_t>(numElemReg) ;
      --x    Real_t *pbvc = Allocate<Real_t>(numElemReg) ;
--        bvc          : Element_Real_Array_Access := new Element_Real_Array (0..numElemReg);
--        compHalfStep : Element_Compression_Array_Access := new Element_Compression_Array (0..numElemReg);
--        compression  : Element_Compression_Array_Access := new Element_Compression_Array (0..numElemReg);
--        delvc        : Element_Volume_Array_Access := new Element_Volume_Array (0..numElemReg);
--        e_new        : Element_Energy_Array_Access := new Element_Energy_Array (0..numElemReg);
--        e_old        : Element_Energy_Array_Access := new Element_Energy_Array (0..numElemReg);
--        p_new        : Element_Pressure_Array_Access := new Element_Pressure_Array (0..numElemReg);
--        p_old        : Element_Pressure_Array_Access := new Element_Pressure_Array (0..numElemReg);
--        pbvc         : Element_Real_Array_Access := new Element_Real_Array (0..numElemReg);
--        q_new        : Element_Pressure_Array_Access := new Element_Pressure_Array (0..numElemReg);
--        q_old        : Element_Pressure_Array_Access := new Element_Pressure_Array (0..numElemReg);
--        ql_old       : Element_Pressure_Array_Access := new Element_Pressure_Array (0..numElemReg);
--        qq_old       : Element_Pressure_Array_Access := new Element_Pressure_Array (0..numElemReg);
--        work         : Element_Real_Array_Access := new Element_Real_Array (0..numElemReg);

      info : EOS_Info_Record (numElemReg);
   begin
      ---    //loop to add load imbalance based on region number
      --x    for(Int_t j = 0; j < rep; j++) {
      for j in 0..rep-1 loop
         ---       /* compress data, minimal set */
         -- #pragma omp parallel
         --       {
         -- #pragma omp for nowait firstprivate(numElemReg)
         --x          for (Index_t i=0; i<numElemReg; ++i) {
         --x             Index_t elem = regElemList[i];
         --x             e_old[i] = domain.e(elem) ;
         --x             delvc[i] = domain.delv(elem) ;
         --x             p_old[i] = domain.p(elem) ;
         --x             q_old[i] = domain.q(elem) ;
         --x             qq_old[i] = domain.qq(elem) ;
         --x             ql_old[i] = domain.ql(elem) ;
         --x          }
         for region_element in regElemList'Range loop
            declare
               element : Element_Record renames
                 domain.elements (regElemList (region_element));
            begin
              info.e_old(region_element)  := element.eenergy;
              info.delvc(region_element)  := element.new_volume_relative_delta;
              info.p_old(region_element)  := element.pressure_static;
              info.q_old(region_element)  := element.pressure_dynamic;
              info.qq_old(region_element) := element.pressure_dynamic_quadratic;
              info.ql_old(region_element) := element.pressure_dynamic_linear;
            end;
         end loop;

         -- #pragma omp for firstprivate(numElemReg)
         --x          for (Index_t i = 0; i < numElemReg ; ++i) {
         --x             Index_t elem = regElemList[i];
         --x             Real_t vchalf ;
         --x             compression[i] = Real_t(1.) / vnewc[elem] - Real_t(1.);
         --x             vchalf = vnewc[elem] - delvc[i] * Real_t(.5);
         --x             compHalfStep[i] = Real_t(1.) / vchalf - Real_t(1.);
         --x          }
         for region_element in regElemList'Range loop
            declare
               element : constant Element_Index := regElemList (region_element);
               vchalf  : Volume;
               function To_Compression (this : in Volume) return Compression_Type is
                  (Compression_Type(1.0 / vnewc(element) - 1.0));
            begin
              info.compression(region_element) := To_Compression(vnewc(element));
              vchalf := vnewc(element) - info.delvc(region_element) * 0.5;
              info.compHalfStep(region_element) := To_Compression(vchalf);
            end;
         end loop;
         ---       /* Check for v > eosvmax or v < eosvmin */
         --x          if ( eosvmin != Real_t(0.) ) {
         -- #pragma omp for nowait firstprivate(numElemReg, eosvmin)
         --x             for(Index_t i=0 ; i<numElemReg ; ++i) {
         --x                Index_t elem = regElemList[i];
         --x                if (vnewc[elem] <= eosvmin) { /* impossible due to calling func? */
         --x                   compHalfStep[i] = compression[i] ;
         --x                }
         --x             }
         --x          }
         --x          if ( eosvmax != Real_t(0.) ) {
         -- #pragma omp for nowait firstprivate(numElemReg, eosvmax)
         --x             for(Index_t i=0 ; i<numElemReg ; ++i) {
         --x                Index_t elem = regElemList[i];
         --x                if (vnewc[elem] >= eosvmax) { /* impossible due to calling func? */
         --x                   p_old[i]        = Real_t(0.) ;
         --x                   compression[i]  = Real_t(0.) ;
         --x                   compHalfStep[i] = Real_t(0.) ;
         --x                }
         --x             }
         --x          }
         if info.eosvmin /= 0.0 then
            for region_element in regElemList'Range loop
               if vnewc(regElemList(region_element))  <= info.eosvmin then ---/* impossible due to calling func? */
                    info.compHalfStep(region_element) := info.compression(region_element);
               end if;
            end loop;
         end if;
         if info.eosvmax /= 0.0 then
            for region_element in regElemList'Range loop
               if vnewc(regElemList(region_element)) >= info.eosvmax then ---/* impossible due to calling func? */
                  info.p_old(region_element)        := 0.0;
                  info.compression(region_element)  := 0.0;
                  info.compHalfStep(region_element) := 0.0;
               end if;
            end loop;
         end if;

         -- #pragma omp for nowait firstprivate(numElemReg)
         --x          for (Index_t i = 0 ; i < numElemReg ; ++i) {
         --x             work[i] = Real_t(0.) ;
         --x          }
         --       }
         --x       CalcEnergyForElems(p_new, e_new, q_new, bvc, pbvc,
         --x                          p_old, e_old,  q_old, compression, compHalfStep,
         --x                          vnewc, work,  delvc, pmin,
         --x                          p_cut, e_cut, q_cut, emin,
         --x                          qq_old, ql_old, rho0, eosvmax,
         --x                          numElemReg, regElemList);
         info.work.all := (others =>0.0);
         CalcEnergyForElems (domain, vnewc, info, regElemList);
         --x    }
      end loop;
      -- #pragma omp parallel for firstprivate(numElemReg)
      --x    for (Index_t i=0; i<numElemReg; ++i) {
      --x       Index_t elem = regElemList[i];
      --x       domain.p(elem) = p_new[i] ;
      --x       domain.e(elem) = e_new[i] ;
      --x       domain.q(elem) = q_new[i] ;
      --x    }
      for region_element in regElemList'Range loop
         declare
            element : Element_Record renames
              domain.elements (regElemList (region_element));
         begin
            element.pressure_static  := info.p_new(region_element);
            element.eenergy          := info.e_new(region_element);
            element.pressure_dynamic := info.q_new(region_element);
         end;
      end loop;
      --x    CalcSoundSpeedForElems(domain,
      --x                           vnewc, rho0, e_new, p_new,
      --x                           pbvc, bvc, ss4o3,
      --x                           numElemReg, regElemList) ;
      CalcSoundSpeedForElems(domain,
                             vnewc, info.rho0, info.e_new, info.p_new,
                             info.pbvc, info.bvc, info.ss4o3,
                            regElemList);
      --x    Release(&pbvc) ;
      --x    Release(&bvc) ;
      --x    Release(&q_new) ;
      --x    Release(&e_new) ;
      --x    Release(&p_new) ;
      --x    Release(&work) ;
      --x    Release(&ql_old) ;
      --x    Release(&qq_old) ;
      --x    Release(&compHalfStep) ;
      --x    Release(&compression) ;
      --x    Release(&q_old) ;
      --x    Release(&p_old) ;
      --x    Release(&delvc) ;
      --x    Release(&e_old) ;
      pragma Warnings (Off, "* modified by call, but value never referenced");
      Release (info.bvc);
      Release (info.compHalfStep);
      Release (info.compression);
      Release (info.delvc);
      Release (info.e_new);
      Release (info.e_old);
      Release (info.p_new);
      Release (info.p_old);
      Release (info.pbvc);
      Release (info.q_new);
      Release (info.q_old);
      Release (info.ql_old);
      Release (info.qq_old);
      Release (info.work);
      pragma Warnings (On, "* modified by call, but value never referenced");
      --x }
   end EvalEOSForElems;

   --- /******************************************/

   --x static inline
   --x void ApplyMaterialPropertiesForElems(Domain& domain, Real_t vnew[])
   --x {
   procedure ApplyMaterialPropertiesForElems
     (domain : in out Domain_Record;
      vnew   : access Element_Volume_Array)
     with Inline
   is
      --x    Index_t numElem = domain.numElem() ;
      numElem : constant Element_Count := domain.numElem;
   begin
      --x   if (numElem != 0) {
      if numElem /=0 then
         declare
            ---     /* Expose all of the variables needed for material evaluation */
            --x     Real_t eosvmin = domain.eosvmin() ;
            --x     Real_t eosvmax = domain.eosvmax() ;
            eosvmin : constant Volume := domain.parameters.eosvmin;
            eosvmax : constant Volume := domain.parameters.eosvmax;
         begin
            -- #pragma omp parallel
            --     {
            ---        // Bound the updated relative volumes with eosvmin/max
            --x        if (eosvmin != Real_t(0.)) {
            -- #pragma omp for firstprivate(numElem)
            --x           for(Index_t i=0 ; i<numElem ; ++i) {
            --x              if (vnew[i] < eosvmin)
            --x                 vnew[i] = eosvmin ;
            --x           }
            --x        }
            if eosvmin /= 0.0 then
               for element in 0..numElem-1 loop
                  if vnew(element) < eosvmin then
                     vnew(element) := eosvmin;
                  end if;
               end loop;
            end if;
            --x        if (eosvmax != Real_t(0.)) {
            -- #pragma omp for nowait firstprivate(numElem)
            --x           for(Index_t i=0 ; i<numElem ; ++i) {
            --x              if (vnew[i] > eosvmax)
            --x                 vnew[i] = eosvmax ;
            --x           }
            --x        }
            if eosvmax /= 0.0 then
               for element in 0..numElem-1 loop
                  if vnew(element) > eosvmax then
                     vnew(element) := eosvmax;
                  end if;
               end loop;
            end if;
            ---        // This check may not make perfect sense in LULESH, but
            ---        // it's representative of something in the full code -
            ---        // just leave it in, please
            -- #pragma omp for nowait firstprivate(numElem)
            --x        for (Index_t i=0; i<numElem; ++i) {
            --x           Real_t vc = domain.v(i) ;
            --x           if (eosvmin != Real_t(0.)) {
            --x              if (vc < eosvmin)
            --x                 vc = eosvmin ;
            --x           }
            --x           if (eosvmax != Real_t(0.)) {
            --x              if (vc > eosvmax)
            --x                 vc = eosvmax ;
            --x           }
            --x           if (vc <= 0.) {
            --x #if USE_MPI
            --x              MPI_Abort(MPI_COMM_WORLD, VolumeError) ;
            --x #else
            --x              exit(VolumeError);
            --x #endif
            --x           }
            --x        }
            for element in 0..numElem-1 loop
               declare
                  vc : Volume := domain.elements(element).volume_relative;
               begin
                  if eosvmin /= 0.0 then
                     if vc < eosvmin then
                        vc := eosvmin;
                     end if;
                  end if;
                  if eosvmax /= 0.0 then
                     if vc > eosvmax then
                        vc := eosvmax;
                     end if;
                  end if;
                  if vc <= 0.0 then
                     if USE_MPI then
                        MPI.Abortt (MPI.COMM_WORLD, VolumeError);
                     else
                        raise VolumeError;
                     end if;
                  end if;
                  end;
            end loop;
            --     }
            --x     for (Int_t r=0 ; r<domain.numReg() ; r++) {
            --x        Index_t numElemReg = domain.regElemSize(r);
            --x        Index_t *regElemList = domain.regElemlist(r);
            --x        Int_t rep;
            ---        //Determine load imbalance for this region
            ---        //round down the number with lowest cost
            --x        if(r < domain.numReg()/2)
            --x 	rep = 1;
            ---        //you don't get an expensive region unless you at least have 5 regions
            --x        else if(r < (domain.numReg() - (domain.numReg()+15)/20))
            --x          rep = 1 + domain.cost();
            ---        //very expensive regions
            --x        else
            --x 	rep = 10 * (1+ domain.cost());
            --x        EvalEOSForElems(domain, vnew, numElemReg, regElemList, rep);
            --x     }
            for region in 0..domain.numReg-1 loop
               declare
                  numElemReg : constant Element_Count :=
                    domain.regions(region).size;
                  regElemList : constant Element_Element_Index_Array_Access :=
                    domain.regions(region).elements;
                  rep         : Cost_Type;
               begin
                  ---        //Determine load imbalance for this region
                  ---        //round down the number with lowest cost
                  if region < domain.numReg/2 then
                     ---        //you don't get an expensive region unless you at least have 5 regions
                     rep := 1;
                  elsif region < (domain.numReg - (domain.numReg+15)/20) then
                     rep := 1 + domain.parameters.imbalance_cost;
                     ---        //very expensive regions
                  else
                     rep := 10 * (1 + domain.parameters.imbalance_cost);
                     EvalEOSForElems(domain, vnew, numElemReg, regElemList, rep);
                  end if;
               end;
            end loop;
            --x   }
         end;
      end if;
      --x }
   end ApplyMaterialPropertiesForElems;

   --- /******************************************/

   --x static inline
   --x void UpdateVolumesForElems(Domain &domain, Real_t *vnew,
   --x                            Real_t v_cut, Index_t length)
   --x {
   procedure UpdateVolumesForElems
     (domain : in out Domain_Record;
      vnew   : access Element_Volume_Array)
     with Inline is
   begin
      --x    if (length != 0) {
      -- #pragma omp parallel for firstprivate(length, v_cut)
      --x       for(Index_t i=0 ; i<length ; ++i) {
      --x          Real_t tmpV = vnew[i] ;
      --x          if ( FABS(tmpV - Real_t(1.0)) < v_cut )
      --x             tmpV = Real_t(1.0) ;
      --x          domain.v(i) = tmpV ;
      --x       }
      --x    }
      for element in vnew'Range loop
         declare
            tmpV : Volume := vnew (element);
         begin
            if abs (tmpV - 1.0) < domain.parameters.volume_relative_tolerance then
               tmpV := 1.0;
            end if;
            domain.elements (element).volume_relative := tmpV;
         end;
      end loop;
      --x    return ;
      --x }
   end UpdateVolumesForElems;

   --- /******************************************/

   --x static inline
   --x void LagrangeElements(Domain& domain, Index_t numElem)
   --x {
   --x   Real_t *vnew = Allocate<Real_t>(numElem) ;  /* new relative vol -- temp */
   procedure LagrangeElements (domain : in out domain_record)
     with inline
   is
      vnew : Element_Volume_Array_Access;
   begin
      vnew := new Element_Volume_Array (0..domain.numElem-1);
      --x   CalcLagrangeElements(domain, vnew) ;
      --x   /* Calculate Q.  (Monotonic q option requires communication) */
      --x   CalcQForElems(domain, vnew) ;
      --x   ApplyMaterialPropertiesForElems(domain, vnew) ;
      --x   UpdateVolumesForElems(domain, vnew,
      --x                         domain.v_cut(), numElem) ;
      CalcLagrangeElements (domain, vnew);
      CalcQForElems (domain, vnew);
      ApplyMaterialPropertiesForElems (domain, vnew);
      UpdateVolumesForElems (domain, vnew);
      --x   Release(&vnew);
      pragma Warnings (Off, "* modified by call, but value never referenced");
      Release (vnew);
      pragma Warnings (On, "* modified by call, but value never referenced");
      --x }
   end LagrangeElements;

   --- /******************************************/

   --x static inline
   --x void CalcCourantConstraintForElems(Domain &domain, Index_t length,
   --x                                    Index_t *regElemlist,
   --x                                    Real_t qqc, Real_t& dtcourant)
   --x {
   procedure CalcCourantConstraintForElems
     (domain : in out Domain_Record;
      region : in     Region_Index)
     with inline
   is
      -- #if _OPENMP
      --    Index_t threads = omp_get_max_threads();
      --    static Index_t *courant_elem_per_thread;
      --    static Real_t *dtcourant_per_thread;
      --    static bool first = true;
      --    if (first) {
      --      courant_elem_per_thread = new Index_t[threads];
      --      dtcourant_per_thread = new Real_t[threads];
      --      first = false;
      --    }
      -- #else
      --x    Index_t threads = 1;
      --x    Index_t courant_elem_per_thread[1];
      --x    Real_t  dtcourant_per_thread[1];
      threads                 : constant Thread_Index := 1;
      courant_elem_per_thread : Thread_Element_Count_Array (0..threads-1);
      dtcourant_per_thread    : Thread_Time_Span_Array (0..threads-1);
      -- #endif

      -- #pragma omp parallel firstprivate(length, qqc)
      --    {
      --x       Real_t   qqc2 = Real_t(64.0) * qqc * qqc ;
      --x       Real_t   dtcourant_tmp = dtcourant;
      --x       Index_t  courant_elem  = -1 ;
      qqc2          : constant Real_Type := 64.0 * domain.parameters.qqc**2;
      dtcourant_tmp : Time_Span := domain.variables.dtcourant;
      -- Picking an unlikely number instead of -1:
      NEVER_SET     : constant Element_Count := 2_000_000_003;
      courant_elem  : Element_Count := NEVER_SET;

      -- #if _OPENMP
      --       Index_t thread_num = omp_get_thread_num();
      -- #else
      --x       Index_t thread_num = 0;
      thread_num : constant Thread_Index := 0;
      -- #endif
      use type ART.Time_Span;
   begin
      -- #pragma omp for
      --x       for (Index_t i = 0 ; i < length ; ++i) {
      --x          Index_t indx = regElemlist[i] ;
      --x          Real_t dtf = domain.ss(indx) * domain.ss(indx) ;
      --x          if ( domain.vdov(indx) < Real_t(0.) ) {
      --x             dtf = dtf
      --x                 + qqc2 * domain.arealg(indx) * domain.arealg(indx)
      --x                 * domain.vdov(indx) * domain.vdov(indx) ;
      --x          }
      --x          dtf = SQRT(dtf) ;
      --x          dtf = domain.arealg(indx) / dtf ;
      for region_element in domain.regions(region).elements'Range loop
         declare
            element : constant Element_Index :=
              domain.regions(region).elements(region_element);
            dtf     : Time_Span :=
              ART.Time_Span_Last;
            ss      : constant Real_Type := Real_Type
              (domain.elements(element).sound_speed);
            vdov    : constant Real_Type := Real_Type
              (domain.elements(element).volume_derivative_over_volume);
            arealg  : constant Real_Type := Real_Type
              (domain.elements(element).characteristic_length);
         begin
            if vdov < 0.0 then
               dtf := To_Time_Span
                 (arealg / SQRT(ss**2 + qqc2 * arealg**2 * vdov**2));
            else
               dtf := To_Time_Span
                 (arealg / abs(ss));
            end if;
            --x          if (domain.vdov(indx) != Real_t(0.)) {
            --x             if ( dtf < dtcourant_tmp ) {
            --x                dtcourant_tmp = dtf ;
            --x                courant_elem  = indx ;
            --x             }
            --x          }
            --x       }
            if vdov /= 0.0 and then dtcourant_tmp > dtf then
               dtcourant_tmp := dtf;
               courant_elem  := element;
            end if;
         end;
      end loop;
      --x       dtcourant_per_thread[thread_num]    = dtcourant_tmp ;
      --x       courant_elem_per_thread[thread_num] = courant_elem ;
      --    }
      dtcourant_per_thread(thread_num)    := dtcourant_tmp;
      courant_elem_per_thread(thread_num) := courant_elem;
      --x    for (Index_t i = 1; i < threads; ++i) {
      --x       if (dtcourant_per_thread[i] < dtcourant_per_thread[0] ) {
      --x          dtcourant_per_thread[0]    = dtcourant_per_thread[i];
      --x          courant_elem_per_thread[0] = courant_elem_per_thread[i];
      --x       }
      --x    }
      for thread in 0..threads loop
         if dtcourant_per_thread(thread) < dtcourant_per_thread(0) then
            dtcourant_per_thread(0)    := dtcourant_per_thread(thread);
            courant_elem_per_thread(0) := courant_elem_per_thread(thread);
         end if;
      end loop;
      --x    if (courant_elem_per_thread[0] != -1) {
      --x       dtcourant = dtcourant_per_thread[0] ;
      --x    }
      if courant_elem_per_thread(0) /= NEVER_SET then
        domain.variables.dtcourant := dtcourant_per_thread(0);
      end if;
      --x    return ;
      --x }
   end CalcCourantConstraintForElems;

   --x /******************************************/

   --x static inline
   --x void ((Domain &domain, Index_t length,
   --x                                  Index_t *regElemlist, Real_t dvovmax, Real_t& dthydro)
   --x {
   procedure CalcHydroConstraintForElems
     (domain : in out Domain_Record;
      region : in     Region_Index)
     with inline
   is
      -- #if _OPENMP
      --    Index_t threads = omp_get_max_threads();
      --    static Index_t *hydro_elem_per_thread;
      --    static Real_t *dthydro_per_thread;
      --    static bool first = true;
      --    if (first) {
      --      hydro_elem_per_thread = new Index_t[threads];
      --      dthydro_per_thread = new Real_t[threads];
      --      first = false;
      --    }
      -- #else
      --x    Index_t threads = 1;
      --x    Index_t hydro_elem_per_thread[1];
      --x    Real_t  dthydro_per_thread[1];
      threads               : constant Thread_Index := 1;
      hydro_elem_per_thread : Thread_Element_Count_Array (0..threads-1);
      dthydro_per_thread    : Thread_Time_Span_Array (0..threads-1);
      -- #endif
      -- #pragma omp parallel firstprivate(length, dvovmax)
      --    {
      --x       Real_t dthydro_tmp = dthydro ;
      --x       Index_t hydro_elem = -1 ;
      dthydro_tmp : Time_Span := domain.variables.dthydro;
      --- Picking an unlikely number instead of -1:
      NEVER_SET   : constant Element_Count := 2_000_000_001;
      hydro_elem  : Element_Count := NEVER_SET;
      -- #if _OPENMP
      --       Index_t thread_num = omp_get_thread_num();
      -- #else
      --x       Index_t thread_num = 0;
      thread_num : constant Thread_Index := 0;
      -- #endif
      use type ART.Time_Span;
   begin
      -- #pragma omp for
      --x       for (Index_t i = 0 ; i < length ; ++i) {
      --x          Index_t indx = regElemlist[i] ;
      --x          if (domain.vdov(indx) != Real_t(0.)) {
      --x             Real_t dtdvov = dvovmax / (FABS(domain.vdov(indx))+Real_t(1.e-20)) ;
      --x             if ( dthydro_tmp > dtdvov ) {
      --x                   dthydro_tmp = dtdvov ;
      --x                   hydro_elem = indx ;
      --x             }
      --x          }
      --x       }
      for region_element in domain.regions(region).elements'Range loop
         declare
            element : constant Element_Index :=
              domain.regions(region).elements(region_element);
            vdov : VDOV_Type renames
              domain.elements(element).volume_derivative_over_volume;
         begin
            if vdov /= 0.0 then
               declare
                  --- dvovmax/vdov=(dv/vmax)/(vd/v)=(dv*v)/(vmax*vd)=dv/vd
                  dtdvov : constant Time_Span :=
                    domain.parameters.dvovmax / (abs(vdov)+1.0e-20);
               begin
                  if dthydro_tmp > dtdvov then
                     dthydro_tmp := dtdvov;
                     hydro_elem := element;
                  end if;
               end;
            end if;
         end;
      end loop;
      --x       dthydro_per_thread[thread_num]    = dthydro_tmp ;
      --x       hydro_elem_per_thread[thread_num] = hydro_elem ;
      --    }
      dthydro_per_thread(thread_num)    := dthydro_tmp;
      hydro_elem_per_thread(thread_num) := hydro_elem;
      --x    for (Index_t i = 1; i < threads; ++i) {
      --x       if(dthydro_per_thread[i] < dthydro_per_thread[0]) {
      --x          dthydro_per_thread[0]    = dthydro_per_thread[i];
      --x          hydro_elem_per_thread[0] =  hydro_elem_per_thread[i];
      --x       }
      --x    }
      for thread in 0..threads loop
         if dthydro_per_thread(thread) < dthydro_per_thread(0) then
            dthydro_per_thread(0)    := dthydro_per_thread(thread);
            hydro_elem_per_thread(0) := hydro_elem_per_thread(thread);
         end if;
      end loop;
      --x    if (hydro_elem_per_thread[0] != -1) {
      --x       dthydro =  dthydro_per_thread[0] ;
      --x    }
      if hydro_elem_per_thread(0) /= NEVER_SET then
        domain.variables.dthydro := dthydro_per_thread(0);
      end if;
   --x    return ;
   --x }
   end CalcHydroConstraintForElems;

   --- /******************************************/

   --x static inline
   --x void CalcTimeConstraintsForElems(Domain& domain) {
   procedure CalcTimeConstraintsForElems (domain : in out Domain_Record)
     with inline is
   begin
      ---    // Initialize conditions to a very large value
      --x    domain.dtcourant() = 1.0e+20;
      --x    domain.dthydro() = 1.0e+20;
      domain.variables.dtcourant := ART.Time_Span_Last;
      domain.variables.dthydro := ART.Time_Span_Last;
      --x    for (Index_t r=0 ; r < domain.numReg() ; ++r) {
      --x       /* evaluate time constraint */
      --x       CalcCourantConstraintForElems(domain, domain.regElemSize(r),
      --x                                     domain.regElemlist(r),
      --x                                     domain.qqc(),
      --x                                     domain.dtcourant()) ;
      --x       /* check hydro constraint */
      --x       CalcHydroConstraintForElems(domain, domain.regElemSize(r),
      --x                                   domain.regElemlist(r),
      --x                                   domain.dvovmax(),
      --x                                   domain.dthydro()) ;
      --x    }
      for region in 0..domain.numReg-1 loop
         CalcCourantConstraintForElems(domain, region);
         CalcHydroConstraintForElems(domain, region);
      end loop;
      --x }
   end CalcTimeConstraintsForElems;

   --- /******************************************/

   --x static inline
   --x void LagrangeLeapFrog(Domain& domain)
   ----------------------
   -- EXPORTED (private):
   ----------------------
   procedure LagrangeLeapFrog
     (domain : in out Domain_Record) is
      --x {
   begin
      -- #ifdef SEDOV_SYNC_POS_VEL_LATE
      --    Domain_member fieldData[6] ;
      -- #endif

      ---    /* calculate nodal forces, accelerations, velocities, positions, with
      ---     * applied boundary conditions and slide surface considerations */
      --x    LagrangeNodal(domain);
      LagrangeNodal(domain);


      -- #ifdef SEDOV_SYNC_POS_VEL_LATE
      -- #endif

      ---    /* calculate element quantities (i.e. velocity gradient & q), and update
      ---     * material states */
      --x    LagrangeElements(domain, domain.numElem());
      LagrangeElements(domain);

      -- #if USE_MPI
      -- #ifdef SEDOV_SYNC_POS_VEL_LATE
      --    CommRecv(domain, MSG_SYNC_POS_VEL, 6,
      --             domain.sizeX() + 1, domain.sizeY() + 1, domain.sizeZ() + 1,
      --             false, false) ;

      --    fieldData[0] = &Domain::x ;
      --    fieldData[1] = &Domain::y ;
      --    fieldData[2] = &Domain::z ;
      --    fieldData[3] = &Domain::xd ;
      --    fieldData[4] = &Domain::yd ;
      --    fieldData[5] = &Domain::zd ;

      --    CommSend(domain, MSG_SYNC_POS_VEL, 6, fieldData,
      --             domain.sizeX() + 1, domain.sizeY() + 1, domain.sizeZ() + 1,
      --             false, false) ;
      -- #endif
      -- #endif

      --x    CalcTimeConstraintsForElems(domain);
      CalcTimeConstraintsForElems(domain);

      --x #if USE_MPI
      --x #ifdef SEDOV_SYNC_POS_VEL_LATE
      --x    CommSyncPosVel(domain) ;
      --x #endif
      --x #endif
      if USE_MPI and then SEDOV_SYNC_POS_VEL_LATE then
        Comm.SyncPosVel (domain);
      end if;
      --x }
   end LagrangeLeapFrog;

end LULESH;
