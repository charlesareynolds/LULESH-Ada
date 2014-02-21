--- /******************************************/
with Ada.Text_IO;
with LULESH.Init;
with LULESH.Util;
with LULESH.Viz;

--x int main(int argc, char *argv[])
--x {
procedure Lulesh.main is
   package ATI renames Ada.Text_IO;

   --x    Domain *locDom ;
   --x    Int_t numRanks ;
   --x    Int_t myRank ;
   --x    struct cmdLineOpts opts;
   locDom   : Access_Domain;
   numRanks : Int_t;
   myRank   : Int_t;
   opts     : LULESH.Util.cmdLineOpts;

   -- from declarations below:
   start         : AC.Time;
   endd          : AC.Time;
   elapsed_time  : AC.Day_Duration;
   elapsed_timeG : AC.Day_Duration;

   -- #if USE_MPI
   --    Domain_member fieldData ;

begin
   --    MPI_Init(&argc, &argv) ;
   --    MPI_Comm_size(MPI_COMM_WORLD, &numRanks) ;
   --    MPI_Comm_rank(MPI_COMM_WORLD, &myRank) ;
   -- #else
   --x    numRanks = 1;
   --x    myRank = 0;
   numRanks := 1;
   myRank   := 0;
   -- #endif

   ---    /* Set defaults that can be overridden by command line opts */
   --x    opts.its = 9999999;
   --x    opts.nx  = 30;
   --x    opts.numReg = 11;
   --x    opts.numFiles = (int)(numRanks+10)/9;
   --x    opts.showProg = 0;
   --x    opts.quiet = 0;
   --x    opts.viz = 0;
   --x    opts.balance = 1;
   --x    opts.cost = 1;
   opts :=
     (its      => 9999999,
      nx       => 30,
      numReg   => 11,
      numFiles => (numRanks+10)/9,
      showProg => False,
      quiet    => False,
      viz      => False,
      balance  => 1,
      cost     => 1);

   --x    ParseCommandLineOptions(argc, argv, myRank, &opts);
   LULESH.Util.ParseCommandLineOptions(myRank, opts);

   --x    if ((myRank == 0) && (opts.quiet == 0)) {
   --x       printf("Running problem size %d^3 per domain until completion\n", opts.nx);
   --x       printf("Num processors: %d\n", numRanks);
   -- #if _OPENMP
   --       printf("Num threads: %d\n", omp_get_max_threads());
   -- #endif
   --x       printf("Total number of elements: %lld\n\n", numRanks*opts.nx*opts.nx*opts.nx);
   --x       printf("To run other sizes, use -s <integer>.\n");
   --x       printf("To run a fixed number of iterations, use -i <integer>.\n");
   --x       printf("To run a more or less balanced region set, use -b <integer>.\n");
   --x       printf("To change the relative costs of regions, use -c <integer>.\n");
   --x       printf("To print out progress, use -p\n");
   --x       printf("To write an output file for VisIt, use -v\n");
   --x       printf("See help (-h) for more options\n\n");
   --x    }
   if (myRank = 0 and not opts.quiet) then
      ATI.Put_Line ("Running problem size " & opts.nx'Img & "^3 per domain until completion");
      ATI.Put_Line ("Num processors: " & numRanks'Img);
      ATI.Put_Line ("Total number of elements: " &
                      Int_t'Image(numRanks*opts.nx*opts.nx*opts.nx));
      ATI.Put_Line ("");
      ATI.Put_Line ("To run other sizes, use -s <integer>.");
      ATI.Put_Line ("To run a fixed number of iterations, use -i <integer>.");
      ATI.Put_Line ("To run a more or less balanced region set, use -b <integer>.");
      ATI.Put_Line ("To change the relative costs of regions, use -c <integer>.");
      ATI.Put_Line ("To print out progress, use -p");
      ATI.Put_Line ("To write an output file for VisIt, use -v");
      ATI.Put_Line ("See help (-h) for more options");
      ATI.Put_Line ("");
   end if;

   ---    // Set up the mesh and decompose. Assumes regular cubes for now
   --x    Int_t col, row, plane, side;
   declare
      col   : Int_t;
      row   : Int_t;
      plane : Int_t;
      side  : Int_t;
   begin
      --x    InitMeshDecomp(numRanks, myRank, &col, &row, &plane, &side);
      LULESH.Init.InitMeshDecomp
        (numRanks,
         myRank,
         col,
         row,
         plane,
         side);
      ---    // Build the main data structure and initialize it
      --x    locDom = new Domain(numRanks, col, row, plane, opts.nx,
      --x                        side, opts.numReg, opts.balance, opts.cost) ;
      locDom := LULESH.Init.Create
        (numRanks => numRanks,
         colLoc   => Index_Type (col),
         rowLoc   => Index_Type (row),
         planeLoc => Index_Type (plane),
         nx       => Index_Type (opts.nx),
         tp       => side,
         nr       => opts.numReg,
         balance  => opts.balance,
         cost     => opts.cost);
   end;


   -- #if USE_MPI
   --    fieldData = &Domain::nodalMass ;

   --    // Initial domain boundary communication
   --    CommRecv(*locDom, MSG_COMM_SBN, 1,
   --             locDom->sizeX() + 1, locDom->sizeY() + 1, locDom->sizeZ() + 1,
   --             true, false) ;
   --    CommSend(*locDom, MSG_COMM_SBN, 1, &fieldData,
   --             locDom->sizeX() + 1, locDom->sizeY() + 1, locDom->sizeZ() +  1,
   --             true, false) ;
   --    CommSBN(*locDom, 1, &fieldData) ;

   --    // End initialization
   --    MPI_Barrier(MPI_COMM_WORLD);
   -- #endif

   --    // BEGIN timestep to solution */
   -- #if USE_MPI
   --    start = MPI_Wtime();
   -- #else
   --x    timeval start;
   --x    gettimeofday(&start, NULL) ;
   -- #endif
   start := Ada.Calendar.Clock;
   --- //debug to see region sizes
   --- //   for(Int_t i = 0; i < locDom->numReg(); i++)
   --- //      std::cout << "region" << i + 1<< "size" << locDom->regElemSize(i) <<std::endl;
   ---    while((locDom->time() < locDom->stoptime()) && (locDom->cycle() < opts.its)) {

   --x       TimeIncrement(*locDom) ;
   --x       LagrangeLeapFrog(*locDom) ;

   --x       if ((opts.showProg != 0) && (opts.quiet == 0) && (myRank == 0)) {
   --x          printf("cycle = %d, time = %e, dt=%e\n",
   --x                 locDom->cycle(), double(locDom->time()), double(locDom->deltatime()) ) ;
   --x       }
   --x    }
   while locDom.variables.time < locDom.variables.stoptime and
     locDom.variables.cycle < opts.its loop

      TimeIncrement (locDom);
      LagrangeLeapFrog (locDom);

      if opts.showProg and not opts.quiet and myRank = 0 then
         ATI.Put_Line ("cycle = " & locDom.variables.cycle'Img
                       & ", time = "& locDom.variables.time'Img
                       & ", dt = " & locDom.variables.deltatime'Img);
      end if;
   end loop;

   --    // Use reduced max elapsed time
   --    double elapsed_time;
   -- #if USE_MPI
   --    elapsed_time = MPI_Wtime() - start;
   -- #else
   --x    timeval end;
   --x    gettimeofday(&end, NULL) ;
   --x    elapsed_time = (double)(end.tv_sec - start.tv_sec) + ((double)(end.tv_usec - start.tv_usec))/1000000 ;
   endd := AC.Clock;
   elapsed_time := AC."-" (endd, start);
   -- #endif
   --x    double elapsed_timeG;
   -- #if USE_MPI
   --    MPI_Reduce(&elapsed_time, &elapsed_timeG, 1, MPI_DOUBLE,
   --               MPI_MAX, 0, MPI_COMM_WORLD);
   -- #else
   --x    elapsed_timeG = elapsed_time;
   elapsed_timeG := elapsed_time;
   -- #endif

   ---    // Write out final viz file */
   --x    if (opts.viz) {
   --x       DumpToVisit(*locDom, opts.numFiles, myRank, numRanks) ;
   --x    }
   if opts.viz then
      LULESH.Viz.DumpToVisit
        (domainn  => locDom,
         numfiles => opts.numFiles,
         myRank   => myRank,
         numRanks => numRanks);
   end if;

   --x    if ((myRank == 0) && (opts.quiet == 0)) {
   --x       VerifyAndWriteFinalOutput(elapsed_timeG, *locDom, opts.nx, numRanks);
   --x    }
   if myRank = 0 and not opts.quiet then
      LULESH.Util.VerifyAndWriteFinalOutput
        (elapsed_timeG, locDom, opts.nx, numRanks);
   end if;

   -- #if USE_MPI
   --    MPI_Finalize() ;
   -- #endif

   --x    return 0 ;
   --x }
end LULESH.main;
