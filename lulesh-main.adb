--- /******************************************/
with Ada.Text_IO;
with LULESH.Comm;
with LULESH.Init;
with LULESH.Util;
with LULESH.Viz;
with MPI;

--x int main(int argc, char *argv[])
--x {
procedure Lulesh.main is
   package ATI renames Ada.Text_IO;

   --x    Domain *locDom ;
   --x    Int_t numRanks ;
   --x    Int_t myRank ;
   --x    struct cmdLineOpts opts;
   locDom   : Domain_Record;
   numRanks : Rank_Type;
   myRank   : Rank_Type;
   opts     : LULESH.Util.cmdLineOpts;

   --- from declarations below:
   start         : ART.Time;
   endd          : ART.Time;
   elapsed_time  : ART.Time_Span;
   elapsed_timeG : ART.Time_Span;
   fieldData     : Comm.Domain_member;
   use type ART.Time;

   function Image (this : in ART.Time) return String
   is
      SC : ART.Seconds_Count;
      TS : ART.Time_Span;
   begin
      ART.Split (this, SC, TS);
      return SC'Img & " " & ART.To_Duration (TS)'Img;
   end Image;

   function Image (this : in ART.Time_Span) return String is
     (ART.To_Duration(this)'Img);

begin
   --x #if USE_MPI
   --x    Domain_member fieldData ;
   --x    MPI_Init(&argc, &argv) ;
   --x    MPI_Comm_size(MPI_COMM_WORLD, &numRanks) ;
   --x    MPI_Comm_rank(MPI_COMM_WORLD, &myRank) ;
   --x #else
   --x    numRanks = 1;
   --x    myRank = 0;
   --x #endif
   if USE_MPI then
     MPI.Init;
     MPI.Comm_size (MPI.COMM_WORLD, numRanks);
     MPI.Comm_rank (MPI.COMM_WORLD, myRank);
   else
      numRanks := 1;
      myRank   := 0;
   end if;
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
     (its         => 9999999,
      side_length => 30,
      numReg      => 11,
      numFiles    => (Integer (numRanks)+10) / 9,
      showProg    => False,
      quiet       => False,
      viz         => False,
      Balance     => 1,
      Cost        => 1);
   --x    ParseCommandLineOptions(argc, argv, myRank, &opts);
   LULESH.Util.ParseCommandLineOptions(myRank, opts);
   --x    if ((myRank == 0) && (opts.quiet == 0)) {
   --x       printf("Running problem size %d^3 per domain until completion\n", opts.nx);
   --x       printf("Num processors: %d\n", numRanks);
   --x #if _OPENMP
   --x       printf("Num threads: %d\n", omp_get_max_threads());
   --x #endif
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
      ATI.Put_Line ("Running problem size " & opts.side_length'Img & "^3 per domain until completion");
      ATI.Put_Line ("Num processors: " & numRanks'Img);
--        if USE_OMP then
--           ATI.Put_Line ("Num threads: " & omp_get_max_threads'Img);
--        end if;
      ATI.Put_Line ("Total number of elements: " &
                      Integer'Image(Integer(numRanks)*Integer(opts.side_length)**3));
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
      domain_col       : Domain_Index;
      domain_row       : Domain_Index;
      domain_plane     : Domain_Index;
      domains_per_side : Domain_Index;
   begin
      --x    InitMeshDecomp(numRanks, myRank, &col, &row, &domain_, &side);
      LULESH.Init.InitMeshDecomp
        (numRanks,
         myRank,
         domain_col,
         domain_row,
         domain_plane,
         domains_per_side);
      ---    // Build the main data structure and initialize it
      --x    locDom = new Domain(numRanks, col, row, plane, opts.nx,
      --x                        side, opts.numReg, opts.balance, opts.cost) ;
      locDom := LULESH.Init.Create
        (numRanks    => numRanks,
         colLoc      => domain_col,
         rowLoc      => domain_row,
         planeLoc    => domain_plane,
         side_length => opts.side_length,
         tp          => domains_per_side,
         nr          => opts.numReg,
         balance     => opts.balance,
         cost        => opts.cost);
   end;
   --x #if USE_MPI
   --x    fieldData = &Domain::nodalMass ;
   --x    // Initial domain boundary communication
   --x    CommRecv(*locDom, MSG_COMM_SBN, 1,
   --x             locDom->sizeX() + 1, locDom->sizeY() + 1, locDom->sizeZ() + 1,
   --x             true, false) ;
   --x    CommSend(*locDom, MSG_COMM_SBN, 1, &fieldData,
   --x             locDom->sizeX() + 1, locDom->sizeY() + 1, locDom->sizeZ() +  1,
   --x             true, false) ;
   --x    CommSBN(*locDom, 1, &fieldData) ;
   --x    // End initialization
   --x    MPI_Barrier(MPI_COMM_WORLD);
   --x #endif
   if USE_MPI then
      --!! pointer to the array of all masses??
     fieldData := Domain.nodalMass ;
     --- // Initial domain boundary communication
      LULESH.Comm.Recv
        (domain     => locDom,
         msgType    => MSG_COMM_SBN,
         xferFields => 1,
         dx         => locDom.parameters.size(X) + 1,
         dy         => locDom.parameters.size(Y) + 1,
         dz         => locDom.parameters.size(Z) + 1,
         doRecv     => true,
         planeOnly  => false);
      LULESH.Comm.Send
        (domain     => locDom,
         msgType    => MSG_COMM_SBN,
         xferFields => 1,
         fieldData  => (0 => FieldData),
         dx         => locDom.parameters.size(X) + 1,
         dy         => locDom.parameters.size(Y) + 1,
         dz         => locDom.parameters.size(Z) + 1,
         doSend     => true,
         planeOnly  => false);
      LULESH.Comm.SBN
        (domain     => locDom,
         xferFields => 1,
         fieldData  => fieldData);
     --- // End initialization
     MPI.Barrier (MPI.COMM_WORLD);
   end if;
   ---    // BEGIN timestep to solution */
   --x #if USE_MPI
   --x    start = MPI_Wtime();
   --x #else
   --x    timeval start;
   --x    gettimeofday(&start, NULL) ;
   --x #endif
   if USE_MPI then
      start := MPI.Wtime;
   else
      start := ART.Clock;
   end if;
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
   while locDom.variables.current_time < locDom.variables.stoptime and
     locDom.variables.cycle < opts.its loop

      TimeIncrement (locDom);
      LagrangeLeapFrog (locDom);

      if opts.showProg and not opts.quiet and myRank = 0 then
         ATI.Put ("cycle = " & locDom.variables.cycle'Img
                  & ", current_time = ");
         MKSIO.Put (locDom.variables.current_time);
         ATI.Put (", dt = ");
         MKSIO.Put (locDom.variables.deltatime); ATI.New_Line;
      end if;
   end loop;
   ---    // Use reduced max elapsed time
   --x    double elapsed_time;
   --x #if USE_MPI
   --x    elapsed_time = MPI_Wtime() - start;
   --x #else
   --x    timeval end;
   --x    gettimeofday(&end, NULL) ;
   --x    elapsed_time = (double)(end.tv_sec - start.tv_sec) + ((double)(end.tv_usec - start.tv_usec))/1000000 ;
   --x #endif
   if USE_MPI then
      elapsed_time := MPI.Wtime - start;
   else
      endd := ART.Clock;
      elapsed_time := endd - start;
   end if;
   --x    double elapsed_timeG;
   --x #if USE_MPI
   --x    MPI_Reduce(&elapsed_time, &elapsed_timeG, 1, MPI_DOUBLE,
   --x               MPI_MAX, 0, MPI_COMM_WORLD);
   --x #else
   --x    elapsed_timeG = elapsed_time;
   --x #endif
   if USE_MPI then
      MPI.Reduce
        (sendbuf  => elapsed_time'Address,
         recvbuf  => elapsed_timeG'Address,
         count    => 1,
         datatype => MPI.DOUBLE,
         op       => MPI.MAX,
         root     => 0,
         comm     => MPI.COMM_WORLD);
   else
      elapsed_timeG := elapsed_time;
   end if;
   ---    // Write out final viz file */
   --x    if (opts.viz) {
   --x       DumpToVisit(*locDom, opts.numFiles, myRank, numRanks) ;
   --x    }
   if opts.viz then
      LULESH.Viz.DumpToVisit
        (domain  => locDom,
         numfiles => opts.numFiles,
         myRank   => myRank,
         numRanks => numRanks);
   end if;
   --x    if ((myRank == 0) && (opts.quiet == 0)) {
   --x       VerifyAndWriteFinalOutput(elapsed_timeG, *locDom, opts.nx, numRanks);
   --x    }
   if myRank = 0 and not opts.quiet then
      LULESH.Util.VerifyAndWriteFinalOutput
        (elapsed_timeG, locDom, opts.side_length, numRanks);
   end if;
   --x #if USE_MPI
   --x    MPI_Finalize() ;
   --x #endif
   if USE_MPI then
     MPI.Finalize;
   end if;
   --x    return 0 ;
   --x }
end LULESH.main;
