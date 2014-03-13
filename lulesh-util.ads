package LULESH.Util is

   --x struct cmdLineOpts {
   --x    Int_t its; // -i
   --x    Int_t nx;  // -s
   --x    Int_t numReg; // -r
   --x    Int_t numFiles; // -f
   --x    Int_t showProg; // -p
   --x    Int_t quiet; // -q
   --x    Int_t viz; // -v
   --x    Int_t cost; // -c
   --x    Int_t balance; // -b
   --x };
   type cmdLineOpts is record
      its         : Int_t;
      side_length : Element_Index;
      numReg      : Region_Index;
      numFiles    : Int_t;
      showProg    : Boolean;
      quiet       : Boolean;
      viz         : Boolean;
      cost        : Cost_Type;
      balance     : Balance_Type;
   end record;

   --- // lulesh-util

   --x void ParseCommandLineOptions(int argc, char *argv[],
   --x                              Int_t myRank, struct cmdLineOpts *opts);
   procedure ParseCommandLineOptions
     (myRank : in MPI.Rank_Type;
      opts   : in out cmdLineOpts);

   --x void VerifyAndWriteFinalOutput(Real_t elapsed_time,
   --x                                Domain& locDom,
   --x                                Int_t nx,
   --x                                Int_t numRanks);
   procedure VerifyAndWriteFinalOutput
     (elapsed_time : in ART.Time_Span;
      locDom       : in out Domain_Record;
      side_length  : in Element_Index;
      numRanks     : in MPI.Rank_Type);


end LULESH.Util;
