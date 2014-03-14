-- This spec is just a stand-in for a real Ada to C MPI binding, to enable
-- compilation of Ada translated from MPI-using C code.
with Ada.Real_Time;
with System;

package MPI is

   package ART renames Ada.Real_Time;

   -- MPI_comm - communicator (handle):
   type Comm_Type is private;
   COMM_WORLD : constant Comm_Type;

   type errorcode_Type is new Integer;
   type Datatype_Type is
     (DOUBLE);

   type Reduce_Op_Type is
     (MAX);

   subtype Rank_Type is Integer;

   -- Environment Management Routines:

   -- int MPI_Init(int *argc, char ***argv)
   procedure Init;

   --  int MPI_Comm_size( MPI_Comm comm, int *size )
   procedure Comm_size
     (comm : in  Comm_Type;
      size : out Rank_Type);

   -- int MPI_Comm_rank( MPI_Comm comm, int *rank )
   procedure Comm_rank
     (comm : in  Comm_Type;
      rank : out Rank_Type);

   -- int MPI_Abort( MPI_Comm comm, int errorcode )
   procedure Abortt
     (comm      : in  Comm_Type;
      errorcode : in errorcode_Type);

   -- double MPI_Wtime()
   function Wtime return ART.Time;

   -- int MPI_Finalize()
   procedure Finalize;

   -- Point to Point Communication Routines:

   -- Collective Communication Routines:

   -- int MPI_Barrier( MPI_Comm comm )
   procedure Barrier
     (comm : in  Comm_Type);

   -- int MPI_Reduce ( void *sendbuf, void *recvbuf, int count,
   --                   MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm )
   procedure Reduce
     (sendbuf  : in System.Address;
      recvbuf  : in System.Address;
      count    : in Natural;
      datatype : in Datatype_Type;
      op       : in Reduce_Op_Type;
      root     : in Integer;
      comm     : in Comm_Type);



private
   type Comm_Type is new Integer;
   COMM_WORLD : constant Comm_Type := 0;
end MPI;
