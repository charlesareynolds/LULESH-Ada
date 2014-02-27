with Ada.Real_Time;
with Ada.Text_IO;
with GNAT.Source_Info;
with Unchecked_Deallocation;

procedure Test_Optimization is
   package ATI renames Ada.Text_IO;
   package ART renames Ada.Real_Time;
   function Location return string
     renames GNAT.Source_Info.Source_Location;
   type Integer_Matrix is array
     (Positive range <>, Positive range <>, Positive range <>) of Integer;
   type Integer_Matrix_Access is access Integer_Matrix;
   procedure Free is new Unchecked_Deallocation
     (Object => Integer_Matrix,
      Name   => Integer_Matrix_Access);

   procedure Test_Matrix (Size : in Positive) is
      Integers : Integer_Matrix_Access;
      Start_Time : ART.Time;
      Elapsed : Duration;
      Time_Per_Capita : Float;
      use type ART.Time;
   begin
      ATI.Put_line ("Testing" & Size'Img & " cubed matrix.");
      Start_Time := ART.Clock;
      Integers := new Integer_Matrix (1..Size, 1..Size, 1..Size);
      for X in Integers'Range(1) loop
         for Y in Integers'Range(2) loop
            for Z in Integers'Range(3) loop
               Integers (X,Y,Z) := X+Y+Z;
            end loop;
         end loop;
      end loop;
      Elapsed := ART.To_Duration(ART.Clock - Start_Time);
      ATI.Put_line ("Elapsed time:" & Elapsed'Img);
      Time_Per_Capita := Float(Elapsed) / Float (Size * Size * Size);
      ATI.Put_line ("Time per element:" & Time_Per_Capita'Img);
      Free (Integers);
   end Test_Matrix;

begin
   ATI.Put_line ("BEGIN");
   for Index in 0..20 loop
      Test_Matrix (2**Index);
   end loop;
   ATI.Put_line ("END");
end Test_Optimization;
