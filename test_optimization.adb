with Ada.Real_Time;
with Ada.Text_IO;
with Unchecked_Deallocation;

procedure Test_Optimization is
   package ATI renames Ada.Text_IO;
   package ART renames Ada.Real_Time;
   type Value_Type is new Short_Short_Integer;
   type Value_Matrix is array
     (Positive range <>,
      Positive range <>,
      Positive range <>) of Value_Type;
   type Value_Matrix_Access is access Value_Matrix;
   procedure Free is new Unchecked_Deallocation
     (Object => Value_Matrix,
      Name   => Value_Matrix_Access);

   procedure Test_Matrix (Size : in Positive) is
      Values          : Value_Matrix_Access;
      Start_Time      : ART.Time;
      Elapsed         : Duration;
      Iterations      : Natural := 0;
      Time_Per_Capita : Float;
      use type ART.Time;
   begin
      ATI.Put_line ("Testing" & Size'Img & " cubed matrix.");
      Start_Time := ART.Clock;
      Values := new Value_Matrix (1..Size, 1..Size, 1..Size);
      for repeats in 1..10 loop
         for X in Values'Range(1) loop
            for Y in Values'Range(2) loop
               for Z in Values'Range(3) loop
                  Iterations := Iterations + 1;
                  Values (X,Y,Z) := Value_Type(55);
               end loop;
            end loop;
         end loop;
      end loop;
      Elapsed := ART.To_Duration(ART.Clock - Start_Time);
      ATI.Put_line ("Iterations:" & Iterations'Img);
      ATI.Put_line ("Elapsed time:" & Elapsed'Img);
      Time_Per_Capita := Float(Elapsed) / Float (Iterations);
      ATI.Put_line ("Time per element:" & Time_Per_Capita'Img);
      Free (Values);
   end Test_Matrix;

begin
   ATI.Put_line ("BEGIN");
   for Index in 0..9 loop
      Test_Matrix (2**Index);
   end loop;
   ATI.Put_line ("END");
end Test_Optimization;
