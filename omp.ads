-- This spec is just a stand-in for a real Ada OMP binding, to enable
-- compilation of Ada translated from OMP-using C code.
package OMP is

   type Thread_Index is new Natural;
   subtype Thread_Count is Thread_Index range 1.. Thread_Index'Last;

   -- https://computing.llnl.gov/tutorials/openMP/#OMP_GET_MAX_THREADS
   -- int omp_get_max_threads(void)
   function Get_Max_Threads return Thread_Count;

   -- https://computing.llnl.gov/tutorials/openMP/#OMP_GET_NUM_THREADS
   -- int omp_get_num_threads(void)
   function Get_Num_Threads return Thread_Count;

   -- https://computing.llnl.gov/tutorials/openMP/#OMP_GET_THREAD_NUM
   -- int omp_get_thread_num(void)
   function Get_Thread_Num return Thread_Index;

end OMP;
