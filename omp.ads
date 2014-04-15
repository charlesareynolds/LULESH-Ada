package OMP is

   type Thread_Count is new Positive;

   -- https://computing.llnl.gov/tutorials/openMP/#OMP_GET_MAX_THREADS
   -- int omp_get_max_threads(void)
   function Get_Max_Threads return Thread_Count;

   -- https://computing.llnl.gov/tutorials/openMP/#OMP_GET_NUM_THREADS
   -- int omp_get_num_threads(void)
   function Get_Num_Threads return Thread_Count;

end OMP;
