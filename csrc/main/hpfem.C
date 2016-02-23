#include "../header/hpfem.h"

//#define DEBUG
 
int main(int argc, char *argv[]) 
{
  int i, j, k;//-- counters
  
  HashTable*   BT_Node_Ptr; 
  HashTable*   BT_Elem_Ptr; 
  
  //-- MPI
  int   myid, numprocs;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

#ifdef DEBUG
  if (myid==0){
    int w;
    printf("type in a number: \n");
    (void) scanf ("%d", &w);
  } 
//  MPI_Barrier(MPI_COMM_WORLD);
#endif 
  
  /* create new MPI datastructures for class objects */
  MPI_New_Datatype();

  /* read original data from serial preprocessing
     code and then initialize element 
     stiffness routines info */
  double epsilon = 1., intfrictang = 1, bedfrictang = 1, gamma = 1; 
  double frict_tiny = 0.1, mu = .0001, rho = 2200, porosity = 1;
  
  MatProps matprops(intfrictang, bedfrictang, porosity, mu, rho, epsilon, 
		    gamma, frict_tiny, (double) 1, (double) 1, (double) 1);

  int max_time_steps = 3000, numoutput = 1, adaptflag;  
  double end_time = 10000.0;  
  /*
   * viz_flag is used to determine which viz output to use
   * viz_flag%2 == 0 means output tecplotxxxx.plt
   * viz_flag%3 == 0 means output mshplotxxxx.plt
   * viz_flag%5 == 0 means output pady's stuff (viz_filenames.out and viz_outputxxx.plt)
   * viz_flag%7 == 0 means output hdf stuff (not implemented yet)
   */
  int viz_flag = 0;
  int order_flag = 0;  //order flag for time stepping scheme -- not used as of 6/19/03
  Read_data(&BT_Node_Ptr, myid, &BT_Elem_Ptr, &matprops, &max_time_steps, 
	    &end_time, &numoutput, &adaptflag, &viz_flag, &order_flag);
  printf("bed friction angle is %e and internal friction angle is %e, epsilon is %e\n",
	 (double) (matprops.bedfrict*180./PI), 
	 (double) (matprops.intfrict*180./PI),
	 (double) (matprops.epsilon));
  printf("METHOD ORDER %d \n",ORDER);

  double dummyt=-100.;
  double maxFluxIntegral=0; //overall maximum of the integral of the 
  //flux on the element boundary   
  int h_c=0;
  //  H_adapt(BT_Elem_Ptr, BT_Node_Ptr, h_c, dummyt, &matprops);
  //H_adapt(BT_Elem_Ptr, BT_Node_Ptr, h_c, dummyt, &matprops);


  if(viz_flag%3==0)
    meshplotter(BT_Elem_Ptr, BT_Node_Ptr, 0,&matprops);
    
  /*
    cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

                  Time Stepping Loop

    cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    */
  
  int time_step = 0;
  double time = 0;
  while(time_step <= max_time_steps && time < end_time)
    {
      //if(myid ==0) 
	//	printf("doing time_step = %d and time is %e\n",time_step,
	//       time*sqrt(matprops.LENGTH_SCALE/matprops.GRAVITY_SCALE));
      /*  
       *  mesh adaption routines 
       */
      double TARGET = 0.005;
      double UNREFINE_TARGET=GEOFLOW_TINY;
      int h_count = 0;
      if (time_step < 0)
	matprops.frict_tiny=0.1;
      else if (time_step >= 0)
	matprops.frict_tiny=0.000001;

      if(adaptflag != 0) {
	if(time_step%4 == 0 ) {
	  H_adapt(BT_Elem_Ptr, BT_Node_Ptr, h_count, TARGET, &matprops, &maxFluxIntegral);
	}
	else if(time_step%4 == 2 ) { 
	  unrefine(BT_Elem_Ptr, BT_Node_Ptr, UNREFINE_TARGET, myid, numprocs, time_step, &matprops);
	}

	//P_adapt(BT_Elem_Ptr, BT_Node_Ptr, TARGET);
	//CN	if(viz_flag%3==0)
	//CN     meshplotter(BT_Elem_Ptr, BT_Node_Ptr, time_step+1,&matprops);
	
	/*
	 *  mesh repartitioning routines
	 */
	if(time_step % 10 == 1 && numprocs > 1) {
	  delete_ghost_elms(BT_Elem_Ptr, myid);
	  repartition(BT_Elem_Ptr, BT_Node_Ptr, time_step);
	  // move_data(numprocs, myid, BT_Elem_Ptr, BT_Node_Ptr);
	}
      }

      step(BT_Elem_Ptr, BT_Node_Ptr, myid, numprocs, end_time, &time,
	   &matprops, time_step,  &maxFluxIntegral,numoutput); 
      //  printf(" maxFluxIntegral %e in hpfem.C \n ",maxFluxIntegral);
      //  exit(0);

      calc_volume(BT_Elem_Ptr, BT_Node_Ptr, myid, numprocs, &matprops);

      /*
       * output results to file 
       */
      if(time_step % numoutput == 0) {
	move_data(numprocs, myid, BT_Elem_Ptr, BT_Node_Ptr);
	if(viz_flag%3==0)
	  meshplotter(BT_Elem_Ptr, BT_Node_Ptr, time_step+1,&matprops);
      }
      
      time_step++;
    }
  move_data(numprocs, myid, BT_Elem_Ptr, BT_Node_Ptr);
  if(viz_flag%3==0)
    meshplotter(BT_Elem_Ptr, BT_Node_Ptr, time_step+1,&matprops);
  printf("%d Finished -- Final simulation time %e\n",myid,time);
  MPI_Finalize();    
  return(0);  
}

