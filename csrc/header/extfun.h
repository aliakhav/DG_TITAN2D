#ifndef EXTFUN_H
#define EXTFUN_H
#define DEBUG_HEADER
const int QUADNODES = 9;

inline int intpow(int i, int j) {
  assert(j>= 0);
  int k, l=1;
  for(k=0;k<j;k++)
    l = l*i;
  return(l);
}

void delete_unused_elements_nodes(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr, 
				  int myid);

void remove_neg_pile_height(double pile_in[], double pile_ave, double *pilex, double *piley);

int calc_num_shape_functions(int);

void unrefine(HashTable* El_Table, HashTable* NodeTable, double target, 
	      int, int, int, MatProps*);

void step(HashTable* El_Table, HashTable* NodeTable, int myid, int nump,
	  double end_time, double* time_ptr, MatProps*, int , double*,int);

extern void Read_data(HashTable** NodeTable, int myid, 
		      HashTable** ElemTable, 
		      MatProps*, int*, double*, int*, int*, int*, int*);
extern void Get_element_stiffness(HashTable*, HashTable*, int);

extern void Assemble_stiffness(int, double*, double*, HashTable*, HashTable*, 
			       int, int*, int*, int*, int*, int*, int, int, int*);

extern int* dofordering(int*, int*, int*, int*, int*, int*,int*, int*, int*, 
			int*, int*, int*,
                        HashTable* ht_node_p, HashTable* ht_elem_p,
			InterFace* infap, NNLink* NNP, VVLink* VVP,
			SSLink* SSP, EELink* EEP, BBLink* BBP,
			int md, int num );

extern void Get_the_solution(int loc_tot_dof, int glob_tot_dof, int nn, 
			     int ss, int vv, int ee, int bb, 
			     int Nn, int Ss, int Vv, int Ee, int Bb,
			     int* ix, double* stiff, double* load, double* u);

extern void pcg_sol(int subdof, int totdof, int, int, int* ix, int nump, int myid, double* u,
	     double* rhs, double* stiff, int* rpntr, int* bpntr, 
	     int* bindx, int* indx, int num_blocks, HashTable* BT_Node_Ptr,
	     HashTable* BT_Elem_Ptr);  //diag prec.

extern void pcg_sol(int subdof, int totdof, int, int, int* ix, int nump, int myid, double* u,
	     double* rhs, double* stiff, int* rpntr, int* bpntr, 
	     int* bindx, int* indx, int num_blocks, HashTable* BT_Node_Ptr,
	     HashTable* BT_Elem_Ptr, int nn, int ss, int vv, int ee, 
	     int Nn, int Ss, int Vv, int Ee);  //Nn block prec

extern void pcg_sol2(int subdof, int totdof, int, int, int* ix, int nump, int myid, double* u,
	      double* rhs, double* stiff, int* rpntr, int* bpntr, 
	      int* bindx, int* indx, int num_blocks, HashTable* BT_Node_Ptr,
	      HashTable* BT_Elem_Ptr, int nn, int ss, int vv, int ee, 
	      int Nn, int Ss, int Vv, int Ee);  //full matrix prec

extern void pcg_sol(int, int, int, int, int*, int, int, double*, double*, double*, int*, int*, 
	     int*, int*, int, HashTable*, HashTable*, int, int, int);  //reduced memory

void gmres_sol(int myid, int subdof, int totdof, int bb, int Bb, 
	       HashTable* BT_Node_Ptr, HashTable* BT_Elem_Ptr, 
	       int num_blocks, int* rpntr, 
	       int* bpntr, int* bindx, int* indx, double* load, 
	       double* stiff, int* ix, double* u);  //full vector storage

extern void gmres_sol(int subdof, int totdof, int bb, int Bb,
	       int* ix, int numprocs, int myid, double* u,
	       double* load, double* stiffness, int* rpntr, int* bpntr, 
	       int* bindx, int* indx, int num_blocks, HashTable* BT_Node_Ptr,
	       HashTable* BT_Elem_Ptr, int NNSS, int nnss, int vvee);//reduced memory


extern void Post_process1(double* u, int* ix, HashTable* HT_Node_Ptr, 
			  HashTable* HT_Elem_Ptr, NNLink*, VVLink*, SSLink*, EELink*, BBLink*, 
			  int myid, int numprocs, int h_count);

extern void H_adapt(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr, int, double TARGET, 
		    MatProps*, double*);

extern void P_adapt(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr, int, double);

extern void nodesinfo(HashTable*, HashTable*, int);//--for debuging
extern void elemsinfo(HashTable*, HashTable*, int);//--for debuging
extern void showE(double* u, int* ix, HashTable* HT_Node_Ptr,//--for debuging
		  EELink* Head, int myid, int numprocs);
extern void showS(double* u, int* ix, HashTable* HT_Node_Ptr,//--for debuging
		  SSLink* Head, int myid, int numprocs);
extern void showN(double* u, int* ix, HashTable* HT_Node_Ptr,//--for debuging
		  NNLink* Head, int myid, int numprocs);
extern void showV(double* u, int* ix, HashTable* HT_Node_Ptr,//--for debuging
		  VVLink* Head, int myid, int numprocs);

extern void Pack_element(Element* sendel, ElemPack** elemptr, 
			 HashTable* HT_Node_Ptr, int s_f);


extern void MPI_New_Datatype();

extern void repartition(HashTable*, HashTable*);
extern void repartition(HashTable*, HashTable*, int);

extern void smooth(HashTable*, HashTable*);

extern void smooth_II(HashTable*, HashTable*);

extern void Delete_Table(HashTable*, HashTable*);

extern void all_check(HashTable* eltab, HashTable* ndtab, int myid, int m);

extern void tecplotter(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr, int which, double, MatProps*);

extern void vizplotter(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr, int which, MatProps*);

void viz_output(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr,
		int time_step, double time, int myid, int numprocs, MatProps*);

void meshplotter(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr, int which, MatProps*);

extern void check_p_order(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr);

extern int howmanyelements(HashTable*);	

extern double random_error(HashTable*);

extern void search_object ( HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr, unsigned* key , int myid );

extern int compare_key (unsigned*, unsigned*);

extern void list_elements(HashTable*, HashTable*, int, int);

extern void get_max_key(HashTable*, HashTable*);

extern void objects_per_slot(HashTable*, HashTable*, int, int);

extern void Assemble_bubble(int, int, int, int, double*, int*, 
		     HashTable*, HashTable*, int);

extern int make_block(int, int**, int**, int**, int**, int*,int*,
		HashTable*, HashTable*, int, int, int, int);

#endif
#ifdef DEBUG_HEADER  //debug stuff is in SGI/SUNOS fortran interface 
extern void view_elm(HashTable* El_Table, HashTable* NodeTable, int myid);
extern void Mat_write(int, double*, int, int, int);
extern void Mat_write_f(int, double*, int, int, int);
extern "C" void mat_view_(double*, int*, int*);
extern void mat_restore(int, int, int*, int*, int*, int*, double*);
#endif
