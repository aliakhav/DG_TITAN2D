#include "../header/hpfem.h"

// input x and output ax are in global numbering
void mat_vec(double* x, double* ax, double* store,double* val, int* rpntr, int* bpntr, 
	      int* bindx, int* indx, int num_blocks, int subdof, 
	     int totdof,int* ix)
{
  int i,j,k,l;
  int one = 1;
  for(i=0;i<totdof;i++)  *(store+i) = 0;

  //matrix vector product
  int numrows, numcols, numblks, val_loc;
  for(i=0;i<num_blocks;i++)  //goes through all block rows
    {
      numrows = *(rpntr+i+1) - *(rpntr+i);  //number of rows in this block
    //  numblks = *(bpntr+i+1) - *(bpntr+i);
      for(j=*(bpntr+i);j<*(bpntr+i+1);j++) //go through non-zero blocks in this row
	{
	  numcols = *(rpntr+*(bindx+j)+1) - *(rpntr+*(bindx+j));
	//  val_loc = *(indx+j);

	//  for(k=*(rpntr+*(bindx+j));k<*(rpntr+*(bindx+j)+1);k++)
	  //  for(l=*(rpntr+i);l<*(rpntr+i+1);l++)
	  
	  for(k=0;k<numcols;k++)
	    for(l=0;l<numrows;l++)
	      *(store+*(ix+l+*(rpntr+i))) += *(val+*(indx+j)+l+k*numrows)* 
		*(x+*(ix+k+*(rpntr+*(bindx+j)))); 
	}
    }
  
  i=MPI_Allreduce(store, ax,totdof, MPI_DOUBLE,
		  MPI_SUM,MPI_COMM_WORLD);

}

// input x and output ax are in subdomain numbering
void mat_vec(double* x, double* ax, double* val, int* rpntr, int* bpntr, 
	      int* bindx, int* indx, int num_blocks, int subdof)
{
  int i,j,k,l;
  int one = 1;
  for(i=0;i<subdof;i++)  *(ax+i) = 0;

  //matrix vector product
  int numrows, numcols, numblks, val_loc;
  for(i=0;i<num_blocks;i++)  //goes through all block rows
    {
      numrows = *(rpntr+i+1) - *(rpntr+i);  //number of rows in this block
    //  numblks = *(bpntr+i+1) - *(bpntr+i);
      for(j=*(bpntr+i);j<*(bpntr+i+1);j++) //go through non-zero blocks in this row
	{
	  numcols = *(rpntr+*(bindx+j)+1) - *(rpntr+*(bindx+j));
	//  val_loc = *(indx+j);

	//  for(k=*(rpntr+*(bindx+j));k<*(rpntr+*(bindx+j)+1);k++)
	  //  for(l=*(rpntr+i);l<*(rpntr+i+1);l++)
	  
	  for(k=0;k<numcols;k++)
	    for(l=0;l<numrows;l++)
	      *(ax+l+*(rpntr+i)) += *(val+*(indx+j)+l+k*numrows)* 
		*(x+k+*(rpntr+*(bindx+j))); 
	}
    }
}
