#ifndef REFINED_NEIGHBOR_INFO_H
#define REFINED_NEIGHBOR_INFO_H



struct refined_neighbor_pack
{
  int orig_gen;  

  unsigned target_element[KEYLENGTH];
  unsigned old_neighbor[KEYLENGTH];
  unsigned sons[2][KEYLENGTH]; /*only those sons which are at that side*/

};

struct refined_neighbor
{
  int target_proc;
  int orig_gen;

  unsigned target_element[KEYLENGTH];
  unsigned old_neighbor[KEYLENGTH];
  unsigned sons[2][KEYLENGTH]; /*only those sons which are at that side*/
  
  refined_neighbor* next;

  refined_neighbor(){next=NULL;};

  void set_parameters(int proc, unsigned* target, unsigned* old, unsigned* new_sons, int* which, int og, int case_flag)
    {
      target_proc=proc;
      orig_gen=og;
 
      for(int i=0; i<KEYLENGTH; i++)
	{
	  target_element[i]=*(target+i);
	  old_neighbor[i]=*(old+i);

	  if(case_flag==1)
	    {
	      sons[0][i]=*(new_sons+(*(which+1)*KEYLENGTH+i));
	      sons[1][i]=*(new_sons+(*which*KEYLENGTH+i));
	    }
	  else if(case_flag==2)
	    {
	      sons[0][i]=*(new_sons+((*(which))*KEYLENGTH+i));
	      sons[1][i]=0;
	    }
	}
      next=NULL;
    };
  

};

#endif
