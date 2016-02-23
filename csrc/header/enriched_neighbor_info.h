#ifndef ENRICHED_NEIGHBOR_INFO_H
#define ENRICHED_NEIGHBOR_INFO_H



struct enriched_neighbor_pack
{
  int order;  

  unsigned target_element[KEYLENGTH];
  unsigned neighbor[KEYLENGTH];


};

struct enriched_neighbor
{
  int target_proc;
  int order;

  unsigned target_element[KEYLENGTH];
  unsigned neighbor[KEYLENGTH];
  
  enriched_neighbor* next;

  enriched_neighbor(){next=NULL;};

  void set_parameters(int proc, unsigned* target, unsigned* old, int tor)
    {
      target_proc=proc;
      order=tor;
 
      for(int i=0; i<KEYLENGTH; i++)
	{
	  target_element[i]=*(target+i);
	  neighbor[i]=*(old+i);
	}
      
      next=NULL;
    };
  

};

#endif
