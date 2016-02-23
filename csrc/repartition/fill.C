#include "../header/hpfem.h"
//#define DEBUG
//for BSFC repartitioning scheme
void Pack_element(Element* sendel, ElemPack* elem, HashTable* HT_Node_Ptr, int destination_proc)
{
  int j,i;
  Node* node;
  
  elem->myprocess = destination_proc;
#ifdef DEBUG
  printf("destination proces %d",destination_proc);
  assert(destination_proc<100000000000000000);
#endif
  elem->generation = sendel->generation;
#ifdef DEBUG
  printf("generation %d",sendel->generation);
  assert(sendel->generation<100000000000000000);
#endif
  elem->material=sendel->material;
#ifdef DEBUG
  printf("material %e \n",sendel->material);
  assert(sendel->material<100000000000000000);
#endif
  for(i=0; i<8; i++)
    {
      elem->neigh_proc[i]=sendel->neigh_proc[i];
      elem->neigh_gen[i]=sendel->neigh_gen[i];
    }
#ifdef DEBUG
 for(i=0; i<8; i++)
   {
     printf("material %d\n",sendel->neigh_proc[i]);
     assert(sendel->neigh_proc[i]<100000000000000000);
     printf("material %d\n",sendel->neigh_gen[i]);
     assert(sendel->neigh_gen[i]<100000000000000000);
   }
#endif
 
 elem->order = sendel->order;
#ifdef DEBUG
 printf("order %d",sendel->order);
 assert(sendel->order<100000000000000000);
#endif
 
 // for(i=0;i<8;i++)
 //   elem->neigh_gen[i] = sendel->neigh_gen[i];
 //#ifdef DEBUG
 // for(i=0;i<8;i++)
 //  {
 //    printf("neigh_gen %d\n",sendel->neigh_gen[i]);
 //    assert(sendel->neigh_gen[i]<100000000000000000);
 //  }
 //#endif

  elem->ndof = sendel->ndof;
#ifdef DEBUG
  printf("elemn ndof %d",sendel->ndof);
  assert(sendel->ndof<100000000000000000);
#endif

  elem->refined = sendel->refined;
#ifdef DEBUG
  printf("sendel->refined %d",sendel->refined);
  assert(sendel->refined<100000000000000000);
#endif

  elem->which_son = sendel->which_son;
#ifdef DEBUG
  printf("sendel->which_son %d",sendel->which_son);
  assert(sendel->which_son<100000000000000000);
#endif

  elem->new_old = sendel->new_old;
#ifdef DEBUG
  printf("sendel->new_old %d \n",sendel->new_old);
  assert(sendel->new_old<100000000000000000);
#endif

  for(i=0; i<KEYLENGTH; i++)
    elem->key[i] = sendel->key[i];
#ifdef DEBUG
  for(i=0; i<KEYLENGTH; i++)
    {
      printf("sendel->key[i] %d \n",sendel->key[i]);
      assert(sendel->key[i]<100000000000000000);
    }
#endif

for(i=0; i<8; i++)
    for(j=0; j<KEYLENGTH; j++)
      {
	elem->node_key[i][j] = sendel->node_key[i][j];
	elem->neighbor[i][j] = sendel->neighbor[i][j];
	if(i<4) {
	  elem->son[i][j] = sendel->son[i][j];
	  elem->brothers[i][j] = sendel->brothers[i][j];
	}
      }
#ifdef DEBUG
  for(i=0; i<8; i++)
    for(j=0; j<KEYLENGTH; j++)
    {
      printf("sendel->nodekey %d",sendel->node_key[i][j]);
      assert(sendel->node_key[i][j]<100000000000000000);
      printf("sendel->neighbor %d\n",sendel->neighbor[i][j]);
      assert(sendel->neighbor[i][j]<100000000000000000);

      if(i<4) 
	{
      printf("sendel->son %d",sendel->son[i][j]);
      assert(sendel->son[i][j]<100000000000000000);
      printf("sendel->brothers %d",sendel->brothers[i][j]);
      assert(sendel->brothers[i][j]<100000000000000000);
	}
    }
#endif

  for(i=0; i<EQUATIONS; i++)
       elem->el_error[i] = sendel->el_error[i];
    //  elem->el_error[i] =0;
#ifdef DEBUG
  for(i=0; i<EQUATIONS; i++)
    {
      printf("sendel->el_error[%d] %e\n",i,elem->el_error[i]);
      assert(sendel->el_error[i]<100000000000000000);
    }
#endif

  for(i=0;i<ELM_DOF;i++) {
    elem->el_solution[i] = sendel->el_solution[i];
    elem->prev_el_solution[i] = sendel->prev_el_solution[i];
  }
#ifdef DEBUG
  for(i=0; i<ELM_DOF; i++)
    {
      printf("sendel->el_solution[%d] %e",i,sendel->el_solution[i]);
      assert(sendel->el_solution[i]<100000000000000000);
      printf("sendel->prev_el_solution[%d] %e",i,sendel->prev_el_solution[i]);
      assert(sendel->prev_el_solution[i]<100000000000000000);
    }
#endif

  elem->lb_weight = sendel->lb_weight;
#ifdef DEBUG
     printf("sendel->lb_weight %e\n",sendel->lb_weight);
     assert(sendel->lb_weight<100000000000000000);
#endif

  elem->opposite_brother_flag = sendel->opposite_brother_flag;
#ifdef DEBUG
     printf("sendel->opposite_brother_flag %d",sendel->opposite_brother_flag);
     assert(sendel->opposite_brother_flag<100000000000000000);
#endif

  for(i=0;i<2;i++)
    elem->elm_loc[i] = sendel->elm_loc[i];

#ifdef DEBUG
  for(i=0;i<2;i++)
    {
      printf("sendel->elm_loc[i] %d",sendel->elm_loc[i]);
      assert(sendel->elm_loc[i]<100000000000000000);
    }
#endif

  //and the node info:
  for(i=0; i<8; i++)
    {
      node = (Node*) HT_Node_Ptr->lookup(elem->node_key[i]);
      assert(node);
      elem->n_info[i] = node->info;
#ifdef DEBUG
     printf("node->info %d \n",node->info);
     assert(node->info<100000000000000000);
#endif

      for(j=0; j<2; j++)
	elem->n_coord[i][j] = node->coord[j];
#ifdef DEBUG
      for(j=0; j<2; j++)
	{
	  printf("node->coord[j] %e",node->coord[j]);
	  assert(node->coord[j]<100000000000000000);
	}
#endif



      elem->node_elevation[i] = node->elevation;
#ifdef DEBUG
      printf("node->elevation %e",node->elevation);
      assert(node->elevation<100000000000000000);
#endif
      
    }
  
  node = (Node*) HT_Node_Ptr->lookup(elem->key);
  assert(node);
  elem->n_info[8] = node->info;
#ifdef DEBUG
     printf("node->info %d",node->info);
     assert(node->info<100000000000000000);
#endif

  for(j=0; j<2; j++)
    elem->n_coord[8][j] = node->coord[j];
#ifdef DEBUG
      for(j=0; j<2; j++)
	{
	  printf("node->coord[j] %e",node->coord[j]);
	  assert(node->coord[j]<100000000000000000);
	}
#endif

  elem->node_elevation[8] = node->elevation;
#ifdef DEBUG
     printf(" node->elevation %e \n",node->elevation);
     printf(" ");
     assert(node->elevation<100000000000000000);
#endif
  
}


void Pack_neighbor(int target_proc, ELinkPtr* EL_head, int* counter, NePtr* packed_neighbor_info)
{

  //creates a packing for the changed neighbor info
  *counter = 0;
  int counter2=0;
  
  ELinkPtr EL_temp = *EL_head;
  NeighborPack* pack_try;

  if(*EL_head)
    {
      while(EL_temp)
	{
	  if(EL_temp->target_proc == target_proc) (*counter)++;
	  EL_temp = EL_temp->next;
	  
	}
      
      if(*counter)
	{
	  pack_try =  new NeighborPack[*counter];
	  EL_temp = *EL_head;
	  
	  
	  while(counter2 < *counter && EL_temp)
	    {
	      
	      if(EL_temp->target_proc == target_proc)
		{	  
		  pack_try[counter2].target_proc = target_proc;
		  pack_try[counter2].new_proc = EL_temp-> new_proc;
		  for(int i=0; i<KEYLENGTH; i++)
		    {
		      pack_try[counter2].elkey[i] = EL_temp->elkey[i];
		      pack_try[counter2].targetkey[i] = EL_temp->targetkey[i];
		    }
		  counter2++;
		}
	      
	      EL_temp = EL_temp->next;
	    }
	}
      *packed_neighbor_info = pack_try;
    }
}
