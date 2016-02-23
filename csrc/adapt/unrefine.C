#include "../header/hpfem.h"

#define MIN_GENERATION 0


void unrefine(HashTable* El_Table, HashTable* NodeTable, double target,
	      int myid, int nump, int time_step, MatProps* matprops_ptr) {
  int i, j, k, unrefined = 0;
  Element* Curr_El;
  //-------------------go through all the elements of the subdomain------------------------
  HashEntryPtr* buck = El_Table->getbucketptr();
  // start unrefinement
  for(i=0; i<El_Table->get_no_of_buckets(); i++)
    if(*(buck+i))
      {
	HashEntryPtr currentPtr = *(buck+i);
	while(currentPtr)
	  {	      
	    Curr_El = (Element*) currentPtr->value;
	    //  need to get currentPtr->next now since currentPtr might get deleted!
	    currentPtr=currentPtr->next;
	    if(!(Curr_El->get_refined_flag()))//if this is a refined element don't involve!!! 
	      {
		// if this if the original element, don't unrefine.  only son 0 checks for unrefinement!
		if(Curr_El->get_gen() > MIN_GENERATION  && Curr_El->get_which_son() == 0) {
		  unsigned* father = Curr_El->getfather();
		  //check to see if currentPtr might get deleted and if it might, find next ptr that won't
		  if(currentPtr != NULL) {
		    int newnext = 0;
		    while(newnext == 0 && currentPtr != NULL) {
		      Element* nextelm = (Element*) currentPtr->value;
		      if(nextelm->get_gen() == 0 || nextelm->get_which_son() == 0)
			newnext = 1;
		      else 
			currentPtr = currentPtr->next;
		    }
		  }
		  unrefined += Curr_El->find_brothers(El_Table, NodeTable,
						      target, myid, matprops_ptr);
		}
	      }
	  }
      }
  //printf("%d elements unrefined on process %d ====================================!!!!!!!!!!!\n", unrefined*4, myid);

  return;
}


/*
 * brother[0] = father, brother[1-4] = sons
 */
//()---new node numbering

//                     side 2
//            3---(14)--6---(15)--2
//            |         |         |
//            |         |         |
//           (11) (3)  (12) (2)  (13)
//            |         |         |
//            |         |         |
//  side 3    7---(9)---E---(10)--5    side 1
//            |         |         |
//            |         |         |
//           (6)  (0)  (7)  (1)  (8)
//            |         |         |
//            |         |         |
//            0---(4)---4---(5)---1
//                     side 0



//only 4 one step because of the info FLAG!!!
//if the new node is on INTERFACE flag will be -1

void unrefine_elements(Element* brothers[], HashTable* El_Table, HashTable* NodeTable) {
  int i, j, max;
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);  
  Node *nd, *nd1, *nd2;
  Element *elm1, *elm2;
  brothers[0]->refined = 0;

  //printf("proc %d is unrefinining to %u %u @@@@@@@@@@@@@@@@@@\n",myid, *(brothers[0]->pass_key()), *(brothers[0]->pass_key()+1));

  //change node info
  for(i=4;i<8;i++) 
    nd = (Node*) NodeTable->lookup(&brothers[0]->node_key[i][0]);

  // bubble nodes
  // bubble order is the maximum of the sons bubble order
  max = 1;
  for(i=1;i<5;i++) {
    nd = (Node*) NodeTable->lookup(brothers[i]->key);
    if(max < brothers[i]->order)
      max = brothers[i]->order;
    NodeTable->remove(nd->key, 0);
    delete nd;
  }

  nd = (Node*) NodeTable->lookup(brothers[0]->key);
  nd->info = BUBBLE;
  // delete interior SIDE nodes
  nd = (Node*) NodeTable->lookup(&brothers[1]->node_key[5][0]);
  NodeTable->remove(&brothers[1]->node_key[5][0], 0);
  delete nd;
  nd = (Node*) NodeTable->lookup(&brothers[1]->node_key[6][0]);
  NodeTable->remove(&brothers[1]->node_key[6][0], 0);
  delete nd;
  nd = (Node*) NodeTable->lookup(&brothers[3]->node_key[4][0]);
  NodeTable->remove(&brothers[3]->node_key[4][0], 0);
  delete nd;
  nd = (Node*) NodeTable->lookup(&brothers[3]->node_key[7][0]);
  NodeTable->remove(&brothers[3]->node_key[7][0], 0);
  delete nd;
  brothers[0]->myprocess = myid;

  // take care of nodes on the edges of the 4 element patch
  //side 0
  if(brothers[1]->neigh_proc[0] < 0) { // there is a boundary on this side
    nd = (Node*) NodeTable->lookup(&brothers[1]->node_key[4][0]);
    NodeTable->remove(&brothers[1]->node_key[4][0], 0);
    delete nd;
    nd = (Node*) NodeTable->lookup(&brothers[2]->node_key[4][0]);
    NodeTable->remove(&brothers[2]->node_key[4][0], 0);
    delete nd;
    nd = (Node*) NodeTable->lookup(&brothers[1]->node_key[1][0]);
    nd->info = SIDE;
    brothers[0]->neigh_proc[0] = -1;
    brothers[0]->neigh_proc[4] = -2; 
    // neigh_gen info might not be necessary
    brothers[0]->neigh_gen[0] = brothers[0]->neigh_gen[4] = brothers[0]->generation;
  }
  else if(brothers[1]->generation == brothers[1]->neigh_gen[0]) {  // currently same generation on this side
    nd1 = (Node*) NodeTable->lookup(&brothers[1]->node_key[4][0]);
    nd1->info = S_S_CON;
    nd2 = (Node*) NodeTable->lookup(&brothers[2]->node_key[4][0]);
    nd2->info = S_S_CON;
    nd = (Node*) NodeTable->lookup(&brothers[1]->node_key[1][0]);
    nd->info = S_C_CON;
    //neighbors
    elm1 = (Element*) El_Table->lookup(&brothers[1]->neighbor[0][0]);
    elm1->change_neigh_info(brothers[0]->key, brothers[1]->key, brothers[0]->generation);
    elm2 = (Element*) El_Table->lookup(&brothers[2]->neighbor[0][0]);
    elm2->change_neigh_info(brothers[0]->key, brothers[2]->key, brothers[0]->generation);
    for(i=0;i<KEYLENGTH;i++) {
      brothers[0]->neighbor[0][i] = elm1->key[i];
      brothers[0]->neighbor[4][i] = elm2->key[i];
    }
    
    brothers[0]->neigh_gen[0] = brothers[0]->neigh_gen[4] = brothers[0]->generation+1;
    brothers[0]->neigh_proc[0] = brothers[0]->neigh_proc[4] = myid;
  }
  else {  // neighbor is currently previous generation
    nd = (Node*) NodeTable->lookup(&brothers[1]->node_key[4][0]);
    NodeTable->remove(&brothers[1]->node_key[4][0], 0);
    delete nd;
    nd = (Node*) NodeTable->lookup(&brothers[2]->node_key[4][0]);
    NodeTable->remove(&brothers[2]->node_key[4][0], 0);
    delete nd;
    nd = (Node*) NodeTable->lookup(&brothers[1]->node_key[1][0]);
    nd->info = SIDE;
    elm1 = (Element*) El_Table->lookup(&brothers[1]->neighbor[0][0]);
    elm1->change_neigh_info(brothers[0]->key, brothers[1]->key, 
			    brothers[0]->generation);
    for(i=0;i<KEYLENGTH;i++) {
      brothers[0]->neighbor[0][i] = elm1->key[i];
      brothers[0]->neighbor[4][i] = elm1->key[i];
    }
    brothers[0]->neigh_gen[0] = brothers[0]->neigh_gen[4] = brothers[0]->generation;  
    brothers[0]->neigh_proc[0] = myid;
    brothers[0]->neigh_proc[4] = -2;
  }
  //side 1
  if(brothers[2]->neigh_proc[1] < 0) {  // boundary on this side
    nd = (Node*) NodeTable->lookup(&brothers[2]->node_key[5][0]);
    NodeTable->remove(&brothers[2]->node_key[5][0], 0);
    delete nd;
    nd = (Node*) NodeTable->lookup(&brothers[3]->node_key[5][0]);
    NodeTable->remove(&brothers[3]->node_key[5][0], 0);
    delete nd;
    nd = (Node*) NodeTable->lookup(&brothers[2]->node_key[2][0]);
    nd->info = SIDE;
    brothers[0]->neigh_proc[1] = -1;
    brothers[0]->neigh_proc[5] = -2;
    // neigh_gen info is probably not necessary
    brothers[0]->neigh_gen[1] = brothers[0]->neigh_gen[5] = brothers[0]->generation;
  }
  else if(brothers[2]->generation == brothers[2]->neigh_gen[1]) { // currently same generation on this side
    nd1 = (Node*) NodeTable->lookup(&brothers[2]->node_key[5][0]);
    nd1->info = S_S_CON;
    nd2 = (Node*) NodeTable->lookup(&brothers[3]->node_key[5][0]);
    nd2->info = S_S_CON;
    nd = (Node*) NodeTable->lookup(&brothers[2]->node_key[2][0]);
    nd->info = S_C_CON;
    //neighbors
    elm1 = (Element*) El_Table->lookup(&brothers[2]->neighbor[1][0]);
    elm1->change_neigh_info(brothers[0]->key, brothers[2]->key, brothers[0]->generation);
    elm2 = (Element*) El_Table->lookup(&brothers[3]->neighbor[1][0]);
    elm2->change_neigh_info(brothers[0]->key, brothers[3]->key, brothers[0]->generation);
    for(i=0;i<KEYLENGTH;i++) {
      brothers[0]->neighbor[1][i] = elm1->key[i];
	brothers[0]->neighbor[5][i] = elm2->key[i];
    }
    brothers[0]->neigh_gen[1] = brothers[0]->neigh_gen[5] = brothers[0]->generation+1;
    brothers[0]->neigh_proc[5] = brothers[0]->neigh_proc[1] = myid;
  }
  else {  // neighbor is currently previous generation
    nd = (Node*) NodeTable->lookup(&brothers[2]->node_key[5][0]);
    NodeTable->remove(&brothers[2]->node_key[5][0], 0);
    delete nd;
    nd = (Node*) NodeTable->lookup(&brothers[3]->node_key[5][0]);
    NodeTable->remove(&brothers[3]->node_key[5][0], 0);
    delete nd;
    nd = (Node*) NodeTable->lookup(&brothers[2]->node_key[2][0]);
    nd->info = SIDE;
    elm1 = (Element*) El_Table->lookup(&brothers[2]->neighbor[1][0]);
    elm1->change_neigh_info(brothers[0]->key, brothers[2]->key, 
			    brothers[0]->generation);
    for(i=0;i<KEYLENGTH;i++) {
      brothers[0]->neighbor[1][i] = elm1->key[i];
      brothers[0]->neighbor[5][i] = elm1->key[i];
    }
    brothers[0]->neigh_gen[1] = brothers[0]->neigh_gen[5] = brothers[0]->generation;      
    brothers[0]->neigh_proc[1] = myid;
    brothers[0]->neigh_proc[5] = -2;
  }
  //side 2
  if(brothers[3]->neigh_proc[2] < 0) {  // this is a boundary
    nd = (Node*) NodeTable->lookup(&brothers[3]->node_key[6][0]);
    NodeTable->remove(&brothers[3]->node_key[6][0], 0);
    delete nd;
    nd = (Node*) NodeTable->lookup(&brothers[4]->node_key[6][0]);
    NodeTable->remove(&brothers[4]->node_key[6][0], 0);
    delete nd;
    nd = (Node*) NodeTable->lookup(&brothers[3]->node_key[3][0]);
    nd->info = SIDE;
    brothers[0]->neigh_proc[2] = -1;
    brothers[0]->neigh_proc[6] = -2;
    // neigh_gen info is probably not necessary
    brothers[0]->neigh_gen[2] = brothers[0]->neigh_gen[6] = brothers[0]->generation;
  } 
  else if(brothers[3]->generation == brothers[3]->neigh_gen[2]) { // currently same generation on this side
    nd1 = (Node*) NodeTable->lookup(&brothers[3]->node_key[6][0]);
    nd1->info = S_S_CON;
    nd2 = (Node*) NodeTable->lookup(&brothers[4]->node_key[6][0]);
    nd2->info = S_S_CON;
    nd = (Node*) NodeTable->lookup(&brothers[3]->node_key[3][0]);
    nd->info = S_C_CON;
    //neighbors
    elm1 = (Element*) El_Table->lookup(&brothers[3]->neighbor[2][0]);
    elm1->change_neigh_info(brothers[0]->key, brothers[3]->key, brothers[0]->generation);
    elm2 = (Element*) El_Table->lookup(&brothers[4]->neighbor[2][0]);
    elm2->change_neigh_info(brothers[0]->key, brothers[4]->key, brothers[0]->generation);
    for(i=0;i<KEYLENGTH;i++) {
      brothers[0]->neighbor[2][i] = elm1->key[i];
      brothers[0]->neighbor[6][i] = elm2->key[i];
    }
    brothers[0]->neigh_gen[2] = brothers[0]->neigh_gen[6] = brothers[0]->generation+1;
    brothers[0]->neigh_proc[2] = brothers[0]->neigh_proc[6] = myid;
  }
  else {  // neighbor is currently previous generation
    nd = (Node*) NodeTable->lookup(&brothers[3]->node_key[6][0]);
    NodeTable->remove(&brothers[3]->node_key[6][0], 0);
    delete nd;
    nd = (Node*) NodeTable->lookup(&brothers[4]->node_key[6][0]);
    NodeTable->remove(&brothers[4]->node_key[6][0], 0);
    delete nd;
    nd = (Node*) NodeTable->lookup(&brothers[3]->node_key[3][0]);
    nd->info = SIDE;
    elm1 = (Element*) El_Table->lookup(&brothers[3]->neighbor[2][0]);
    elm1->change_neigh_info(brothers[0]->key, brothers[3]->key, 
			   brothers[0]->generation);
    for(i=0;i<KEYLENGTH;i++) {
      brothers[0]->neighbor[2][i] = elm1->key[i];
      brothers[0]->neighbor[6][i] = elm1->key[i];
    }
    brothers[0]->neigh_gen[2] = brothers[0]->neigh_gen[6] = brothers[0]->generation; 
    brothers[0]->neigh_proc[2] = myid;
    brothers[0]->neigh_proc[6] = -2;
  }
  //side 3
  if(brothers[4]->neigh_proc[3] < 0) {  // this is a boundary
    nd = (Node*) NodeTable->lookup(&brothers[4]->node_key[7][0]);
    NodeTable->remove(&brothers[4]->node_key[7][0], 0);
    delete nd;
    nd = (Node*) NodeTable->lookup(&brothers[1]->node_key[7][0]);
    NodeTable->remove(&brothers[1]->node_key[7][0], 0);
    delete nd;
    nd = (Node*) NodeTable->lookup(&brothers[4]->node_key[1][0]);
    nd->info = SIDE;
    brothers[0]->neigh_proc[3] = -1;
    brothers[0]->neigh_proc[7] = -2;
    // neigh_gen info is probably not necessary
    brothers[0]->neigh_gen[3] = brothers[0]->neigh_gen[7] = brothers[0]->generation;      
  }
  else if(brothers[1]->generation == brothers[1]->neigh_gen[3]) { //  else if(brothers[4]->generation == brothers[4]->neigh_gen[3]) { // currently same generation on this side
    nd1 = (Node*) NodeTable->lookup(&brothers[4]->node_key[7][0]);
    nd1->info = S_S_CON;
    nd2 = (Node*) NodeTable->lookup(&brothers[1]->node_key[7][0]);
    nd2->info = S_S_CON;
    nd = (Node*) NodeTable->lookup(&brothers[4]->node_key[0][0]);
    nd->info = S_C_CON;
    //neighbors
    elm1 = (Element*) El_Table->lookup(&brothers[4]->neighbor[3][0]);
    elm1->change_neigh_info(brothers[0]->key, brothers[4]->key, brothers[0]->generation);
    elm2 = (Element*) El_Table->lookup(&brothers[1]->neighbor[3][0]);
    elm2->change_neigh_info(brothers[0]->key, brothers[1]->key, brothers[0]->generation);
    for(i=0;i<KEYLENGTH;i++) {
      brothers[0]->neighbor[3][i] = elm1->key[i];
      brothers[0]->neighbor[7][i] = elm2->key[i];
    }
    brothers[0]->neigh_gen[3] = brothers[0]->neigh_gen[7] = brothers[0]->generation+1;
    brothers[0]->neigh_proc[7] = brothers[0]->neigh_proc[3] = myid;
  }
  else {  // neighbor is currently previous generation
    nd = (Node*) NodeTable->lookup(&brothers[4]->node_key[7][0]);
    NodeTable->remove(&brothers[4]->node_key[7][0], 0);
    delete nd;
    nd = (Node*) NodeTable->lookup(&brothers[1]->node_key[7][0]);
    NodeTable->remove(&brothers[1]->node_key[7][0], 0);
    delete nd;
    nd = (Node*) NodeTable->lookup(&brothers[4]->node_key[0][0]);
    nd->info = SIDE;
    elm1 = (Element*) El_Table->lookup(&brothers[4]->neighbor[3][0]);
    elm1->change_neigh_info(brothers[0]->key, brothers[4]->key, 
			    brothers[0]->generation);
    for(i=0;i<KEYLENGTH;i++) {
      brothers[0]->neighbor[3][i] = elm1->key[i];
      brothers[0]->neighbor[7][i] = elm1->key[i];
    }
    brothers[0]->neigh_gen[3] = brothers[0]->neigh_gen[7] = brothers[0]->generation;      
    brothers[0]->neigh_proc[3] = myid;
    brothers[0]->neigh_proc[7] = -2;
  } 
  // project solution back to the parent element
  for(i=0;i<4;i++)
    for(j=0;j<KEYLENGTH;j++)
      brothers[0]->son[i][j] = *(brothers[i+1]->pass_key()+j);
  brothers[0]->project_sol(NodeTable, El_Table, 1);
  // remove elements
  for(i=1;i<5;i++) {
    //printf(" element %u %u deleted!!!\n",brothers[i]->key[0],brothers[i]->key[1]);
    El_Table->remove(brothers[i]->key, 1);
    delete (brothers[i]);
  }
  

  return;
}


void Element::change_neigh_info(unsigned* fth_key, unsigned* ng_key, int ng_gen) {
  int i,j, which_side = -1, same;
  i = 0;
  while(i<8 && which_side == -1) {
    if(neigh_proc[i] >= 0) {
      j = 0;
      same = 1;
      while(j<KEYLENGTH && same == 1) {
	if(neighbor[i][j] != ng_key[j]) 
	  same = 0;
	j++;
      }
      if(same == 1) {
	if(i<4)
	  which_side = i;
	else
	  which_side = i -4;
      }
    }
    i++;
  }
  assert(which_side >= 0);
  for(i=0;i<KEYLENGTH;i++) {
    neighbor[which_side][i] = fth_key[i];
    neighbor[which_side+4][i] = fth_key[i];
  }

  neigh_gen[which_side] = ng_gen;
  neigh_gen[which_side+4] = ng_gen;
  neigh_proc[which_side+4] = -2; 
  
  return;
} 
