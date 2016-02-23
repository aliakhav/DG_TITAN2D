#include"../header/hpfem.h"

extern void fhsfc2d_(double, unsigned, unsigned);
extern void hsfc2d(unsigned* , unsigned* , unsigned* );
extern void create_new_node(int, int, int, HashTable*, Node*[], 
			    unsigned[][2], int, int*, int, MatProps*);


//()---new node numbering


//  3---(14)--6---(15)--2
//  |         |         |
//  |         |         |
// (11) (3)  (12) (2)  (13)
//  |         |         |
//  |         |         |
//  7---(9)---E---(10)--5
//  |         |         |
//  |         |         |
// (6)  (0)  (7)  (1)  (8)
//  |         |         |
//  |         |         |
//  0---(4)---4---(5)---1


//max & min coordinates need to be passed later now only for L-shape!!!
//only 4 one step because of the info FLAG!!!
//if the new node is on INTERFACE flag will be -1

void refine(Element* EmTemp, HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr, MatProps* matprops_ptr)
{ 
  //printf("refining element %u %u \n",*(EmTemp->pass_key()), *(EmTemp->pass_key()+1));
  int which;
  Node *n1, *n2, *n3, *n4;
  
  unsigned* KeyTemp;
  Node* NodeTemp[9];
  unsigned NewNodeKey[16][KEYLENGTH];
  Element* Quad9P;
  int numprocs, myid, i;
  Element* neigh_elm;
  unsigned* neigh_node_key;  
  int RefinedNeigh=0;
  int info;
  int other_proc=0;
  int boundary;
  int order;

  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  
  
  KeyTemp = EmTemp->getNode();
  
  for(i=0; i<8; i++)//-- corners and sides
    {
      NodeTemp[i]=(Node*) HT_Node_Ptr->lookup(KeyTemp+i*KEYLENGTH);
      assert(NodeTemp[i] );
    }
  
  order=EmTemp->get_order();
  
  /*filling up the new order array
    str: side orders remain;
    newsides get the higher order of the already existing sides
    bubbles get the order of the old bubble
  */

  NodeTemp[8]=(Node*) HT_Node_Ptr->lookup(EmTemp->pass_key());//-- bubble

  info=S_S_CON;

  //SIDE 0
  if(*(EmTemp->get_neigh_proc())==-1) 
    boundary=1;
  else 
    boundary=0;

  if(boundary == 1 || *(EmTemp->get_neigh_gen()) <= EmTemp->get_gen()) {
    RefinedNeigh = 0;
    info=S_S_CON;
    if(*(EmTemp->get_neigh_proc())!=myid) 
      {
	other_proc=1;
	info=-1;
      }

    which=4;
    //---Fourth new node---
    n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*4);
    n2 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode());
    
    create_new_node(which, 0, 4, HT_Node_Ptr, NodeTemp, NewNodeKey, info,
		    &RefinedNeigh, boundary, matprops_ptr);
    
    //---Fourth old node---
    if(RefinedNeigh || boundary) 
      NodeTemp[4]->putinfo(CORNER);
    else if(other_proc) 
      NodeTemp[4]->putinfo(-1);
    else 
      NodeTemp[4]->putinfo(S_C_CON);
    
    //---Fifth new node---
    which=5;
    n2 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH);
    
    create_new_node(which, 1, 4, HT_Node_Ptr, NodeTemp, NewNodeKey, info, &RefinedNeigh,
		    boundary, matprops_ptr);
  }
  else {
    // fourth new node
    neigh_elm = (Element*) HT_Elem_Ptr->lookup(EmTemp->get_neighbors());
    i = 0;
    which = -1;
    while(i<4 && which == -1) {
      if(compare_key(neigh_elm->get_neighbors()+i*KEYLENGTH, EmTemp->pass_key()))
	which = i;
      i++;
    }
    assert(which != -1);
    neigh_node_key = neigh_elm->getNode();
    for(i=0;i<KEYLENGTH;i++)
      NewNodeKey[4][i] = neigh_node_key[i+(which+4)*KEYLENGTH];
    n1 = (Node*) HT_Node_Ptr->lookup(NewNodeKey[4]);
    if(neigh_elm->get_refined_flag() == 0 || neigh_elm->get_refined_flag() == GHOST)
      n1->putinfo(SIDE);
    else
      n1->putinfo(S_C_CON);
    //fourth old node
    NodeTemp[4]->putinfo(CORNER);
    if(other_proc) 
      NodeTemp[4]->putinfo(-1);
    // fifth new node
    neigh_elm = (Element*) HT_Elem_Ptr->lookup(EmTemp->get_neighbors()+4*KEYLENGTH);
    i = 0;
    which = -1;
    while(i<4 && which == -1) {
      if(compare_key(neigh_elm->get_neighbors()+i*KEYLENGTH, EmTemp->pass_key()))
	which = i;
      i++;
    }
    assert(which != -1);
    neigh_node_key = neigh_elm->getNode();
    for(i=0;i<KEYLENGTH;i++)
      NewNodeKey[5][i] = neigh_node_key[i+(which+4)*KEYLENGTH];
    n1 = (Node*) HT_Node_Ptr->lookup(NewNodeKey[5]);
    if(neigh_elm->get_refined_flag() == 0 || neigh_elm->get_refined_flag() == GHOST)
      n1->putinfo(SIDE);
    else
      n1->putinfo(S_C_CON);
  }

//+++++++++++++++++++++++++++SIDE1

  if(*(EmTemp->get_neigh_proc()+1)==-1) 
    boundary=1;
  else 
    boundary=0;

  if(boundary == 1 || *(EmTemp->get_neigh_gen()+1) <= EmTemp->get_gen()) {
    RefinedNeigh = 0;
    info=S_S_CON;  
    if(*(EmTemp->get_neigh_proc()+1)!=myid)// && *(EmTemp->get_neigh_proc()+1)>0) 
      {
	other_proc=1;
	info=-1;
      }
    else 
      other_proc=0;

    //---Eight new node---
    which=8;
    n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*5);
    n2 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH);
    
    create_new_node(which, 1, 5, HT_Node_Ptr, NodeTemp, NewNodeKey, info,
		    &RefinedNeigh, boundary, matprops_ptr);
    
    //---Fifth old node---
    if(RefinedNeigh || boundary)
      NodeTemp[5]->putinfo(CORNER);
    else { 
      if(other_proc) 
	NodeTemp[5]->putinfo(info);
      else 
	NodeTemp[5]->putinfo(S_C_CON);
    }
    
    //---Thirteenth new node---
    which=13;
    n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*5);
    n2 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*2);
    
    create_new_node(which, 2, 5, HT_Node_Ptr, NodeTemp, NewNodeKey, info, 
		    &RefinedNeigh, boundary, matprops_ptr);
  }
  else {
    // eighth new node
    neigh_elm = (Element*) HT_Elem_Ptr->lookup(EmTemp->get_neighbors()+KEYLENGTH);
    i = 0;
    which = -1;
    while(i<4 && which == -1) {
      if(compare_key(neigh_elm->get_neighbors()+i*KEYLENGTH, EmTemp->pass_key()))
	which = i;
      i++;
    }
    assert(which != -1);
    neigh_node_key = neigh_elm->getNode();
    for(i=0;i<KEYLENGTH;i++)
      NewNodeKey[8][i] = neigh_node_key[i+(which+4)*KEYLENGTH];
    n1 = (Node*) HT_Node_Ptr->lookup(NewNodeKey[8]);
    if(neigh_elm->get_refined_flag() == 0 || neigh_elm->get_refined_flag() == GHOST)
      n1->putinfo(SIDE);
    else
      n1->putinfo(S_C_CON);
    // fifth old node
    NodeTemp[5]->putinfo(CORNER);
    if(other_proc) 
      NodeTemp[5]->putinfo(-1);
    // thirteenth new node
    neigh_elm = (Element*) HT_Elem_Ptr->lookup(EmTemp->get_neighbors()+5*KEYLENGTH);
    i = 0;
    which = -1;
    while(i<4 && which == -1) {
      if(compare_key(neigh_elm->get_neighbors()+i*KEYLENGTH, EmTemp->pass_key()))
	which = i;
      i++;
    }
    assert(which != -1);
    neigh_node_key = neigh_elm->getNode();
    for(i=0;i<KEYLENGTH;i++)
      NewNodeKey[13][i] = neigh_node_key[i+(which+4)*KEYLENGTH];
    n1 = (Node*) HT_Node_Ptr->lookup(NewNodeKey[13]);
    if(neigh_elm->get_refined_flag() == 0 || neigh_elm->get_refined_flag() == GHOST)
      n1->putinfo(SIDE);
    else
      n1->putinfo(S_C_CON);
  }

  //+++++++++++++++++++++++++++SIDE2

  if(*(EmTemp->get_neigh_proc()+2)==-1) 
    boundary=1;
  else 
    boundary=0;
  
  if(boundary == 1 || *(EmTemp->get_neigh_gen()+2) <= EmTemp->get_gen()) {
    info=S_S_CON;
    
    if(*(EmTemp->get_neigh_proc()+2)!=myid)// && *(EmTemp->get_neigh_proc()+2)>0) 
      {
	other_proc=1;
	info=-1;
      }
    else 
      other_proc=0;

    RefinedNeigh=0;
    
    //---Fourteenth new node---
    which=14;
    n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*3);
    n2 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*6);
    
    create_new_node(which, 3, 6, HT_Node_Ptr, NodeTemp, NewNodeKey, info,
		    &RefinedNeigh, boundary, matprops_ptr);
    
    //---Sixth old node---
    if(RefinedNeigh || boundary) 
      NodeTemp[6]->putinfo(CORNER);
    else if(other_proc) 
      NodeTemp[6]->putinfo(-1);
    else
      NodeTemp[6]->putinfo(S_C_CON);
    
    //---Fifteenth new node---
    which=15;
    // geoflow info
    n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*6);
    n2 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*2);
    
    create_new_node(which, 2, 6, HT_Node_Ptr, NodeTemp, NewNodeKey, info,
		    &RefinedNeigh, boundary, matprops_ptr);
  }
  else {
    // fourteenth new node
    neigh_elm = (Element*) HT_Elem_Ptr->lookup(EmTemp->get_neighbors()+6*KEYLENGTH);
    i = 0;
    which = -1;
    while(i<4 && which == -1) {
      if(compare_key(neigh_elm->get_neighbors()+i*KEYLENGTH, EmTemp->pass_key()))
	which = i;
      i++;
    }
    assert(which != -1);
    neigh_node_key = neigh_elm->getNode();
    for(i=0;i<KEYLENGTH;i++)
      NewNodeKey[14][i] = neigh_node_key[i+(which+4)*KEYLENGTH];
    n1 = (Node*) HT_Node_Ptr->lookup(NewNodeKey[14]);
    if(neigh_elm->get_refined_flag() == 0 || neigh_elm->get_refined_flag() == GHOST)
      n1->putinfo(SIDE);
    else
      n1->putinfo(S_C_CON);
    // sixth old node
    NodeTemp[6]->putinfo(CORNER);
    if(other_proc) 
      NodeTemp[6]->putinfo(-1);
    // fifteenth new node
    neigh_elm = (Element*) HT_Elem_Ptr->lookup(EmTemp->get_neighbors()+2*KEYLENGTH);
    i = 0;
    which = -1;
    while(i<4 && which == -1) {
      if(compare_key(neigh_elm->get_neighbors()+i*KEYLENGTH, EmTemp->pass_key()))
	which = i;
      i++;
    }
    assert(which != -1);
    neigh_node_key = neigh_elm->getNode();
    for(i=0;i<KEYLENGTH;i++)
      NewNodeKey[15][i] = neigh_node_key[i+(which+4)*KEYLENGTH];
    n1 = (Node*) HT_Node_Ptr->lookup(NewNodeKey[15]);
    if(neigh_elm->get_refined_flag() == 0 || neigh_elm->get_refined_flag() == GHOST)
      n1->putinfo(SIDE);
    else
      n1->putinfo(S_C_CON);
  }

  //+++++++++++++++++++++++++++SIDE 3

  if(*(EmTemp->get_neigh_proc()+3)==-1) 
    boundary=1;
  else 
    boundary=0;

  if(boundary == 1 || *(EmTemp->get_neigh_gen()+3) <= EmTemp->get_gen()) {
    info=S_S_CON;
    
    if(*(EmTemp->get_neigh_proc()+3)!=myid)  //&& *(EmTemp->get_neigh_proc()+3)>0) 
      {
	other_proc=1;
	info=-1;
      }
    else 
      other_proc=0;  
    
    RefinedNeigh=0;
    
    //---Sixth new node----
    which=6;
    n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*7);
    n2 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode());
    
    create_new_node(which, 0, 7, HT_Node_Ptr, NodeTemp, NewNodeKey, info,
		    &RefinedNeigh, boundary, matprops_ptr);
    
    //---Seventh old node---
    if(RefinedNeigh || boundary) 
      NodeTemp[7]->putinfo(CORNER);
    else if(other_proc) 
      NodeTemp[7]->putinfo(-1);
    else 
      NodeTemp[7]->putinfo(S_C_CON);
    
    //---Eleventh new node---
    which=11;
    n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*7);
    n2 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*3);
    
    create_new_node(which, 3, 7, HT_Node_Ptr, NodeTemp, NewNodeKey, info,
		    &RefinedNeigh, boundary, matprops_ptr);
  }
  else {
    // sixth new node
    neigh_elm = (Element*) HT_Elem_Ptr->lookup(EmTemp->get_neighbors()+7*KEYLENGTH);
    i = 0;
    which = -1;
    while(i<4 && which == -1) {
      if(compare_key(neigh_elm->get_neighbors()+i*KEYLENGTH, EmTemp->pass_key()))
	which = i;
      i++;
    }
    assert(which != -1);
    neigh_node_key = neigh_elm->getNode();
    for(i=0;i<KEYLENGTH;i++)
      NewNodeKey[6][i] = neigh_node_key[i+(which+4)*KEYLENGTH];
    n1 = (Node*) HT_Node_Ptr->lookup(NewNodeKey[6]);
    if(neigh_elm->get_refined_flag() == 0 || neigh_elm->get_refined_flag() == GHOST)
      n1->putinfo(SIDE);
    else
      n1->putinfo(S_C_CON);
    // seventh old node
    NodeTemp[7]->putinfo(CORNER);
    if(other_proc) 
      NodeTemp[7]->putinfo(-1);
    // eleventh new node
    neigh_elm = (Element*) HT_Elem_Ptr->lookup(EmTemp->get_neighbors()+3*KEYLENGTH);
    i = 0;
    which = -1;
    while(i<4 && which == -1) {
      if(compare_key(neigh_elm->get_neighbors()+i*KEYLENGTH, EmTemp->pass_key()))
	which = i;
      i++;
    }
    assert(which != -1);
    neigh_node_key = neigh_elm->getNode();
    for(i=0;i<KEYLENGTH;i++)
      NewNodeKey[11][i] = neigh_node_key[i+(which+4)*KEYLENGTH];
    n1 = (Node*) HT_Node_Ptr->lookup(NewNodeKey[11]);
    if(neigh_elm->get_refined_flag() == 0 || neigh_elm->get_refined_flag() == GHOST)
      n1->putinfo(SIDE);
    else
      n1->putinfo(S_C_CON);
  }




  //++++++++++++++++INTERNAL SIDE NODES 7, 8OLD, 12, 9, 10

  RefinedNeigh=0;
  boundary=0;
  info=SIDE;

  //---Seventh new node---

  which=7;
  // geoflow info
  n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*4);
  n2 = (Node*) HT_Node_Ptr->lookup(EmTemp->pass_key());
 
  create_new_node(which, 4, 8, HT_Node_Ptr, NodeTemp, NewNodeKey, info, 
		  &RefinedNeigh, boundary, matprops_ptr);

  NodeTemp[8]->putinfo(CORNER);//changing the old bubble


  //---Twelwth new node---

  which=12;
  // geoflow info
  n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*6);

  create_new_node(which, 6, 8, HT_Node_Ptr, NodeTemp, NewNodeKey, info,
		  &RefinedNeigh, boundary, matprops_ptr);


  //---Ninth new node---

  which=9;
  // geoflow info
  n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*7);

  create_new_node(which, 7, 8, HT_Node_Ptr, NodeTemp, NewNodeKey, info,
		  &RefinedNeigh, boundary, matprops_ptr);

  //---Tenth new node---

  which=10;
  // geoflow info
  n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*5);

  create_new_node(which, 5, 8, HT_Node_Ptr, NodeTemp, NewNodeKey, info,
		  &RefinedNeigh, boundary, matprops_ptr);


  //+++++++++++++++++++THE NEW BUBBLES 0, 1, 2, 3

  info=BUBBLE;

  //---0th new node---

  NodeTemp[0]= (Node*) HT_Node_Ptr->lookup(NewNodeKey[6]);
  NodeTemp[1]= (Node*) HT_Node_Ptr->lookup(NewNodeKey[7]);
  which=0;
  n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode());
  n3 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*4);
  n4 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*7);

  create_new_node(which, 0, 1, HT_Node_Ptr, NodeTemp, NewNodeKey, info,
		  &RefinedNeigh, boundary, matprops_ptr);
  


  //---1st new node---

  NodeTemp[0]= (Node*) HT_Node_Ptr->lookup(NewNodeKey[7]);
  NodeTemp[1]= (Node*) HT_Node_Ptr->lookup(NewNodeKey[8]);
  which=1;
  n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH);
  n4 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*5);

  create_new_node(which, 0, 1, HT_Node_Ptr, NodeTemp, NewNodeKey, info,
		  &RefinedNeigh, boundary, matprops_ptr);



  //---2nd new node---

  NodeTemp[0]= (Node*) HT_Node_Ptr->lookup(NewNodeKey[12]);
  NodeTemp[1]= (Node*) HT_Node_Ptr->lookup(NewNodeKey[13]);
  which=2;
  n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*2);
  n3 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*5);
  n4 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*6);

  create_new_node(which, 0, 1, HT_Node_Ptr, NodeTemp, NewNodeKey, info,
		  &RefinedNeigh, boundary, matprops_ptr);


  //---3rd new node---

  NodeTemp[0]= (Node*) HT_Node_Ptr->lookup(NewNodeKey[11]);
  NodeTemp[1]= (Node*) HT_Node_Ptr->lookup(NewNodeKey[12]);
  which=3;
  n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*6);
  n3 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*3);
  n4 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*7);

  create_new_node(which, 0, 1, HT_Node_Ptr, NodeTemp, NewNodeKey, info, 
		  &RefinedNeigh, boundary, matprops_ptr);


  //---NEW ELEMENTS---

  unsigned nodes[9][KEYLENGTH];
  unsigned neigh[8][KEYLENGTH];
  unsigned* orig_neighbors=EmTemp->get_neighbors();
  int* orig_neigh_proc=EmTemp->get_neigh_proc();
  int neigh_proc[8];
  unsigned father[KEYLENGTH];
  int generation=EmTemp->get_gen()+1;
  int* orig_neigh_gen=EmTemp->get_neigh_gen();
  int neigh_gen[4];
  int material=EmTemp->get_material();

  for(i=0; i<KEYLENGTH; i++)
    father[i]=*(EmTemp->pass_key()+i);
  Element* father_elm = EmTemp;
  double coord[DIMENSION];
  //---0th new element---


  //the nodes
  
  for(i=0; i<KEYLENGTH; i++)
    {
      nodes[0][i]=*(KeyTemp+i);
      nodes[1][i]=*(KeyTemp+4*KEYLENGTH+i);
      nodes[2][i]=*((EmTemp->pass_key())+i);
      nodes[3][i]=*(KeyTemp+7*KEYLENGTH+i);
      nodes[4][i]=NewNodeKey[4][i];
      nodes[5][i]=NewNodeKey[7][i];
      nodes[6][i]=NewNodeKey[9][i];
      nodes[7][i]=NewNodeKey[6][i];
      nodes[8][i]=NewNodeKey[0][i];
    }
  n1 = (Node*) HT_Node_Ptr->lookup(&(nodes[8][0]));
  for(i=0;i<DIMENSION;i++)
    coord[i] = *(n1->get_coord()+i);
  //neighbors
 
  for(i=0; i<KEYLENGTH; i++)
    {
      neigh[0][i]=neigh[4][i]=*(orig_neighbors+i);
      neigh[1][i]=neigh[5][i]=NewNodeKey[1][i];
      neigh[2][i]=neigh[6][i]=NewNodeKey[3][i];
      if(*(EmTemp->get_neigh_proc()+7)!=-2)
	neigh[3][i]=neigh[7][i]=*(orig_neighbors+7*KEYLENGTH+i);
      else
	neigh[3][i]=neigh[7][i]=*(orig_neighbors+3*KEYLENGTH+i);
    }

  //process of the neighbors

  neigh_proc[0]=*(orig_neigh_proc);
  neigh_proc[1]=myid;
  neigh_proc[2]=myid;
  if(*(orig_neigh_proc+7)!=-2) neigh_proc[3]=*(orig_neigh_proc+7);//depending if the neighboring element is already refined
  else neigh_proc[3]=*(orig_neigh_proc+3);

  neigh_proc[4]=neigh_proc[5]=neigh_proc[6]=neigh_proc[7]=-2;

  neigh_gen[0]=*orig_neigh_gen;
  neigh_gen[1]=generation;
  neigh_gen[2]=generation;
  neigh_gen[3]=*(orig_neigh_gen+3);
 
  double err = (*(EmTemp->get_el_error()))*.5; //added by jp oct11
  double sol = (*(EmTemp->get_el_solution()))*.5;//added by jp oct11
 
  int elm_loc[2] = {*(EmTemp->get_elm_loc()) *2, *(EmTemp->get_elm_loc()+1) *2};
  Quad9P=new Element(nodes, neigh, neigh_proc, generation, 
		     order, neigh_gen, material, EmTemp,
		     HT_Elem_Ptr, HT_Node_Ptr, myid, matprops_ptr, elm_loc);
  Quad9P->put_which_son(0);//--by jp, 0 means son 0
  Quad9P->put_order(order);
  assert(order < 100000000000);
  Quad9P->update_ndof();

  // delete an existing element if it has the same key...
  Element* old_elm = (Element*) HT_Elem_Ptr->lookup(Quad9P->pass_key());
  if(old_elm != NULL) {
    HT_Elem_Ptr->remove(old_elm->pass_key());
    delete old_elm;
  }
    
  HT_Elem_Ptr->add(nodes[8], Quad9P);
  Quad9P->project_sol(HT_Node_Ptr, HT_Elem_Ptr, 0);


  //---1st new element---


  //the nodes
  
  for(i=0; i<KEYLENGTH; i++)
    {
      nodes[0][i]=*(KeyTemp+4*KEYLENGTH+i);
      nodes[1][i]=*(KeyTemp+1*KEYLENGTH+i);
      nodes[2][i]=*(KeyTemp+5*KEYLENGTH+i);
      nodes[3][i]=*((EmTemp->pass_key())+i);
      nodes[4][i]=NewNodeKey[5][i];
      nodes[5][i]=NewNodeKey[8][i];
      nodes[6][i]=NewNodeKey[10][i];
      nodes[7][i]=NewNodeKey[7][i];
      nodes[8][i]=NewNodeKey[1][i];
    }
  n1 = (Node*) HT_Node_Ptr->lookup(&(nodes[8][0]));
  for(i=0;i<DIMENSION;i++)
    coord[i] = *(n1->get_coord()+i);

  //neighbors
 
  for(i=0; i<KEYLENGTH; i++)
    {
      if(*(EmTemp->get_neigh_proc()+4)!=-2)
	neigh[0][i]=neigh[4][i]=*(orig_neighbors+4*KEYLENGTH+i);
      else
	neigh[0][i]=neigh[4][i]=*(orig_neighbors+0*KEYLENGTH+i);
      neigh[1][i]=neigh[5][i]=*(orig_neighbors+1*KEYLENGTH+i);
      neigh[2][i]=neigh[6][i]=NewNodeKey[2][i];
      neigh[3][i]=neigh[7][i]=NewNodeKey[0][i];
     
    }


  //process of the neighbors

  neigh_proc[0]=(*(orig_neigh_proc+4)!=-2) ? *(orig_neigh_proc+4) : *(orig_neigh_proc);
  neigh_proc[1]=*(orig_neigh_proc+1);
  neigh_proc[2]=myid;
  neigh_proc[3]=myid;

  neigh_proc[4]=neigh_proc[5]=neigh_proc[6]=neigh_proc[7]=-2;

  neigh_gen[0]=*orig_neigh_gen;
  neigh_gen[1]=*(orig_neigh_gen+1);
  neigh_gen[2]=generation;
  neigh_gen[3]=generation;

  elm_loc[0] += 1;
  Quad9P=new Element(nodes, neigh, neigh_proc, generation, 
		     order, neigh_gen, material, EmTemp, 
		     HT_Elem_Ptr, HT_Node_Ptr, myid, matprops_ptr, elm_loc);
  Quad9P->put_which_son(1);//--by jp
  Quad9P->put_order(order);
  assert(order< 100000000);
  Quad9P->update_ndof();

  //Quad9P->fth_orig_proc = EmTemp->my_orig_proc;
  // printf("myid=%d  son_1\n",myid);
   Quad9P->put_myprocess(myid);
  old_elm = (Element*) HT_Elem_Ptr->lookup(Quad9P->pass_key());
  if(old_elm != NULL) {
    HT_Elem_Ptr->remove(old_elm->pass_key());
    delete old_elm;
  }

  HT_Elem_Ptr->add(nodes[8], Quad9P);
  Quad9P->project_sol(HT_Node_Ptr, HT_Elem_Ptr, 0);


  //---2nd new element---


  //the nodes
  
  for(i=0; i<KEYLENGTH; i++)
    {
      nodes[0][i]=*((EmTemp->pass_key())+i);
      nodes[1][i]=*(KeyTemp+5*KEYLENGTH+i);
      nodes[2][i]=*(KeyTemp+2*KEYLENGTH+i);
      nodes[3][i]=*(KeyTemp+6*KEYLENGTH+i);
      nodes[4][i]=NewNodeKey[10][i];
      nodes[5][i]=NewNodeKey[13][i];
      nodes[6][i]=NewNodeKey[15][i];
      nodes[7][i]=NewNodeKey[12][i];
      nodes[8][i]=NewNodeKey[2][i];
    }
  n1 = (Node*) HT_Node_Ptr->lookup(&(nodes[8][0]));
  for(i=0;i<DIMENSION;i++)
    coord[i] = *(n1->get_coord()+i);

  //neighbors
  

 
  for(i=0; i<KEYLENGTH; i++)
    {
      neigh[0][i]=neigh[4][i]=NewNodeKey[1][i];
      if(*(EmTemp->get_neigh_proc()+5)!=-2)
	neigh[1][i]=neigh[5][i]=*(orig_neighbors+5*KEYLENGTH+i);
      else
	neigh[1][i]=neigh[5][i]=*(orig_neighbors+1*KEYLENGTH+i);
      neigh[2][i]=neigh[6][i]=*(orig_neighbors+2*KEYLENGTH+i);
      neigh[3][i]=neigh[6][i]=NewNodeKey[3][i];

    }


  //process of the neighbors

  neigh_proc[0]=myid;
  neigh_proc[1]=(*(orig_neigh_proc+5)!=-2) ? *(orig_neigh_proc+5) : *(orig_neigh_proc+1);
  neigh_proc[2]=*(orig_neigh_proc+2);
  neigh_proc[3]=myid;

  neigh_proc[4]=neigh_proc[5]=neigh_proc[6]=neigh_proc[7]=-2;

  neigh_gen[0]=generation;
  neigh_gen[1]=*(orig_neigh_gen+1);
  neigh_gen[2]=*(orig_neigh_gen+2);
  neigh_gen[3]=generation;

  elm_loc[1]+= 1;
  Quad9P=new Element(nodes, neigh, neigh_proc, generation, 
		     order, neigh_gen, material, EmTemp, 
		     HT_Elem_Ptr, HT_Node_Ptr, myid, matprops_ptr, elm_loc);
  Quad9P->put_which_son(2);//--by jp
  Quad9P->put_order(order);
  assert(order< 100000000);
  Quad9P->update_ndof();

  //Quad9P->fth_orig_proc = EmTemp->my_orig_proc;
  //printf("myid=%d  son_2\n",myid);
  Quad9P->put_myprocess(myid);
  old_elm = (Element*) HT_Elem_Ptr->lookup(Quad9P->pass_key());
  if(old_elm != NULL) {
    HT_Elem_Ptr->remove(old_elm->pass_key());
    delete old_elm;
  }

  HT_Elem_Ptr->add(nodes[8], Quad9P);
  Quad9P->project_sol(HT_Node_Ptr, HT_Elem_Ptr, 0);

  //---3rd new element---

  //the nodes
  
  for(i=0; i<KEYLENGTH; i++)
    {
      nodes[0][i]=*(KeyTemp+7*KEYLENGTH+i);
      nodes[1][i]=*((EmTemp->pass_key())+i);
      nodes[2][i]=*(KeyTemp+6*KEYLENGTH+i);
      nodes[3][i]=*(KeyTemp+3*KEYLENGTH+i);
      nodes[4][i]=NewNodeKey[9][i];
      nodes[5][i]=NewNodeKey[12][i];
      nodes[6][i]=NewNodeKey[14][i];
      nodes[7][i]=NewNodeKey[11][i];
      nodes[8][i]=NewNodeKey[3][i];
    }
  n1 = (Node*) HT_Node_Ptr->lookup(&(nodes[8][0]));
  for(i=0;i<DIMENSION;i++)
    coord[i] = *(n1->get_coord()+i);

  //neighbors
  for(i=0; i<KEYLENGTH; i++)
    {
      neigh[0][i]=neigh[4][i]=NewNodeKey[0][i];
      neigh[1][i]=neigh[5][i]=NewNodeKey[2][i];
      if(*(EmTemp->get_neigh_proc()+6)!=-2)
	neigh[2][i]=neigh[6][i]=*(orig_neighbors+6*KEYLENGTH+i);
      else
	neigh[2][i]=neigh[6][i]=*(orig_neighbors+2*KEYLENGTH+i);
      
      neigh[3][i]=neigh[6][i]=*(orig_neighbors+3*KEYLENGTH+i);    
      
    }


  //process of the neighbors

  neigh_proc[0]=myid;
  neigh_proc[1]=myid;
  neigh_proc[2]=(*(orig_neigh_proc+6)!=-2) ? *(orig_neigh_proc+6) : *(orig_neigh_proc+2);
  neigh_proc[3]=*(orig_neigh_proc+3);
 
  neigh_proc[4]=neigh_proc[5]=neigh_proc[6]=neigh_proc[7]=-2;

  neigh_gen[0]=generation;
  neigh_gen[1]=generation;
  neigh_gen[2]=*(orig_neigh_gen+2);
  neigh_gen[3]=*(orig_neigh_gen+3);
 
  elm_loc[0] += -1;
  Quad9P=new Element(nodes, neigh, neigh_proc, generation,
		     order, neigh_gen, material, EmTemp, 
		     HT_Elem_Ptr, HT_Node_Ptr, myid, matprops_ptr, elm_loc);
  Quad9P->put_which_son(3);//--by jp
  Quad9P->put_order(order);
  assert(order< 100000000);
  Quad9P->update_ndof();

  //Quad9P->fth_orig_proc = EmTemp->my_orig_proc;
  //printf("myid=%d  son_3\n",myid);
  Quad9P->put_myprocess(myid);
  old_elm = (Element*) HT_Elem_Ptr->lookup(Quad9P->pass_key());
  if(old_elm != NULL) {
    HT_Elem_Ptr->remove(old_elm->pass_key());
    delete old_elm;
  }

  HT_Elem_Ptr->add(nodes[8], Quad9P);
  Quad9P->project_sol(HT_Node_Ptr, HT_Elem_Ptr, 0);



  //---CHANGING THE FATHER---

  EmTemp->putson(&NewNodeKey[0][0]);
  // putting in brother info
  for(i=0;i<4;i++) {
    EmTemp = (Element*) HT_Elem_Ptr->lookup(&NewNodeKey[i][0]);
    EmTemp->put_myprocess(myid);
    EmTemp->putbrothers(&NewNodeKey[0][0]);  //was  EmTemp->putbrothers(&NewNodeKey[i][0]);
  }
  
  return;
}


void create_new_node(int which, int Node1, int Node2, HashTable* HT_Node_Ptr, 
		     Node* NodeTemp[], unsigned NewNodeKey[][KEYLENGTH], int info,
		     int* RefNe, int boundary, MatProps* matprops_ptr)
{
  double NewNodeCoord[2];
  double norm_coord[2];
  unsigned u_norm_coord[2];
  unsigned nkey=2;
  unsigned key[KEYLENGTH];
  Node* NewNode;
  Node* p;
  static double XRange[2]; 
  static double YRange[2];
  int i;

  for(i=0; i<2; i++)
    {
      XRange[i]=*(HT_Node_Ptr->get_Xrange()+i);
      YRange[i]=*(HT_Node_Ptr->get_Yrange()+i);
    }
 
  for(i=0; i<2; i++)    
    NewNodeCoord[i]=(*(NodeTemp[Node1]->get_coord()+i) + *(NodeTemp[Node2]->get_coord()+i))*.5;      
  
  norm_coord[0]=(NewNodeCoord[0]-XRange[0])/(XRange[1]-XRange[0]);
  norm_coord[1]=(NewNodeCoord[1]-YRange[0])/(YRange[1]-YRange[0]);

  fhsfc2d_(norm_coord, &nkey, key);

  for(i=0; i<KEYLENGTH; i++)

    NewNodeKey[which][i]=key[i];


  p=(Node*) HT_Node_Ptr->lookup(key);
  
  if(!p)
    
    {
      NewNode = new Node(key, NewNodeCoord, info, matprops_ptr);

      HT_Node_Ptr->add(key, NewNode);

      p=NewNode;

    }

  else if(*(p->get_coord())!=NewNodeCoord[0] || *(p->get_coord()+1)!=NewNodeCoord[1])
   {
     short same_key=0;
     assert(same_key);
   }

  else {
    //printf("trying to create node %u %u again\n", key[0], key[1]);
    *RefNe=1;      
  }
  
  if(*RefNe || boundary)
    p->putinfo(SIDE);
      
  return;
}






