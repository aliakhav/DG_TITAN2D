#include "../header/hpfem.h"
#include "../header/refined_neighbor_info.h"

extern void update_neighbor_interprocessor(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr,
					   refined_neighbor* refined_start, int myid, int numprocs);

void update_neighbor_info(HashTable* HT_Elem_Ptr, Element* refined[], 
			  int count, int myid, int numprocs, 
			  HashTable* HT_Node_Ptr, int h_count )/*-h_count for  debugging-*/
{
  int son_order[2];
  Element* EmTemp;
  Element* Neighbor;
  Element* SonTemp;
  Node*    NdTemp;
  unsigned* orig_neighbors;
  int* orig_neigh_proc;
  int which_side;
  unsigned NewNeighbor[2][KEYLENGTH];
  
  int SIDE_SONS[4][2]={{0, 1}, {1, 2}, {2, 3}, {3, 0}};
 
  refined_neighbor* refined_start=new refined_neighbor();
  refined_neighbor* refined_current=refined_start;
  refined_neighbor* refined_new;

  /*cout<<"H_ADAPT  IN  SUBDOMAIN "<<myid<<" ELEMENTS REFINED: "<<count<<endl<<flush;*/
 // printf("H_ADAPT  IN  SUBDOMAIN %d  ELEMENTS REFINED: %d\n",myid, count);  

  for(int i=0; i<count; i++)
    {
      
      EmTemp=refined[i];//-- element ready for refinement

      orig_neighbors=EmTemp->get_neighbors();
      orig_neigh_proc=EmTemp->get_neigh_proc();
      unsigned* Sons=EmTemp->getson();
      unsigned* NeighSons;
      int MyGeneration=EmTemp->get_gen();
      int NeighGeneration;
      int NeighRefined;
      unsigned* Mykey=EmTemp->pass_key();
      int reg;
      int a, b, k;

      for(int j = 0; j<8; j++)
	{
	  if((*(orig_neigh_proc+j)>=0)&&(*(orig_neigh_proc+j)==myid)) //neighbor is at the same proc.
	    {
	      Neighbor = (Element*) HT_Elem_Ptr->lookup(orig_neighbors+j*KEYLENGTH);	
	      assert(Neighbor);
	      NeighGeneration=Neighbor->get_gen();
	      NeighRefined=Neighbor->get_refined_flag();

	      /*-------------------------------------------------
		`  case 1: neighbor has the same gen. and was not refined in this step
	          case 2: neighbor has the same gen. and was refined in this step
		  case 3: neighbor is older and was refined
		  case 4: neighbor is younger and was refined
		  case 5: neighbor is younger and was not refined
		  other : wrong request
	       *-----------------------------------------------*/
	      if(!NeighRefined)
		{
		  assert(NeighGeneration>=MyGeneration);
		  which_side=Neighbor->which_neighbor(Mykey);
		  assert(which_side<4);

		  if(NeighGeneration == MyGeneration)//--case 1
		    {	
		      a=j+1;
		      if (a==4) a=0; 
		      b = j;
		      reg = 1;
		    }
		  else//Neighbor is 1 smaller
		    {
		      a = EmTemp->which_neighbor(orig_neighbors+j*KEYLENGTH);//--case 5
		      if(a == 7) a = 0;
		      else if(a>=4) a = a-3;
		      b = a;
		      reg = 6;
		    }

		  for(k=0; k<KEYLENGTH; k++)
		    {		      		     
		      NewNeighbor[0][k]=*(Sons+a*KEYLENGTH+k);
		      NewNeighbor[1][k]=*(Sons+b*KEYLENGTH+k);
		    }
		  Neighbor->change_neighbor(&NewNeighbor[0][0], which_side, myid, reg);  
		}

	      else//-- neighbor is refined
		{	      
		  if(NeighGeneration == (MyGeneration +1))//--case 4
		    {
		      
		      a=j;/*-- a: my son's number--*/
		      if(j>=4) a=j-3;
		      if(j==7) a=0;
		      
		      which_side=Neighbor->which_neighbor(Mykey); 
		      NdTemp = (Node*)(HT_Node_Ptr->lookup(Neighbor->getNode()+(which_side+4)*KEYLENGTH));
		      NdTemp->putinfo(S_C_CON);
		      //NdTemp->putorder
		      
		      int b1 = which_side;
		      int b2=which_side+1;/*-- b1, b2: neighbor son's number--*/
		      if(b2==4) b2=0;

		      for(int k=0; k<KEYLENGTH; k++)
			NewNeighbor[0][k]=NewNeighbor[1][k]=*(Sons+a*KEYLENGTH+k);
		      reg = 6;//changed from 4
		      
		      unsigned* NeiSon = Neighbor->getson();

		      Neighbor = (Element*)(HT_Elem_Ptr->lookup(NeiSon+b1*KEYLENGTH));
		      Neighbor->change_neighbor(&NewNeighbor[0][0], which_side, myid, reg);

		      Neighbor = (Element*)(HT_Elem_Ptr->lookup(NeiSon+b2*KEYLENGTH));
		      Neighbor->change_neighbor(&NewNeighbor[0][0], which_side, myid, reg);
		      
		    }

		  else if(NeighGeneration == MyGeneration )//-- case 2
		    {
		      which_side=Neighbor->which_neighbor(Mykey);
		      NeighSons=Neighbor->getson();

		      a=j+1;
		      if(a==4) a=0;
		      b=which_side+1;
		      if(b==4) b=0;

		      for(int k=0; k<KEYLENGTH; k++)
			NewNeighbor[0][k]=NewNeighbor[1][k]=*(NeighSons+which_side*KEYLENGTH+k);
		      
		      Element* EmSon= (Element*) HT_Elem_Ptr->lookup(Sons+a*KEYLENGTH);
		      
		      EmSon->change_neighbor(&NewNeighbor[0][0], j, myid, 6);//changed to 6 from 2 lend
		  
		      for(k=0; k<KEYLENGTH; k++)
			NewNeighbor[0][k]=NewNeighbor[1][k]=*(NeighSons+b*KEYLENGTH+k);

		      EmSon= (Element*) HT_Elem_Ptr->lookup(Sons+j*KEYLENGTH);
		      
		      EmSon->change_neighbor(&NewNeighbor[0][0], j, myid, 6);//changed to 6 from 2 lend
		    }

		  else if(NeighGeneration<MyGeneration)//-- case 3  
		    { 
		      if(j>=4) {
			printf("neigh_proc[%d] = %d\n",j,orig_neigh_proc[j]);
			j = j-4;
		      }
			
		      which_side=Neighbor->which_neighbor(Mykey);
		      NeighSons=Neighbor->getson();

		      a=which_side;
		      if(which_side>=4) a=which_side-3;
		      if(which_side==7) a=0;

		      b = j+1;
		      if(b==4) b = 0;
		      Neighbor = (Element*)(HT_Elem_Ptr->lookup(NeighSons+a*KEYLENGTH));

		      int edge_of_son;
		      if(a == which_side) edge_of_son = which_side;
		      else edge_of_son = which_side - 4;
		      
		      for(int k=0; k<KEYLENGTH; k++)
			{		      		     
			  NewNeighbor[0][k]=*(Sons+b*KEYLENGTH+k);
			  NewNeighbor[1][k]=*(Sons+j*KEYLENGTH+k);
			}
		      Neighbor->change_neighbor(&NewNeighbor[0][0], edge_of_son, myid, 3);  
		      
		    }		  		  
		}
	      
	    }//end of if same proc
	  
	  else if(*(orig_neigh_proc+j)!=myid && *(orig_neigh_proc+j) >= 0)
	    {
	      refined_new=new refined_neighbor();
	      refined_current->next=refined_new;
	      
	      NeighGeneration=*(EmTemp->get_neigh_gen()+j);

	      if(NeighGeneration==MyGeneration)
		{
		  assert(j<4);
		  refined_new->set_parameters(*(orig_neigh_proc+j), (orig_neighbors+j*KEYLENGTH), 
					      Mykey, Sons, SIDE_SONS[j], MyGeneration, 1);
		}
	      
	      else
		refined_new->set_parameters(*(orig_neigh_proc+j), (orig_neighbors+j*KEYLENGTH),
					    Mykey, Sons, &SIDE_SONS[j%4][j/4], MyGeneration, 2);//the neighbor is younger...

	      refined_current=refined_new;	    

	      /*for the neighbor information...assume that the neighbor was not 
		refined...if yes then it will be changed in update_inter*/
	      if(j<4)
		{
		  if(NeighGeneration==MyGeneration)
		    {
		      NdTemp=(Node*)HT_Node_Ptr->lookup(EmTemp->getNode()+(j+4)*KEYLENGTH);
		      NdTemp->putinfo(S_C_CON);
		      for(int k=0; k<2; k++)
			{
			  SonTemp=(Element*)HT_Elem_Ptr->lookup(EmTemp->getson()+KEYLENGTH*SIDE_SONS[j][k]);
			  NdTemp=(Node*)HT_Node_Ptr->lookup(SonTemp->getNode()+(j+4)*KEYLENGTH);
			  assert(NdTemp);
			  NdTemp->putinfo(S_S_CON);
			}
		    }

		  else//NeighGen < MyGen
		    {
		      NdTemp=(Node*)HT_Node_Ptr->lookup(EmTemp->getNode()+(j+4)*KEYLENGTH);
		      NdTemp->putinfo(CORNER);
		      for(int k=0; k<2; k++)
			{
			  SonTemp=(Element*)HT_Elem_Ptr->lookup(EmTemp->getson()+KEYLENGTH*SIDE_SONS[j][k]);
			  NdTemp=(Node*)HT_Node_Ptr->lookup(SonTemp->getNode()+(j+4)*KEYLENGTH);
			  assert(NdTemp);
			  NdTemp->putinfo(SIDE);
			}
		      
		    }
		}

	    }//end of is different proc
	  
	} 
    } 

  MPI_Barrier(MPI_COMM_WORLD);
  update_neighbor_interprocessor(HT_Elem_Ptr, HT_Node_Ptr, refined_start, myid, numprocs);
  delete refined_start;
  //cout<<"ready with update: "<<myid<<"\n\n"<<flush;
}
