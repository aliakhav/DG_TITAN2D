#include "../header/hpfem.h"

extern void  nn_ss(HashTable*,unsigned*, int, int, int, int*, int, int,
		   Element*, NNLink*, SSLink*, int*);
extern void  vv_ee(HashTable*, unsigned*, int, int, int, int*, int, int,
		   Element*, VVLink*, EELink*, int*);
extern void  bb_bb(HashTable*, HashTable*, int*, BBLink*, int*);
extern void  rm_surplus_dof(NNLink*, SSLink*, HashTable*, HashTable*, int, int, int*);
extern void  gldof(HashTable*, HashTable*, int*, int, int, NNLink*, SSLink*, 
		   VVLink*, EELink*, BBLink*, int*, int*, int*, int*, int*);


/*---numbering dof----------

  input::
 
  HT_Elem_Ptr
  HT_Node_Ptr
  InterFaceHead
  myid, numprocs

  output:
  nn, ss, vv, vv, bb
  total dof of subdomain
  local and global ordering

  ------------------------*/

/*======================
         <3>
     3----6----2
     |         |
     |         |
 <4> 7    8    5 <2>
     |         |
     |         |
     0----4----1
         <1> 
======================*/


int* dofordering(int* nn, int* ss, int* vv, int* ee, int* bb,
		 int* Nnn, int* Sss, int* Vvv, int* Eee, int* Bbb, 
		 int* sub_dof, int* gl_dof, 
		 HashTable* ht_node_ptr,
		 HashTable* ht_elem_ptr,
		 InterFace* infaptr,
		 NNLink* NNHead, VVLink* VVHead, SSLink* SSHead,
		 EELink* EEHead, BBLink* BBHead,
		 int myid, int numprocs)

{
  int            i, j, k, mi, ijl;
  InterFacePtr   InTemp;
  Element*       EmTemp;
  Node*          NdTemp;
  unsigned       KeyTemp[KEYLENGTH];
  HashEntryPtr   entryp;
  unsigned*      keyP;
  void*          p;
  int*           assocP;
  int*           dofP;
  int            power;
  int            interf;
  int            NN_total;
  int*           neigh_gen;
  int            order;

  int  on_off;
  int  Nn=0, Ss=0, Vv=0, Ee=0, Bb=0;
  int  Npre;

  int  DofCounter = 0; 

  int  NNtype = 1;
  int  SStype = 2;
  int  VVtype = 3;
  int  EEtype = 4;
  int  BBtype = 5;
  int  nsvebtype = 6;
  int  keytype = 7;
  int  rmtype = 8;

  static int  CallingCount;
  CallingCount = 0;

  int   Adpt_Counting = 1;//-- used for counting the times of h refinement, pass from calling routine

  //int   Neigh[4];//-- used for indicating how many neighbors on the sides along the interface

  NNPtr  NN_new;
  NNPtr  NN_old;
  SSPtr  SS_new;
  SSPtr  SS_old;
  VVPtr  VV_new;
  VVPtr  VV_old;
  EEPtr  EE_new;
  EEPtr  EE_old;
  BBPtr  BB_new;
  BBPtr  BB_old;

  NN_old = NNHead;
  SS_old = SSHead;
  VV_old = VVHead;
  EE_old = EEHead;
  BB_old = BBHead;

  unsigned k_buf[KEYLENGTH+1];//--the last digit is used as a counter
  int rm = -1;//--indicating rm numbers 
  int nsveb[6];
  nsveb[0] = myid;
  for(i=1;i<6;i++) nsveb[i]=0;

  int* ix;//--ix: mapping (local->global) index

  static int CallCounter;
  CallCounter = 0;

  int elements = ht_elem_ptr->get_no_of_buckets(); 
  int nodes    = ht_node_ptr->get_no_of_buckets();

  /*-------------first hashtable scan: ordering NN-----------------------------------*/
  /*--only free corners will be assigned a corner dof--*/

  for(i=0;i<elements;i++)
    {
      entryp = *(ht_elem_ptr->getbucketptr() + i);
      while(entryp)
	{
	  //for(mi=0;mi<4;mi++) Neighbor_Indicator[mi] = 0;
	  int interface_flag = 0;
	  EmTemp = (Element*)(entryp->value);

	  if(!EmTemp->get_refined_flag())//-- if it was refined we only need to consider its sons
	    {

	      assocP = EmTemp->getassoc();
	      keyP = EmTemp->getNode();

	      /*--reconize NN and SS. local NN dof, global NN link and SS link--*/
	      for(mi=0;mi<4;mi++)/*--scan the sides--*/  
		{
		  if((*(assocP+mi) != -2)&&(*(assocP+mi) != -1)&&(*(assocP+mi)!=myid))
		    {

		      switch(mi)//-- 4 edges
			{
			case 0: /*--edge 1, node 0, 1, 4--*/
			  nn_ss(ht_node_ptr, keyP, 0, 1, 4, nsveb, myid, *(assocP+mi),
				EmTemp, NNHead, SSHead, &DofCounter);
			  break;
			      
			case 1: /*--edge 2, node 1, 2, 5--*/
			  nn_ss(ht_node_ptr, keyP, 1, 2, 5, nsveb, myid, *(assocP+mi),
				EmTemp, NNHead, SSHead, &DofCounter);
			  break;
			  
			case 2:/*--edge 3, node 2, 3, 6--*/
			  nn_ss(ht_node_ptr, keyP, 2, 3, 6, nsveb, myid, *(assocP+mi),
				EmTemp, NNHead, SSHead, &DofCounter);
			  break;
			  
			case 3:/*--edge 4, node 3, 0, 7--*/
			  nn_ss(ht_node_ptr, keyP, 3, 0, 7, nsveb, myid, *(assocP+mi),
				EmTemp, NNHead, SSHead, &DofCounter);			  
			  break;
			}
		    }
		}
	    }
	  entryp = entryp->next;
	}
    }

  *nn = DofCounter;

  /*------------------NN  SS  reconization  finishied, local NN ordering finished----*/
/*
    if(myid==2)
        {

           ofstream dout("hash3.dat", ios::out);
           HashEntryPtr  entryp2;
           Element* dent;
           entryp2 = *(ht_elem_ptr->getbucketptr() + 499);
           while(entryp2)
            {
             dent=(Element*)entryp2->value;
             dout<<entryp2->key[0]<<" "<<entryp2->key[1]<<setw(30);
             dout<<*(dent->pass_key())<<" "<<*(dent->pass_key()+1)<<endl<<flush;
             entryp2=entryp2->next;
            }
        }



*/

  /*------------------ordering local SS----------------------------------------------*/
  for(i=0;i<nodes;i++)
    {
      entryp = *(ht_node_ptr->getbucketptr() + i);
      while(entryp)
	{ 
	  NdTemp = (Node*)(entryp->value);
	  if(((NdTemp->getinfo()==SIDE)||(NdTemp->getinfo()==S_C_CON))&&(*NdTemp->getdof()==-3))
	    //if((NdTemp->getinfo()==SIDE)||(NdTemp->getinfo()==S_C_CON))
	    {
	      order = NdTemp->get_order();
	      NdTemp->putdof(DofCounter, DofCounter+(order-1)*EQUATIONS-1);
	      DofCounter = DofCounter+(order-1)*EQUATIONS;
	    }
	  entryp = entryp->next;
	}
    }

  *ss = DofCounter - *nn;

  /*------------------local SS numbering finished-------------------------------------*/
/*
    if(myid==2)
        {

           ofstream dout("hash3.dat", ios::out);
           HashEntryPtr  entryp2;
           Element* dent;
           entryp2 = *(ht_elem_ptr->getbucketptr() + 499);
           while(entryp2)
            {
             dent=(Element*)entryp2->value;
             dout<<entryp2->key[0]<<" "<<entryp2->key[1]<<setw(30);
             dout<<*(dent->pass_key())<<" "<<*(dent->pass_key()+1)<<endl<<flush;
             entryp2=entryp2->next;
            }
        }

*/


  /*------------------VV  EE reconization and local VV ordering-----------------------*/
  //if(CallCounter == 0) {

   //HashEntryPtr* debug;
   int debug_counter2=0;
  int debug_counter=0;
   int got_it=0;

  for(i=0;i<elements;i++)
    {
      entryp = *(ht_elem_ptr->getbucketptr() + i);

      while(entryp)
	{
	  int interface_flag = 0;
	  EmTemp = (Element*)(entryp->value); 
	  if(!EmTemp->get_refined_flag())
	    {

	      debug_counter2++;
	      assocP = EmTemp->getassoc();
	      keyP = EmTemp->getNode();
	      for(mi=0;mi<4;mi++)//--scan the sides--  
		{
		  if((*(assocP+mi) == -1)||(*(assocP+mi) == myid))
		    {
		      switch(mi)//-- 4 edges
			{
			case 0: /*--edge 1, node 0, 1, 4--*/
			  vv_ee(ht_node_ptr, keyP, 0, 1, 4, nsveb, myid, *(assocP+mi),
				EmTemp, VVHead, EEHead, &DofCounter);
			  break;
			      
			case 1: /*--edge 2, node 1, 2, 5--*/ 
			  vv_ee(ht_node_ptr, keyP, 1, 2, 5, nsveb, myid, *(assocP+mi),
				EmTemp, VVHead, EEHead, &DofCounter);
			  break;
			  
			case 2:/*--edge 3, node 2, 3, 6--*/ 
			  vv_ee(ht_node_ptr, keyP, 2, 3, 6, nsveb, myid, *(assocP+mi),
				EmTemp, VVHead, EEHead, &DofCounter);
			  break;
			  
			case 3:/*--edge 4, node 3, 0, 7--*/ 
			  vv_ee(ht_node_ptr, keyP, 3, 0, 7, nsveb, myid, *(assocP+mi),
				EmTemp, VVHead, EEHead, &DofCounter);			  
			  break;
			}
		    }
		}
	    }
	  entryp = entryp->next;
	}
    }

  *vv = DofCounter - *nn - *ss;

  /*------------------VV EE    reconization  finishied, and local VV ordering finished---- */

/*
    if(myid==2)
        {

           ofstream dout("hash4.dat", ios::out);
           HashEntryPtr  entryp2;
           Element* dent;
           entryp2 = *(ht_elem_ptr->getbucketptr() + 499);
           while(entryp2)
            {
             dent=(Element*)entryp2->value;
             dout<<entryp2->key[0]<<" "<<entryp2->key[1]<<setw(30);
             dout<<*(dent->pass_key())<<" "<<*(dent->pass_key()+1)<<endl<<flush;
             entryp2=entryp2->next;
            }
        }
*/

  /*------------------ordering local EE---------------------------------------------------*/
  
  EE_old = EEHead->next;
  while(EE_old)
    {
      NdTemp = (Node*)(ht_node_ptr->lookup(EE_old->key));
      order = NdTemp->get_order();
      NdTemp->putdof(DofCounter, DofCounter+(order-1)*EQUATIONS-1);
      DofCounter = DofCounter+(order-1)*EQUATIONS;
      
      EE_old = EE_old->next;
    }
  
  *ee = DofCounter - *nn - *ss - *vv;
  
  /*------------------BB reconization, and local BB ordering------------------------------*/

  bb_bb(ht_node_ptr, ht_elem_ptr, nsveb, BBHead, &DofCounter);
  *bb = DofCounter - *nn - *ss - *vv - *ee;

  /*------------------BB reconization  finishied, and local VV ordering finished----------*/

  
  /*----------123456789876543212345678987654321------------------------------
    eliminate the supplus numbering in the following case, for example...   
    +---+
    | 0 |
    +---*---+
    | 2 | 1 |
    +---+---+
     point * will be numbered twice. thus we need to eliminate the second 
     numering in proc 1
     statagy: every proc other than 0 send its global num to 0 and let 0 
     check out the supplus ordering.
     ----------------------------------------------------------------------*/
  rm_surplus_dof(NNHead, SSHead, ht_elem_ptr, ht_node_ptr, myid, numprocs, nsveb);

  
  
  /*-------------global numering------------------------------------------
    NN was finished in duplication elimination 
    only need to do SS, VV, EE and BB
    the statagy is:
    1. each proc broacast its global info got in previous steps
       taking use these pieces of info to finish part of the ordering work
    2. the left are NN and SS. procs exchange their info to complete
       ------------------------------------------------------------------*/

  gldof(ht_elem_ptr, ht_node_ptr, nsveb, myid, numprocs, NNHead, SSHead, 
	VVHead, EEHead, BBHead, &Nn, &Ss, &Vv, &Ee, &Bb);
  
  //added by acbauer
  //need to put in dofP array for nodes that are S_S_CON
  //need dofP[0] = -1 for these though so that they don't 
  //get counted or ordered in the ix array
  k = ht_node_ptr->get_no_of_buckets();  
  for(i=0;i<k;i++)
    {
      entryp = *(ht_node_ptr->getbucketptr() + i);
      while(entryp)
	{  
	  NdTemp = (Node*)(entryp->value);
	  if(NdTemp->getinfo()==S_S_CON)
	    {
	      order = NdTemp->get_order();
	      NdTemp->putdof(INIT,(order-1)*EQUATIONS+INIT-1);
	    }	  
	  entryp = entryp->next;
	}
	
    }
  //acbauer done

  /*----creat ix, for each proc (local ordering->global numbering mapping, index 
    is local, value is global) for parallel solver----------------------------
    stratagy: scan node hashtable...
              if the dof is not equal to the initial value, the node must 
	      be shared by this proc. then dof -> ix index, and glnum -> ix value
	      else ignore the node--------------------------------------------*/
  

  k = ht_node_ptr->get_no_of_buckets();
  ix = new int [DofCounter];//---DofCounter: total dof of the subdomain

  for(i=0;i<k;i++)
    {
      entryp = *(ht_node_ptr->getbucketptr() + i);
      while(entryp)
	{  
	  NdTemp = (Node*)(entryp->value);
	  dofP = NdTemp->getdof();
	  
	  if(*dofP != INIT) {
	    for(j=*dofP;j<=*(dofP+1);j++)
	      {
		ix[j] = NdTemp->getglnum() + j-*dofP -1; //in C++ format
	      }
	  }
	  entryp = entryp->next;
	}
    }
  

  *sub_dof = DofCounter-*bb;  //sub_dof:=nn+ss+vv+ee
  *gl_dof  = Nn+Ss+Vv+Ee;  
  *Nnn     = Nn;
  *Sss     = Ss;
  *Vvv     = Vv;
  *Bbb     = Bb;
  *Eee     = Ee;

  CallCounter++;
  return ix;

}







