
#include "../header/hpfem.h"		      


void nn_ss(HashTable* ht_node_ptr, unsigned* keyP, int start, 
	   int end, int mid, int* nsveb, int myid, int assoc,  
	   Element* EmTemp, NNLink* NNHead, SSLink* SSHead,  
	   int* DofCounter)
{
  void*    p;
  Node*    NdTemp;
  NNLink*  NN_new;
  SSLink*  SS_new;
  int      order;
  int      which_son;

  NNLink*  NNTail = NNHead;
  SSLink*  SSTail = SSHead;
  while(NNTail->next)
      NNTail = NNTail->next;
  while(SSTail->next)
      SSTail = SSTail->next;

  NNLink*  NN_old = NNTail;
  SSLink*  SS_old = SSTail;

  p = ht_node_ptr->lookup(keyP+start*KEYLENGTH);//--starting node of this edge
  NdTemp = (Node*)p;
  if(NdTemp->getinfo() == INIT)//-- for the original mesh
    {
      NdTemp->putinfo(CORNER);
      NdTemp->put_order(1); 
    }
  if((*NdTemp->getdof() == INIT)&&(NdTemp->getinfo() == CORNER))
    {
      NdTemp->put_order(1);
      if(assoc>myid)//--avoid surplus numbering by two procs
	{
	  *(nsveb+1) += EQUATIONS;
	  NN_new = new NNLink(keyP+start*KEYLENGTH);
	  NN_old->next = NN_new;    
	  NN_new->pre  = NN_old;    
	  NN_old = NN_new;
	}
      NdTemp->putdof(*DofCounter, *DofCounter+EQUATIONS-1);
      *DofCounter = *DofCounter+EQUATIONS;
    }
  if((*NdTemp->getdof() == INIT)&&(NdTemp->getinfo() == S_C_CON))
    {
      order = NdTemp->get_order();
      if(assoc>myid)
	{
	  *(nsveb+2) += (order-1)*EQUATIONS;
	  SS_new = new SSLink(keyP+start*KEYLENGTH);//-- start+4 changed to start
	  SS_old->next = SS_new;  
	  SS_new->pre  = SS_old;
	  SS_old = SS_new; 
	}
      NdTemp->putdof(-3, -3);
    }


  p = ht_node_ptr->lookup(keyP+end*KEYLENGTH);//--ending node of this edge
  NdTemp = (Node*)p;
  if(NdTemp->getinfo() == INIT)//-- for the original mesh 
    {
      NdTemp->putinfo(CORNER);
      NdTemp->put_order(1); 
    }
  if((*NdTemp->getdof() == INIT)&&(NdTemp->getinfo() == CORNER))
    {
      NdTemp->put_order(1); 
      if(assoc>myid)//--avoid surplus numbering by two procs 
	{
	  *(nsveb+1) += EQUATIONS;
	  NN_new = new NNLink(keyP+end*KEYLENGTH);
	  NN_old->next = NN_new;    
	  NN_new->pre  = NN_old;    
	  NN_old = NN_new;
	}
      NdTemp->putdof(*DofCounter, *DofCounter+EQUATIONS-1);
      *DofCounter = *DofCounter+EQUATIONS;
    }
  if((*NdTemp->getdof() == INIT)&&(NdTemp->getinfo() == S_C_CON))
    {
      order = NdTemp->get_order();
      if(assoc>myid)
	{
	  *(nsveb+2) += (order-1)*EQUATIONS;
	  SS_new = new SSLink(keyP+end*KEYLENGTH);//-- end+4 changed to end
	  SS_old->next = SS_new;  
	  SS_new->pre  = SS_old;
	  SS_old = SS_new;
	}
      NdTemp->putdof(-3, -3);
    }

      
  p = ht_node_ptr->lookup(keyP+mid*KEYLENGTH);//--middle node of this edge
  NdTemp = (Node*)p;
  if(NdTemp->getinfo() == INIT)
    {
      NdTemp->putinfo(SIDE);
      NdTemp->put_order(*(EmTemp->get_order()+start));
    }
  if((*NdTemp->getdof() == INIT)&&((NdTemp->getinfo() == S_C_CON)||(NdTemp->getinfo() == SIDE)))
    {
      order = NdTemp->get_order();
      if(assoc>myid)
	{
	  *(nsveb+2) += (order-1)*EQUATIONS;
	  SS_new = new SSLink(keyP+mid*KEYLENGTH);//-- start+4 changed to mid
	  SS_old->next = SS_new;  
	  SS_new->pre  = SS_old;
	  SS_old = SS_new;  
	}
      NdTemp->putdof(-3, -3);
    }
  
}
