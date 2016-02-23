
#include "../header/hpfem.h"		      


void vv_ee(HashTable* ht_node_ptr, unsigned* keyP, int start, 
	   int end, int mid, int* nsveb, int myid, int assoc,  
	   Element* EmTemp, VVLink* VVHead, EELink* EEHead,  
	   int* DofCounter)
{
  void*    p;
  Node*    NdTemp;
  VVLink*  VV_new;
  EELink*  EE_new;
  int      order;
  int      which_son;

  VVLink*  VVTail = VVHead;
  EELink*  EETail = EEHead;
  while(VVTail->next)
      VVTail = VVTail->next;
  while(EETail->next)
      EETail = EETail->next;

  VVLink*  VV_old = VVTail;
  EELink*  EE_old = EETail;

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
      *(nsveb+3) += EQUATIONS;
      VV_new = new VVLink(keyP+start*KEYLENGTH);
      VV_old->next = VV_new;    
      VV_new->pre  = VV_old;    
      VV_old = VV_new;
      NdTemp->putdof(*DofCounter, *DofCounter+EQUATIONS-1);
      *DofCounter = *DofCounter+EQUATIONS;
    }
  else if((*NdTemp->getdof() == INIT)&&(NdTemp->getinfo() == S_C_CON))
    {
      order = NdTemp->get_order();
      //order = *(EmTemp->get_order()+start);
      //NdTemp->put_order(order);
      *(nsveb+4) += (order-1)*EQUATIONS;
      EE_new = new EELink(keyP+start*KEYLENGTH);//-- start+4 changed to start
      EE_old->next = EE_new;  
      EE_new->pre  = EE_old;
      EE_old = EE_new; 
      NdTemp->putdof(-3, -3);//-- for avoiding surplus numbering
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
      *(nsveb+3) += EQUATIONS;
      VV_new = new VVLink(keyP+end*KEYLENGTH);
      VV_old->next = VV_new;    
      VV_new->pre  = VV_old;    
      VV_old = VV_new;
      NdTemp->putdof(*DofCounter, *DofCounter+EQUATIONS-1);
      *DofCounter = *DofCounter+EQUATIONS;
    }
  else if((*NdTemp->getdof() == INIT)&&(NdTemp->getinfo() == S_C_CON))
    {
      order = NdTemp->get_order();
      //order = *(EmTemp->get_order()+start);
      //NdTemp->put_order(order);
      *(nsveb+4) += (order-1)*EQUATIONS;
      EE_new = new EELink(keyP+end*KEYLENGTH);//-- end+4 changed to end
      EE_old->next = EE_new;  
      EE_new->pre  = EE_old;
      EE_old = EE_new; 
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
      //order = *(EmTemp->get_order()+start);
      //NdTemp->put_order(order);      
      *(nsveb+4) += (order-1)*EQUATIONS;
      EE_new = new EELink(keyP+mid*KEYLENGTH);//-- start+4 changed to mid
      EE_old->next = EE_new;  
      EE_new->pre  = EE_old;
      EE_old = EE_new; 
      NdTemp->putdof(-3, -3); 
    }
  
}






