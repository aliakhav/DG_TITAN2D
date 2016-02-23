#include "../header/hpfem.h"

void deletingIF(InterFace* headPtr)
{

  InterFace* currentPtr=headPtr->next;
  InterFace* tempPtr;

  while(currentPtr)
    {
      tempPtr=currentPtr;
      currentPtr=currentPtr->next;
      delete tempPtr;
    }


}

void deletingVV(VVLink* headPtr)
{

  VVLink* currentPtr=headPtr->next;
  VVLink* tempPtr;

  while(currentPtr)
    {
      tempPtr=currentPtr;
      currentPtr=currentPtr->next;
      delete tempPtr;
    }


}

void deletingEE(EELink* headPtr)
{

  EELink* currentPtr=headPtr->next;
  EELink* tempPtr;

  while(currentPtr)
    {
      tempPtr=currentPtr;
      currentPtr=currentPtr->next;
      delete tempPtr;
    }


}

void deletingBB(BBLink* headPtr)
{

  BBLink* currentPtr=headPtr->next;
  BBLink* tempPtr;

  while(currentPtr)
    {
      tempPtr=currentPtr;
      currentPtr=currentPtr->next;
      delete tempPtr;
    }


}
void deletingSS(SSLink* headPtr)
{

  SSLink* currentPtr=headPtr->next;
  SSLink* tempPtr;

  while(currentPtr)
    {
      tempPtr=currentPtr;
      currentPtr=currentPtr->next;
      delete tempPtr;
    }


}
void deletingNN(NNLink* headPtr)
{

  NNLink* currentPtr=headPtr->next;
  NNLink* tempPtr;

  while(currentPtr)
    {
      tempPtr=currentPtr;
      currentPtr=currentPtr->next;
      delete tempPtr;
    }


}
void deletingdof(dof* headPtr)
{

  dof* currentPtr=headPtr->next;
  dof* tempPtr;

  while(currentPtr)
    {
      tempPtr=currentPtr;
      currentPtr=currentPtr->next;
      delete tempPtr;
    }


}

void deletingRecv(Recv* headPtr)
{
  
  Recv* currentPtr=headPtr->next;
  Recv* tempPtr;

  while(currentPtr)
    {
      tempPtr=currentPtr;
      currentPtr=currentPtr->next;
      delete tempPtr;
    }

}
