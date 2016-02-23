#include <iostream>
#include <fstream>
#include <assert.h>
#include "../header/hashtab.h"
#include <mpi.h>
#include <limits.h>

#define HASHTABLE_EXTENDER 1000000
#define MaxBits ( sizeof(unsigned) * CHAR_BIT )
#define IScale  ((unsigned)((MaxBits <= 32) ? ~(0u) : (0xffffffff << (MaxBits - 32))))

HashTable:: HashTable(unsigned* min, unsigned* max, int size, int prime)   
{
  MinKey[0] = *min;
  MinKey[1] = *(min+1); 
  MaxKey[0] = *max;
  MaxKey[1] = *(max+1);
  // extend the hashtable bounds a little bit to make it more efficient for adaptive meshes
  unsigned hashtable_extender = HASHTABLE_EXTENDER;
  if(MinKey[0] >= hashtable_extender)
    MinKey[0] -= hashtable_extender;
  else
    MinKey[0] = 0;
  unsigned umax = IScale;
  if((hashtable_extender/2 + MaxKey[0]/2) <= (umax/2))
    MaxKey[0] += hashtable_extender;
  else 
    MaxKey[0] = umax;
  
  NBUCKETS = size;
  PRIME  = prime;
  Range  = *(MaxKey);
  bucket = new HashEntryPtr[NBUCKETS];  

  for(int i = 0; i< NBUCKETS;i++)
    *(bucket+i) = 0;

  /*  MaxMinX[0]=MaxMinY[0]=-1;
  
  MaxMinX[1]=MaxMinY[1]=1;*/
  
}


HashTable::HashTable(unsigned* min, unsigned* max, int size, int prime, double  XR[], double YR[])
{
  int i;
  MinKey[0] = *min;
  MinKey[1] = *(min+1);
  MaxKey[0] = *max;
  MaxKey[1] = *(max+1);
  // extend the hashtable bounds a little bit to make it more efficient for adaptive meshes
  unsigned hashtable_extender = HASHTABLE_EXTENDER;
  if(MinKey[0] >= hashtable_extender)
    MinKey[0] -= hashtable_extender;
  else
    MinKey[0] = 0;
  unsigned umax = IScale;
  if((hashtable_extender/2 + MaxKey[0]/2) <= (umax/2))
    MaxKey[0] += hashtable_extender;
  else 
    MaxKey[0] = umax;

  NBUCKETS = size;
  PRIME  = prime;
  Range  = *(MaxKey);
  bucket = new HashEntryPtr[NBUCKETS];

  for(i = 0; i< NBUCKETS;i++)
    *(bucket+i) = 0;

  for(i=0; i<2; i++)
   {
     Xrange[i]=XR[i];
     Yrange[i]=YR[i];
   }

}


HashTable:: ~HashTable()              //evacuate the table
{
  for (int i = 0; i<NBUCKETS; i++){
    HashEntryPtr p = *(bucket+i); 
    while (p){ 
       HashEntryPtr p_next = p->next; 
       delete p; 
       p = p_next;  
    } 
  }
  delete [] bucket;
}


HashEntryPtr
HashTable:: searchBucket(HashEntryPtr p, unsigned* keyi){
  int i;
  while(p) {
    for(i=0;i<KEYLENGTH;i++) {
      if(p->key[i] != *(keyi+i))
      { 
        p = p->next;                     
        break;                          //not found, check next element
      } 
      else if(i == KEYLENGTH -1)            
        return p;                       //found, return element pointer
    }   
  }
  return p;
}


HashEntryPtr
HashTable:: addElement(int entry, unsigned key[]){  
  
  HashEntryPtr p = new HashEntry(key);  
  
  if(*(bucket+entry)) //this place is already occupied
    {
      HashEntryPtr currentPtr = *(bucket+entry);
      while(currentPtr!=0 && (key[0] > currentPtr->key[0] ))
	{
	  p->pre=currentPtr;
	  currentPtr=currentPtr->next;
	}
      
       if(currentPtr!=0 && key[0] == currentPtr->key[0])
	 {
            while(currentPtr!=0 && (key[1] > currentPtr->key[1] ))
         	{
	         p->pre=currentPtr;
	         currentPtr=currentPtr->next;
        	}

         }
	 
      if(currentPtr) currentPtr->pre=p;
      p->next = currentPtr;
      currentPtr = p->pre;
      //assert(currentPtr);
      if(currentPtr) currentPtr->next=p;
      else bucket[entry]=p;

    }			

  //  p->next = *(bucket+entry);        //add the bucket to the head
  else bucket[entry] = p;                //else eliminate it

  return p;
}

void*
HashTable:: lookup(unsigned* key)
{

  int entry = hash (key);
  
  HashEntryPtr p = searchBucket(
                    *(bucket+entry),
                    key);
  

  if(!p) 
    return NULL;                      //if not found, return 0

  return p->value;                    //if found return a pointer  
}

void 
HashTable:: add (unsigned* key, void* value){

  int entry = hash(key);

  HashEntryPtr p = searchBucket (*(bucket+entry), key);
  if (p == NULL){  //was (!p)
     p = addElement(entry, key);
     p->value = value;
  }

  return;
}

void HashTable:: remove (unsigned* key)
{
  //  printf(" removing an unknown object with key %u %u\n",key[0], key[1]);
  int entry = hash(key);

  HashEntryPtr p = searchBucket (*(bucket+entry), key);

  if(!p)
    return;

  if(*(bucket+entry) == p) {
    *(bucket+entry) = p->next;
    delete p;
  }

  else{
    if(!(p->next))
       delete p;    
    else{
      (p->pre)->next = p->next;
      (p->next)->pre = p->pre;
      delete p;
    }
  }

}
// for debugging...
void HashTable:: remove (unsigned* key, int whatflag)
{
/*  if(whatflag == 1) 
    printf(" removing an element with key %u %u\n",key[0], key[1]);
  else {
    printf(" removing a node with key %u %u\n",key[0], key[1]);
  }*/
  int entry = hash(key);

  HashEntryPtr p = searchBucket (*(bucket+entry), key);

  if(!p) {
    printf(" this object has already been deleted ****************************\n");
    return;
  }

  if(*(bucket+entry) == p) {
    *(bucket+entry) = p->next;
    delete p;
  }

  else{
    if(!(p->next))
       delete p;    
    else{
      (p->pre)->next = p->next;
      (p->next)->pre = p->pre;
      delete p;
    }
  }
}

int 
HashTable:: hash(unsigned* key)
{ 
  unsigned numer = *key;
  float    igaze = (float)numer/Range;
  int      igazee = igaze*NBUCKETS+.5;

  if(igazee>=NBUCKETS) igazee = NBUCKETS-1;

  return (igazee);
}  


void HashTable::print_out(int type)
{

  /*int count=0;
  char filename[6]="tab_x";
  if(type==0) filename[4]='n';
  else filename[4]='e';
    ofstream out2(filename, ios::out);
  out2<<"furcsa"<<'\n';
    for(int i=0; i<NBUCKETS; i++)
      {
      out2<<i<<"   ";
      if(bucket[i])
	{
	  out2<<bucket[i]->key[0]<<' '<<bucket[i]->key[1]<<"   ";
	  count++;
	  HashEntryPtr currentPtr=bucket[i];
	  while(currentPtr)
	    {	   	     
	      out2<<currentPtr->key[0]<<' '<<currentPtr->key[1]<<"   ";
	      count++;
	      currentPtr=currentPtr->next;
	      // if(currentPtr) 
	    }
	}
      out2<<'\n';
      }
    out2<<count;
    out2.close();
  */
}


HashEntryPtr* HashTable::getbucketptr()
{
  return bucket;
}

int HashTable::get_no_of_buckets()
{
  return NBUCKETS;
}

double* HashTable::get_Xrange()
{
  return &Xrange[0];
}

double* HashTable::get_Yrange()
{
  return &Yrange[0];
}



























































