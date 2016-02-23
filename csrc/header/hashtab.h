/* Hash table */
/* Every table can process NBUCKETS entries */
/* All buckets with the same entry are linked together*/
/* The pointer to the first bucket is stored in the table */
/*---------------------------------------------------------*/

#ifndef HASHTABLE_H
#define HASHTABLE_H

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "constant.h"


struct HashEntry{
  unsigned   key[KEYLENGTH];  //key: object key word
  void*      value;   //value: poiter to record
  HashEntry* pre;
  HashEntry* next;    //pre, next: objects with same entry will be stored in a two-way link

  HashEntry(unsigned* keyi)
  {
    int i;
    for(i=0;i<KEYLENGTH; i++)
      key[i] = *(keyi+i);
    next   = NULL;
    pre    = NULL;
  }

  HashEntry()
    {
      value  = NULL;
      next   = NULL;
      pre    = NULL;
    }

  ~HashEntry(){         //keep the follower when deleting an object
     if(next)
       next->pre = pre;
     if(pre)
       pre->next = next;
     
  }

};


typedef HashEntry* HashEntryPtr;


class HashTable{

    friend int hash(unsigned* keyi);
    friend class Element;

  protected:
    unsigned  MinKey[2];
    unsigned  MaxKey[2];
    unsigned   Range;
    double   Xrange[2];
    double   Yrange[2];

    HashEntryPtr* bucket;
    int           NBUCKETS;
    int           PRIME;
    int           ENTRIES;
     
    int hash(unsigned* key);
    HashEntryPtr addElement(int entry, unsigned* key);
    HashEntryPtr searchBucket(HashEntryPtr p, unsigned* key);

  public:
    HashTable(unsigned*, unsigned*, int, int);
    HashTable(unsigned*, unsigned*, int, int, double* XR, double* YR);
    ~HashTable();
    
    void   add(unsigned* key, void* value);
    void*  lookup(unsigned* key);
    void   remove(unsigned* key);
    void   remove(unsigned* key, int);  //for debugging
    void   print_out(int);
    int get_no_of_buckets();
//    void   get_element_stiffness(HashTable*);
    HashEntryPtr*   getbucketptr();
    void* get_value();

    double* get_Xrange();
    double* get_Yrange();

    /*   double* getXrange();
    double*  getYrange();*/
    int get_no_of_entries();

};

#endif
