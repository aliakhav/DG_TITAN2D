#ifndef INTERFACE
#define INTERFACE
struct InterFace{
  unsigned      key[KEYLENGTH]; /*--key of elem on interface--*/
  InterFace*    pre;
  InterFace*    next;
  
  InterFace(unsigned* keyi)
  {
    int i;
    for(i=0;i<KEYLENGTH; i++)
      key[i] = *(keyi+i);
    next   = NULL;
    pre    = NULL;
  }
  InterFace(){
    next   = NULL;
    pre    = NULL;
  }
  ~InterFace(){
     if(next)
       next->pre = pre;
     if(pre)
       pre->next = next;
  }
};

typedef InterFace* InterFacePtr;

struct dof{ 
  int value;
  dof* next;

  dof(int fod)
    {
      value=fod;
      next = 0;

    }
  /*  ~dof()
    {
      dof* temp;

      while(current)
	{
	  temp=current;
	  current = current->next;
	  delete temp;
	}
      
	}*/


};

typedef dof* dofPtr;

struct NNLink{
  unsigned   key[KEYLENGTH]; 
  NNLink*    pre;
  NNLink*    next;
  
  NNLink(unsigned* keyi)
    {
      int i;
      for(i=0;i<KEYLENGTH; i++)
	key[i] = *(keyi+i);
      next   = NULL;
      pre    = NULL;
    }
  NNLink(){
    next   = NULL;
    pre    = NULL;
  }
  ~NNLink(){
    if(next)
      next->pre = pre;
    if(pre)
      pre->next = next;
  }
};

typedef NNLink* NNPtr;

struct SSLink{
  unsigned      key[KEYLENGTH]; 
  SSLink*    pre;
  SSLink*    next;
  
  SSLink(unsigned* keyi)
    {
      int i;
      for(i=0;i<KEYLENGTH; i++)
	key[i] = *(keyi+i);
      next   = NULL;
      pre    = NULL;
    }
  SSLink(){
    next   = NULL;
    pre    = NULL;
  }
  ~SSLink(){
    if(next)
      next->pre = pre;
    if(pre)
      pre->next = next;
  }
};

typedef SSLink*  SSPtr;

struct FFLink{
  unsigned      key[KEYLENGTH]; 
  FFLink*    pre;
  FFLink*    next;
  
  FFLink(unsigned* keyi)
    {
      int i;
      for(i=0;i<KEYLENGTH; i++)
	key[i] = *(keyi+i);
      next   = NULL;
      pre    = NULL;
    }
  FFLink(){
    next   = NULL;
    pre    = NULL;
  }
  ~FFLink(){
    if(next)
      next->pre = pre;
    if(pre)
      pre->next = next;
  }
};

typedef FFLink* FFPtr;

struct VVLink{
  unsigned      key[KEYLENGTH]; 
  VVLink*    pre;
  VVLink*    next;
  
  VVLink(unsigned* keyi)
    {
      int i;
      for(i=0;i<KEYLENGTH; i++)
	key[i] = *(keyi+i);
      next   = NULL;
      pre    = NULL;
    }
  VVLink(){
    next   = NULL;
    pre    = NULL;
  }
  ~VVLink(){
    if(next)
      next->pre = pre;
    if(pre)
      pre->next = next;
  }
};

typedef VVLink* VVPtr;

struct EELink{
  unsigned      key[KEYLENGTH]; 
  EELink*    pre;
  EELink*    next;
  
  EELink(unsigned* keyi)
    {
      int i;
      for(i=0;i<KEYLENGTH; i++)
	key[i] = *(keyi+i);
      next   = NULL;
      pre    = NULL;
    }
  EELink(){
    next   = NULL;
    pre    = NULL;
  }
  ~EELink(){
    if(next)
      next->pre = pre;
    if(pre)
      pre->next = next;
  }
};

typedef EELink* EEPtr;

struct BBLink{
  unsigned      key[KEYLENGTH]; 
  BBLink*    pre;
  BBLink*    next;
  
  BBLink(unsigned* keyi)
    {
      int i;
      for(i=0;i<KEYLENGTH; i++)
	key[i] = *(keyi+i);
      next   = NULL;
      pre    = NULL;
    }
  BBLink(){
    next   = NULL;
    pre    = NULL;
  }
  ~BBLink(){
    if(next)
      next->pre = pre;
    if(pre)
      pre->next = next;
  }
};

typedef BBLink* BBPtr;

struct CCLink{
  unsigned      key[KEYLENGTH]; 
  CCLink*    pre;
  CCLink*    next;
  
  CCLink(unsigned* keyi)
    {
      int i;
      for(i=0;i<KEYLENGTH; i++)
	key[i] = *(keyi+i);
      next   = NULL;
      pre    = NULL;
    }
  CCLink(){
    next   = NULL;
    pre    = NULL;
  }
  ~CCLink(){
    if(next)
      next->pre = pre;
    if(pre)
      pre->next = next;
  }
};

#endif
