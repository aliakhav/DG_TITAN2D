#include "../header/hpfem.h"

// compare key compares 2 keys and returns 0 if they are not the same and 1 if they are
int compare_key(unsigned* key1, unsigned* key2)
{
  int i;
 
  for ( i = 0; i < KEYLENGTH; i++ )
    if ( *(key1 + i) != *(key2 + i) )
      return 0;

  
  return 1;
}
