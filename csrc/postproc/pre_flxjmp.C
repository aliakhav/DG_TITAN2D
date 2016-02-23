#include"../header/hpfem.h"
#include"../header/post.h"

#define NODE_NUM 8   //for 3d


extern "C" void flxjmp_(int* matid,double xy[][2],double neig_Xnod[][18],int Neig[][2],int* Nside,int* Norder,double ssol[][121],double tsol[][242],int* nequ,double* Rhs,int  norder[][5]);


 //extern "C" void flxerrelt_(double xnod[][2],int *norder,double sol[][121],double* fjmp,int *matvalue,int* bcvalue,int* nelb,int* icon,double* err,double* elsol,int* nequ);

extern void send_neigh_sol(HashTable* ht_elem_ptr,HashTable* ht_node_ptr,Neigh_Sol* neigh_list,Neigh_Sol_Pack** rec_neigh_sol,int *elem_num,int* rec_elem_num);


extern Neigh_Sol* prep_neigh_list(HashTable* ht_elem_ptr,HashTable* ht_node_ptr,Neigh_Sol* neigh_list,int* elem_num,double* u,int* ix);


void add_sol(int* assocP,int k,Neigh_Sol_Pack** rec_neigh_sol,int* rec_elem_num,int neigh[][2],double tsol [][242],int* neigh_gen,int mygen,int myid,unsigned* neigh_keyP,int* nside,int** neigh_order,double neig_Xnod[][18]);   //to add solutions from neigbors in different proc.

extern double* elem_sol(double* u, int* ix, HashTable* HT_Node_Ptr, HashTable* HT_Elem_Ptr,int myid, int numprocs,unsigned* key_el, double* Utemp);

void get_el_cord(HashTable* ht_node_ptr,Element* EmTemp,double* Xnod);


void pre_flex(HashTable* ht_elem_ptr,HashTable* ht_node_ptr,double* u,int* ix)
  // to compute the error we need to calculate the fluxjmp, we require the neighbor solution,neighbor side on which the current element is,information about the neigh gen. this code computes all this information
{

  int            i,j,k,mi,l;
  int            ii,jj;
  int            flag;
  int            myid ,numprocs;
  Element*       EmTemp;
  Element*       neighbor;
  Node*          NdTemp;
  unsigned       KeyTemp[KEYLENGTH];
  HashEntryPtr   entryp;
  HashEntry**    temp_entryp;
  unsigned*      keyP;
  void*          p;
  int*           assocP;
  int*           neigh_gen;
  int            mygen; 
  int            nside[4];          //indicates side of element on which the neighbor is situated

  int*           neigh_order[8];   //to store the neighbor order
  
  double         Xnod[9][2]; 
  double*        coord;
  
  double         tsol[8][242];
  unsigned*      neigh_keyP=NULL;
  unsigned*      tempkey =NULL;
  int*           norder;

  double         Utemp[242];    //used to store temporary element solution


  double neig_Xnod[8][18];



  

  MPI_Status status;
  MPI_Request request;
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  char  filename[11] = "debugx.txt";
  // filename[7] = myid/10;
  filename[7] = 48 + myid%10;
 FILE *fpt; 
  if((fpt = fopen(filename,"w"))==NULL)
printf("error opening file\n");


  //if the neighbor is on  a different proc,solution needs to be communicated between processors. the strategy used is, on each subdomain we travel along the interface  make a list of sol that needs to be communicated to each subdomain and then do the interprocessor communication
  Neigh_Sol*  neigh_list;      
  int*        elem_num;  //keeps count of no of element to be sent  to each proc.

  neigh_list = new Neigh_Sol[numprocs];
  elem_num = new int[numprocs];

  neigh_list = prep_neigh_list(ht_elem_ptr,ht_node_ptr,neigh_list,elem_num,u,ix); //prepares a list of neigh sol that need to be communicated
  MPI_Barrier(MPI_COMM_WORLD); 

 //datastruct to receive neigh solution from all procs
  int* rec_elem_num;
  rec_elem_num = new int[numprocs];

  Neigh_Sol_Pack** rec_neigh_sol;
  rec_neigh_sol = new Neigh_Sol_Pack*[numprocs];

  send_neigh_sol(ht_elem_ptr,ht_node_ptr,neigh_list,rec_neigh_sol,elem_num,rec_elem_num); //module which manages the interprocessor communication



int elements = ht_elem_ptr->get_no_of_buckets(); 
temp_entryp = ht_elem_ptr->getbucketptr();
  for(i=0;i<elements;i++)
    { 
      entryp = *(temp_entryp + i);
      while(entryp)
	{
           EmTemp = (Element*)(entryp->value);
           for(ii=0;ii<4;ii++)
             nside[ii] =-1;
           for(ii=0;ii<8;ii++)
             neigh_order[ii] =NULL;
          

	   if(!EmTemp->get_refined_flag())
	     {
      
               fprintf(fpt,"\n \n the element is %u  %u  %d",*(EmTemp->pass_key                       ()+0),*(EmTemp->pass_key()+1),myid);       

	       //get the node coord
                 get_el_cord(ht_node_ptr,EmTemp,Xnod[0]);

	        //neighbour information

                 assocP    = EmTemp->getassoc();
                 neigh_gen = EmTemp->get_neigh_gen();
                 mygen     = EmTemp->get_gen();
                 neigh_keyP= EmTemp->get_neighbors();
                 keyP      = EmTemp->pass_key();
                 norder    = EmTemp->get_order();

		 

                for(int ii=0;ii<8;ii++)
                  for(int jj=0;jj<121;jj++)
                     tsol[ii][jj] = 0.0;
                
             
               int  neigh[4][2] = {0,0,0,0,0,0,0,0}; //stores informaion about neighbor generation
               for(k=0;k<4;k++)
	        {  
                 if(neigh_gen[k] < mygen && assocP[k]>=0)
		   {neigh[k][0] = -1;           //negative no if neighbor is a bigger element
                    neigh[k][1] = 0;
                   }
                 
                 if(neigh_gen[k]==mygen && assocP[k]>= 0)
		   {neigh[k][0] = 1;            //positive if of the same size
                   neigh[k][1] = 0;
                  } 
                 if(neigh_gen[k]> mygen && assocP[k]>=0)
                  {
		    neigh[k][0] = 1;             //positive in both  if neighbor is smaller   
                   neigh[k][1] = 1;
                  }
                }
           
              int p_flag = 0;
              int flag =0;       
              int kk   = 0;

              //traverse the element sides
              for(k=0;k<4;k++)
		{
                  p_flag=0;

                  if(myid==assocP[k])
                   neighbor =(Element*) ht_elem_ptr->lookup(&neigh_keyP[k*KEYLENGTH]);
                  else
                   p_flag =1;  


                 if(p_flag==0)
		   {
                    tempkey =  neighbor->get_neighbors();

                    for(mi=0;mi<8;mi++)
		      { 
                       kk=0;
                       for(l=0;l<KEYLENGTH;l++)
                          {
                         if(tempkey[mi*KEYLENGTH + l]  ==*(EmTemp->pass_key()+l))kk++;
	        	  }

                       if(kk==KEYLENGTH)break;
                      }

		 assert(mi<8);
		 nside[k] = mi>3 ? mi = mi-4:mi=mi;
                
                   
                 elem_sol(u,ix,ht_node_ptr,ht_elem_ptr,myid,numprocs,neighbor->pass_key(),Utemp);  
                    for(int ii=0;ii<2;ii++)
                       for(int jj=0;jj<121;jj++)
                         {
                           tsol[k][ii*121 + jj] = Utemp[ii*121 +jj];
                         }
                     
                  neigh_order[k] = neighbor->get_order();
                  get_el_cord(ht_node_ptr,neighbor,neig_Xnod[k]);
                
                   

                   if(neigh_gen[k]> mygen)
		      {  neighbor =(Element*) ht_elem_ptr->lookup(&neigh_keyP[k+4]); 
                      elem_sol(u,ix,ht_node_ptr,ht_elem_ptr,myid,numprocs,neighbor->pass_key(),Utemp);  
                      for(ii=0;ii<2;ii++)
                          for(jj=0;jj<121;jj++)
                            {
		                  tsol[k+4][ii*121 + jj] = Utemp[ii*121 +jj];
                                  
                            }
                         neigh_order[k+4] = neighbor->get_order();
                         get_el_cord(ht_node_ptr,neighbor,neig_Xnod[k+4]);

                              

                      }// end of   if(neigh_gen[k]> mygen) 

                  
		   }//end of if(p_flag=0) 
                  else
		    {   
                   
                     add_sol(assocP,k,rec_neigh_sol,rec_elem_num,neigh,tsol,neigh_gen,mygen,myid,neigh_keyP,nside,neigh_order,neig_Xnod);
                    }//end of else
               

                }//end of traversing the element sides

               double Rhs[242];
               int nequ = EQUATIONS;
               double ssol[2][121];
               int material = EmTemp->get_material();


             


               for(ii=0;ii<8;ii++)
                 for(jj=0;jj<2;jj++)
                   {
                fprintf(fpt,"\n the neighbor key  is %u ",neigh_keyP[ii*2 + jj]);
                   }


                elem_sol(u,ix,ht_node_ptr,ht_elem_ptr,myid,numprocs,EmTemp->pass_key(),ssol[0]);  

               for(ii=0;ii<2;ii++)
                 for(jj=0;jj<121;jj++)
                     {
	                  fprintf(fpt,"\n the element  sol is %f %d %d ",ssol[ii][jj],ii,jj);
                     }


  
              for(ii=0;ii<8;ii++)
                  for(int jj=0;jj<242;jj++)
                   {
                      fprintf(fpt,"\n the element  neigh sol is %f %d  %d ",tsol[ii][jj],ii,jj);
                   }


              

            

		 // neighbor orders
	int Neigh_poly_order[8][5];
	for(ii=0;ii<8;ii++)
	   for (jj=0;jj<5;jj++)
	     {  if(neigh_order[ii]!=NULL)
                 Neigh_poly_order[ii][jj] = *(neigh_order[ii]+jj);
             }
	    //computing the flxjmp    
          
           flxjmp_(&material,Xnod,neig_Xnod,neigh,nside,norder,ssol,tsol,&nequ,Rhs,Neigh_poly_order);



           //computing the error

           int nelb[4];
           int nicon[4];
           int type[4]={0,0,0,0};
  
            if(EmTemp->get_bcptr()->type)
	      {
            for(ii=0;ii<4;ii++)
               type[ii]= *(EmTemp->get_bcptr()->type+ii);
              }

           EmTemp->get_nelb_icon(ht_node_ptr,ht_elem_ptr,nelb,nicon);

         
//flxerrelt_(Xnod,norder,ssol,Rhs,&material,type,nelb,nicon,EmTemp->get_el_error(),EmTemp->get_el_solution(),&nequ);  
            
           }//end of if(!EmTemp->get_refined_flag())

            entryp =entryp->next;
	}//end of while

    }//end of outer for

   fclose(fpt);   
  if (myid==0){
    int w;

    printf("type ina number\n\r");
    scanf("%d", &w);
  } 

 MPI_Barrier(MPI_COMM_WORLD); 

}//end of func
















void add_sol(int* assocP,int k,Neigh_Sol_Pack** rec_neigh_sol,int* rec_elem_num,int neigh[][2],double tsol [][242],int* neigh_gen,int mygen,int myid,unsigned* neigh_keyP,int* nside,int **neigh_order,double neig_Xnod[][18])
    {  int i_el;  
       int kk;
       int ll;
       int ii,jj;

      if(assocP[k]!= -1 && assocP[k]!= -2 && assocP[k]!=myid)
           {
                       
              for(i_el=0;i_el<rec_elem_num[assocP[k]];i_el++)
                 { ll=0;
                   for( kk=0; kk<KEYLENGTH;kk++)
		    if(rec_neigh_sol[assocP[k]][i_el].key[kk]== neigh_keyP[k*KEYLENGTH + kk ]){ll++;}
                    
                   if(ll == KEYLENGTH)break;
              
                 }

              assert(i_el<rec_elem_num[assocP[k]]);              
              nside[k] = rec_neigh_sol[assocP[k]][i_el].nside;

             
             

              for(ii=0;ii<2;ii++)
                 for(jj=0;jj<121;jj++)
                   {
                     tsol[k][ii*121 + jj] =  rec_neigh_sol[assocP[k]][i_el].solu[ii][jj];
                   }
	        neigh_order[k] =  rec_neigh_sol[assocP[k]][i_el].norder;
                
                  for(ii=0;ii<18;ii++)
	         	{
                  neig_Xnod[k][ii] = rec_neigh_sol[assocP[k]][i_el].Xnod[ii];
	        	} 




         if(neigh_gen[k]> mygen)
	    {   

                for(i_el=0;i_el<rec_elem_num[assocP[k]];i_el++)
                 { ll=0;
                   for( kk=0; kk<KEYLENGTH;kk++)
		    if(rec_neigh_sol[assocP[k]][i_el].key[kk]== neigh_keyP[(k+4                                             )*KEYLENGTH + kk ]){ll++;}
                    
                   if(ll == KEYLENGTH)break;
              
                 }

              assert(i_el<rec_elem_num[assocP[k]]);
	      nside[k] = rec_neigh_sol[assocP[k]][i_el].nside;

              for(int ii=0;ii<2;ii++)
                 for(int jj=0;jj<121;jj++)
                   {
                     tsol[k +4][ii*121 + jj] =  rec_neigh_sol[assocP[k]][i_el].solu[ii][jj];
                   }

               neigh_order[k+4] =  rec_neigh_sol[assocP[k]][i_el].norder;
             
                   for(ii=0;ii<18;ii++)
	          	{
                  neig_Xnod[k+4][ii] = rec_neigh_sol[assocP[k]][i_el].Xnod[ii];
		       }    
                

             }//end of  if(neigh_gen[k]> mygen)
           


   }


}//end of func
 


void get_el_cord(HashTable* ht_node_ptr,Element* EmTemp,double* Xnod)
{
int i,j;
unsigned* keyP;
Node* NdTemp;
double* coord;

     keyP = EmTemp->getNode();
         for(j=0;j<8;j++)
           {  
	      NdTemp = (Node*) ht_node_ptr->lookup(keyP+j*KEYLENGTH);
	      coord = NdTemp->get_coord();
             for(i=0;i<DIMENSION;i++)
	       {
	        Xnod[j*DIMENSION + i] = *(coord +i);
               }
           }

           //bubblenode
           NdTemp = (Node*) ht_node_ptr->lookup(EmTemp->pass_key());
           coord = NdTemp->get_coord();


	      for(i=0;i<DIMENSION;i++)
	       {
	        Xnod[j*DIMENSION + i] = *(coord +i);
               }

}//end of func

