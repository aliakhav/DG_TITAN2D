#include<iostream.h>
#include<fstream.h>
#include<math.h>
#include<stdlib.h>
#include "element.h"
#include "node.h"
#include "boundary.h"
//load has to be applied on the middle node of the face!!!
// executable must include the number of processors on the command line -- e.g. ./preprocess 1    for 1 processor!

void Read_no_of_objects(int*, int*, int*, int*, int*, long*);
void Read_node_data(int*, Node*, long*);
void Read_element_data(int*, Node*, Element*, long*);
void Read_boundary_data(int*, int*, Node*, Boundary*, long*);
void Read_material_data(int, double*, double*, long*);
void Write_data(int, int, int, int, int, Node*, Element**, Boundary*, 
		unsigned*, unsigned*, double*, double*, double* , double*);
void Determine_neighbors(int, Element*, int, Node*);

const int material_length = 80;

int compare_key_fn(const void* elem1, const void* elem2)
{
  Element** em1 = (Element**) elem1;
  Element** em2 = (Element**) elem2;
  if(*((*em1)->pass_key()) < *((*em2)->pass_key()))
    return(-1);
  else if(*((*em1)->pass_key()) > *((*em2)->pass_key()))
    return(1);
  else if(*((*em1)->pass_key()+1) < *((*em2)->pass_key()+1))
    return(-1);
  else if(*((*em1)->pass_key()+1) > *((*em2)->pass_key()+1))
    return(1);
  else if(*((*em1)->pass_key()+1) ==  *((*em2)->pass_key()+1))
    return(0);
      
  cout<<"something wrong in the qsort function compare key!\n";
  return(0);
}

int main(int argc, char** argv)
{
 int node_count, i;
 int element_count;
 int force_count;
 int constraint_count;
 int material_count;
 long location;
 unsigned minkey[2]={0, 0};
 unsigned maxkey[2]={0, 0};
 int p;

 if(argc != 2) {
   p = 1;
   cout<<"Need to input the number of processors as a command line argument if you want more than 1 processor!"<<endl<<flush;
 }

 Read_no_of_objects(&node_count, &element_count, &force_count, 
		    &constraint_count, &material_count, &location);

 Node* node=new Node[node_count];
 Element* element=new Element[element_count];
 Boundary* boundary=new Boundary[force_count+constraint_count];
 Element** ordering=new Element*[element_count];
 double* lambda = new double[material_count];
 double* mu = new double[material_count];

 Read_node_data(&node_count, node, &location);
 Read_element_data(&element_count, node, element, &location);
 Read_boundary_data(&force_count, &constraint_count, node, boundary, &location);
 Read_material_data(material_count, lambda, mu, &location);

 /* for(int i=0; i<element_count; i++)
    element[i].case5();*/

 double max[2]={0, 0};
 double min[2]={0, 0};

 element[0].create_m_node(max, min);

 min[0]=max[0]=*((*(element[0].get_element_node()))->get_node_coord());
 min[1]=max[1]=*((*(element[0].get_element_node()))->get_node_coord()+1);

 for(i=0; i<element_count; i++)
   element[i].create_m_node(max, min);

 unsigned nkey=2;

 for(i=0; i<element_count; i++)
   element[i].find_boundary(constraint_count, force_count, boundary);

 Determine_neighbors(element_count, element, node_count, node);

 for(i=0; i<node_count; i++)
   node[i].determine_max_min(max, min);

 node[0].determine_the_key(nkey, max, min, maxkey, minkey);
 for(i=0; i<2; i++)
   maxkey[i]=minkey[i]=*(node[0].get_key()+i);

 for(i=1; i<node_count; i++)
   node[i].determine_the_key(nkey, max, min, maxkey, minkey);

 for(i=0; i<element_count; i++)
   (*(element[i].get_element_node()+8))->determine_the_key(nkey, max, min, maxkey, minkey);

 for(i=0;i<element_count;i++)
   ordering[i] = &(element[i]);

 // before doing qsort, switch the first element with the middle element
 Element* EmTemp = ordering[0];
 ordering[0] = ordering[element_count/2];
 ordering[element_count/2] = EmTemp;

 qsort(ordering, element_count, sizeof(Element*), compare_key_fn);
 
/* cout<<"Number of processors: ";
 cin>>p;
 while(p <=0)
   {
     cout<<"Number of processors must be greater than 0\n";
     cout<<"Number of processors: ";
     cin>>p;
   } */
 if(argc == 2)
   p = atoi(argv[1]); // the number of processors -- this parameter is passed in
 
 for(i=0; i<element_count; i++)
  ( ordering[i])->myproc(p, i, element_count);
 
 Write_data(p, node_count, element_count, (force_count+constraint_count), material_count, 
	    node, ordering, boundary, maxkey, minkey, min, max, lambda, mu);
 

 delete []lambda;
 delete []mu;
 delete []node;
 delete []element;
 delete []boundary;
 delete []ordering;
 return(0);
 
}




//*****************************FUNCTIONS*******************************

//****************************DATA READ IN*****************************

void Read_no_of_objects(int* nc, int* ec, int* fc, int* cc, int* mc, long* loc)
{
  char endline;;
  ifstream inDatafile("funky.dat", ios::in);
  if(inDatafile.fail()) cout<<"file not found\n";
  inDatafile>>*nc;
  endline = '1';
  while(endline != '\n')
    inDatafile.get(endline);
  inDatafile>>*ec;
  endline = '1';
  while(endline != '\n')
    inDatafile.get(endline);
  inDatafile>>*cc;
  endline = '1';
  while(endline != '\n')
    inDatafile.get(endline);
  inDatafile>>*fc;
  endline = '1';
  while(endline != '\n')
    inDatafile.get(endline);
  inDatafile>>*mc;
  endline = '1';
  while(endline != '\n')
    inDatafile.get(endline);

  *loc=inDatafile.tellg();
  inDatafile.close();
  // cout<<"done"<<flush;

}


void Read_node_data(int* nc, Node n[], long* loc)
{

  int id;
  double node_coordinates[2];

  ifstream inDatafile("funky.dat", ios::in);
  if(!inDatafile) cout<<"file not found\n";
  inDatafile.seekg(*loc);
  for(int i=0; i<*nc; i++)
    {
      inDatafile>>id;
      inDatafile>>node_coordinates[0]; 
      inDatafile>>node_coordinates[1]; 

      n[i].setparameters(id, node_coordinates);

    }
  *loc=inDatafile.tellg();
  inDatafile.close();
}


void Read_element_data(int* ec, Node n[], Element e[], long* loc )

{

  int id;
  int element_nodes[8];
  int material;
  int elm_loc[2];

  Node* address[8];

  ifstream inDatafile("funky.dat", ios::in);
  if(!inDatafile) cout<<"file not found\n"<<flush;
  inDatafile.seekg(*loc);

  for(int i=0; i<*ec; i++)
    {
      inDatafile>>id;
      for(int j=0; j<8; j++)
       inDatafile>>element_nodes[j];

      for(int ii=0; ii<8; ii++)
	{
	  int w=0;
 
	  w = element_nodes[ii] - 1;
	  if(n[w].get_nodeid()!=element_nodes[ii]) {
	    w = 0;
	    while(n[w].get_nodeid()!=element_nodes[ii]) 
	      w++;
	    cout<<"found the node the long way\n"<<flush;
	  }
	  address[ii]=&n[w];

	}

      inDatafile>>material;
      inDatafile>>elm_loc[0];
      inDatafile>>elm_loc[1];
      e[i].setparameters(id, address, material-1, elm_loc);


    }
  *loc=inDatafile.tellg();
  inDatafile.close();
}

void Read_boundary_data(int* fc, int* cc, Node n[], Boundary b[], long* loc )

{

  int id, i;
  double xcomp;
  double ycomp;

  ifstream inDatafile("funky.dat", ios::in);
  if(!inDatafile) cout<<"file not found\n"<<flush;
  inDatafile.seekg(*loc);

  for(i=0; i<*cc; i++) //first constraints, then forces
    {
      inDatafile>>id;
      inDatafile>>xcomp;
      inDatafile>>ycomp;
      int w=0;
      while(n[w].get_nodeid()!=id) w++;
      b[i].setparameters(&n[w], xcomp, ycomp, -2);
    }

  for(i=*cc; i<(*fc+*cc); i++)
    {
      inDatafile>>id;
      inDatafile>>xcomp;
      inDatafile>>ycomp;
      int w=0;
      while(n[w].get_nodeid()!=id) w++;
      b[i].setparameters(&n[w], xcomp, ycomp, -3);
    }

  *loc=inDatafile.tellg();
  inDatafile.close();
}

void Read_material_data(int material_count, double* lambda, double* mu, long* loc)
{
  int i;
  char material_name[material_length];
  ifstream inDatafile("funky.dat", ios::in);
  if(inDatafile.fail()) cout<<"file not found\n";
  inDatafile.seekg(*loc);
  
  char endline = '1';
  while(endline != '\n')
    inDatafile.get(endline);

  for(i=0;i<material_count;i++)
    {
      inDatafile.getline(material_name, material_length);
      ifstream inD2("frict.data", ios::in);
      inD2 >> lambda[i];
      inD2 >> mu[i];
      inD2.close();
/*      cout<<"Input material properties for material "<<i+1<<endl;
      cout<<"Input internal friction angle for "<<material_name<<endl;
      cin>>lambda[i];
      while(lambda[i] <= 0)
	{
	  cout<<"Lambda must be greater than 0\n";
	  cout<<"Input lambda for "<<material_name<<endl;
	  cin>>lambda[i];	 
	} 
      cout<<"Input bed friction angle for "<<material_name<<endl;
      cin>>mu[i];
      while(mu[i] <= 0)
	{
	  cout<<"Mu must be greater than 0\n";
	  cout<<"Input mu for "<<material_name<<endl;
	  cin>>lambda[i];	 
	} */
    }

  *loc=inDatafile.tellg();
  inDatafile.close();
  // cout<<"done"<<flush;

}


//**************************FINDING THE NEIGHBORS******************************
void Determine_neighbors(int element_count, Element* element, int node_count, Node* node) 
{
  int i, j;
  for(i=0;i<node_count;i++)
    node[i].put_element_array_loc(-1);

  for(i=0;i<element_count;i++)
    for(j=4;j<8;j++)
      (*(element[i].get_element_node()+j))->put_element_array_loc(i);
        
  for(i=0; i<element_count; i++)
    element[i].determine_neighbors(i, element);

  for(i=0;i<element_count;i++)
    element[i].determine_opposite_brother();


  return;
}

//**************************DATA OUTPUT******************************

void Write_data(int np, int nc, int ec, int bc, int mc, Node n[], 
		Element* o[], Boundary b[], unsigned maxk[], 
		unsigned mink[], double min[], double max[], 
		double* lambda, double* mu)
{
 char filename[14]="funkyxxxx.inp";
 int el_per_proc=ec/np; 

 int i, j;
 for(i=0; i<np; i++)
   {
     int subdomain_nodes=0;
     int written=0;
     filename[5] = 48 + i/1000;  
     filename[6] = 48 + (i%1000)/100;
     filename[7] = 48 + (i%100)/10;  
     filename[8] = 48 + i%10;

     
     ofstream outDatafile(filename, ios::out);
     if(!outDatafile) cout<<"Could not be created!!!"<<'\n';
     int x, c;
     
     c=i*el_per_proc;
     if(i!=np-1) x=(i+1)*el_per_proc;
     else x=ec;
     
     if(i>0)
       for(j=0; j<nc; j++)
	 n[j].clear_written_flag();
     
     for(j=c; j<x; j++)
       {
	 
	 for(int k=0; k<9; k++)
	   {
	     
	     written=(*((o[j])->get_element_node()+k))->get_written_flag();
	     if(written==0)
	       {
		 subdomain_nodes++;
		 (*((o[j])->get_element_node()+k))->set_written_flag();
	       }
	     
	   }
	 
       }
     

     outDatafile<<subdomain_nodes<<' '<<mink[0]<<' '<<mink[1]<<' '<<' '<<maxk[0]<<' '<<maxk[1]<<'\n';
     
     outDatafile<<min[0]<<" "<<max[0]<<' '<<min[1]<<' '<<max[1]<<endl;
     
     for(j=0; j<nc; j++)
       n[j].clear_written_flag();
     
     /*    for(j=0; j<nc; j++)
	  n[j].write_node_data(&outDatafile);*/
     
     for(j=c; j<x; j++)
       {
	 
	 for(int k=0; k<9; k++)
	   {
	     if(k==8)
	       (*((o[j])->get_element_node()+k))->clear_written_flag();
	     (*((o[j])->get_element_node()+k))->write_node_data(&outDatafile);      	
	  }
	 
       }
     
     //element data start here

     outDatafile<<x-c<<' '<<*((*((o[c])->get_element_node()+8))->get_key())<<' '<<*((*((o[c])->get_element_node()+8))->get_key()+1)<<'\n';
     
     outDatafile<<*((*((o[x-1])->get_element_node()+8))->get_key())<<' '<<*((*((o[x-1])->get_element_node()+8))->get_key()+1)<<'\n';
     

     for(j=c; j<x; j++)
       {
	 (o[j])->write_element_data(&outDatafile);
	 outDatafile<<'\n';
       }

     //material properties start here
     outDatafile<<mc<<endl;
     for(j=0;j<mc;j++)
       outDatafile<<lambda[j]<<"  "<<mu[j]<<endl;
     
    /*
    outDatafile<<'\n';

    outDatafile<<bc<<'\n';

    for(int k=0; k<bc; k++)
      b[k].write_b_data(&outDatafile);
    outDatafile<<'\n';


    outDatafile.close();*/
    
   }

}




