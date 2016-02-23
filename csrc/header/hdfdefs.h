#define DATASETNAME "ExtendibleArray1" 
#define DATASETNAME1 "ExtendibleArray2" 
#define RANK        3
#define LENGTH 1  
#define HEIGHT 10000
#define WIDTH 10
#define HEIGHT1    2500
#define aHEIGHT 10
#define aWIDTH 4
#define aHEIGHT1    2

float      data[1][HEIGHT][WIDTH];
float      data1[1][HEIGHT1][WIDTH];
float      adata[1][aHEIGHT][aWIDTH];
float      adata1[1][aHEIGHT1][aWIDTH];
hsize_t    maxdim[3] = {LENGTH,H5S_UNLIMITED, WIDTH}; 
hsize_t    amaxdim[3] = {LENGTH,H5S_UNLIMITED,aWIDTH};
hsize_t    dim[3] = {1, HEIGHT, WIDTH};   
hsize_t    adim[3] = {1, aHEIGHT, aWIDTH};   
hid_t      cparms,acparms;
hsize_t    newsize[3]={1, HEIGHT, WIDTH};
hsize_t    anewsize[3]={1, aHEIGHT, aWIDTH};
void       *tbuf = NULL;
  

int w[100],aw[100],r[100],ar[100],l[100],al[100],al1[100],l1[100];
