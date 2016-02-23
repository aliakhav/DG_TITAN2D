#ifdef SUNOS
extern "C" void dcopy_(int*, double*, int*, double*, int*);
extern "C" double ddot_(int*, double*, int*, double*, int*);
extern "C" void tri_(double*, int*, int*);
extern "C" void rhsub_(double*,double*, int*, int*, double*);
extern "C" void schur5_(double*, double*, double* , double*, int*, int*,
		       int*, int*, double*, double*, double*, double*);
#endif
#ifdef IBMSP
extern "C" void dcopy(int*, double*, int*, double*, int*);
extern "C" double ddot(int*, double*, int*, double*, int*);
extern "C" void tri(double*, int*, int*);
extern "C" void rhsub(double*,double*, int*, int*, double*);
extern "C" void schur5(double*, double*, double* , double*, int*, int*,
		       int*, int*, double*, double*, double*, double*);
#endif

#ifdef CRAY
extern "C" void DCOPY(int*, double*, int*, double*, int*);
extern "C" double DDOT(int*, double*, int*, double*, int*);
extern "C" void TRI(double*, int*, int*);
extern "C" void RHSUB(double*,double*, int*, int*, double*);
extern "C" void SCHUR5(double*, double*, double* , double*, int*, int*,
		       int*, int*, double*, double*, double*, double*);
#endif

