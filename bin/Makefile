SYSOS   := $(shell uname)
CURFWD  := $(shell pwd)
DATE    := $(shell date)
WHO     := $(shell whoami)
HOST    := $(shell hostname)
SHELL   := /bin/sh
#BITS is also used in source.h to define whether code is 32 bit or 64 bit
#if code is 64 bit, must use BITS = -64
BITS = -64

include ../source.h

#First clear off, then specify
.SUFFIXES:
.SUFFIXES: .C .f .o

DIRNAME = $(notdir $(CURPWD))


#default config
CC          = mpicxx 
CCC         = mpicxx
F77         = mpif90
OPTFLAGS    = -g -O0 
CFLAGS      = $(OPTFLAGS)
FFLAGS      = $(OPTFLAGS)
INCLUDES    =
SYSLIBS	    = -lblas -lz 
LDFLAGS	    = $(OPTFLAGS)
LINKER	    = $(CCC) -o hpfem 
#/usr/bin/libtool --mode=link $(CCC) -o hpfem
CCC_SUFFIXES = 
F77_SUFFIXES = 


CCFLAGS     = $(CFLAGS)

CDIRLIST = $(addprefix :$(CSRCDIR),$(CUSERDIR))
FDIRLIST = $(addprefix :$(FSRCDIR),$(FUSERDIR))

VPATH = .$(CDIRLIST):$(FDIRLIST)

COBJ  = $(CSRC:.C=_$(SYSOS)_C.o)
cOBJ  = $(cSRC:.c=_$(SYSOS)_c.o)
FOBJ  = $(FSRC:.f=_$(SYSOS)_f.o)
RCOBJ = $(addprefix $(COBJDIR), $(COBJ))
RcOBJ = $(addprefix $(cOBJDIR), $(cOBJ))
RFOBJ = $(addprefix $(FOBJDIR), $(FOBJ))

export CCC
export CC
export F77
export CCFLAGS
export CFLAGS
export FFLAGS
export SYSOS
export COBJ
export cOBJ
export FOBJ
export VPATH
export BITS
export CCC_SUFFIXES 
export CC_SUFFIXES 
export F77_SUFFIXES
export INCLUDES

hpfem: Allobj
	$(LINKER) $(LDFLAGS) \
	$(RCOBJ) \
	$(RcOBJ) \
	$(RFOBJ) \
	$(SYSLIBS) 

Allobj:	
	(cd $(FOBJDIR); $(MAKE))
	(cd $(COBJDIR); $(MAKE))
	(cd $(cOBJDIR); $(MAKE))

cleanC:   
	rm -f hpfem
	rm -rf ../cobj/*_$(SYSOS)_C.o
cleanf:
	rm -rf ../fobj/*_$(SYSOS)_f.o  
cleanc:
	rm -rf ../ccobj/*_$(SYSOS)_c.o
clean:
	rm -f hpfem
	rm -rf ../cobj/*_$(SYSOS)_C.o
	rm -rf ../fobj/*_$(SYSOS)_f.o  
	rm -rf ../ccobj/*_$(SYSOS)_c.o
	rm -rf ../csrc/*/*~
	rm -rf ../fsrc/*/*~
