SHELL = /bin/sh 
SYSOS := $(shell uname)
HOST := $(shell hostname)

################################
# system specific information
################################

# for the linux boxes in the trailer...
ifeq "$(HOST)" "popo" 
override HOST = Trailer
endif
ifeq "$(HOST)" "colima" 
override HOST = Trailer
endif
ifeq "$(HOST)" "elchichon" 
override HOST = Trailer
endif
ifeq "$(HOST)" "rainier"
override HOST = Trailer
endif
ifeq "$(HOST)" "sthelens"
override HOST = Trailer
endif
ifeq "$(HOST)" "casita" 
override HOST = Trailer
endif

# Linux OS in the trailer
ifeq "$(HOST)" "Trailer"
libgrasslibdir = /usr/local/lib
libgrassincdir = /usr/local/include

CC = gcc 
LIBTOOL = /usr/bin/libtool 
LIBS = -lnsl -lz -lm -lg2c -lc 
INCLUDES = -I. -I$(libgrassincdir) -I../../csrc/header
DEFS = -DHAVE_CONFIG_H 
CCLD = $(CC) 
endif
#Linux standalone
ifeq "$(HOST)" "localhost"
libgrasslibdir = /usr/local/lib
libgrassincdir = /usr/local/include

CC = gcc
LIBTOOL = /usr/bin/libtool
LIBS = -lnsl -lz -lm -lg2c -lc
INCLUDES = -I. -I$(libgrassincdir)
DEFS = -DHAVE_CONFIG_H
CCLD = $(CC)
endif

# nash.ccr.buffalo.edu
ifeq "$(HOST)" "nash.ccr.buffalo.edu"
libgrasslibdir = /oldutil/grass/libgrass5/lib
libgrassincdir = /oldutil/grass/libgrass5/include

CC = gcc
LIBTOOL = /usr/bin/libtool 
LIBS = -lnsl -lz -lm -lg2c -lc 
INCLUDES = -I. -I$(libgrassincdir) -I../../csrc/header
DEFS = -DHAVE_CONFIG_H 
CCLD = $(CC)
endif

# joplin.ccr.buffalo.edu
ifeq "$(HOST)" "joplin.ccr.buffalo.edu"
libgrasslibdir = /util/grass/libgrass5-1.0.0/lib
libgrassincdir = /util/grass/libgrass5-1.0.0/include

CC = gcc
LIBTOOL = /usr/bin/libtool 
LIBS = -lnsl -lz -lm -lg2c -lc 
INCLUDES = -I. -I$(libgrassincdir) -I../../csrc/header
DEFS = -DHAVE_CONFIG_H 
CCLD = $(CC)
endif

# crosby.ccr.buffalo.edu
ifeq "$(HOST)" "crosby.ccr.buffalo.edu"
libgrasslibdir = /util/grass/libgrass5/lib
libgrassincdir = /util/grass/libgrass5/include

CC = cc -w
LIBTOOL = /usr/freeware/bin/libtool 
LIBS = -lnsl -lz -lm -lc 
INCLUDES = -I. -I$(libgrassincdir) -I../../csrc/header
DEFS = -DHAVE_CONFIG_H 
CCLD = $(F77) 
endif

################################
# everything below shouldn't
# need to be changed
################################
OPT_FLAGS = -O2
GrassApi_SOURCES = create_funky.c 
GrassApi_LDADD = $(libgrasslibdir)/libgrass5.la 

CPPFLAGS = 
LDFLAGS = 
GrassApi_OBJECTS = create_funky.o GrassApi.o  preproc.o
GrassApi_DEPENDENCIES = $(libgrasslibdir)/libgrass5.la 
GrassApi_LDFLAGS = 
CFLAGS = -O3
COMPILE = $(CC) $(DEFS) $(INCLUDES) $(CPPFLAGS) $(CFLAGS) 
LTCOMPILE = $(CC) $(DEFS) $(INCLUDES) $(CPPFLAGS) $(CFLAGS) 
LINK = $(LIBTOOL) --mode=link $(CCLD) $(CFLAGS) $(LDFLAGS) -o $@ 

SOURCES = $(GrassApi_SOURCES) 
OBJECTS = $(GrassApi_OBJECTS) 

.c.o: 
	$(COMPILE) $(OPT_FLAGS) -c $< 

.f.o:
	$(F77) $(OPT_FLAGS) -c $<

preproc.x: $(GrassApi_OBJECTS) $(GrassApi_DEPENDENCIES) 
	@rm -f 
	$(LINK) $(GrassApi_LDFLAGS) $(GrassApi_OBJECTS) $(GrassApi_LDADD) $(LIBS) 
	@rm -rf .libs

GisApi.o: ../../csrc/geoflow/GisApi.C ../../csrc/header/GisApi.h
	 $(COMPILE) $(OPT_FLAGS) -c ../../csrc/geoflow/GisApi.C

create_funky.o: create_funky.c GisApi.o ../../csrc/header/GisApi.h preproc.o

preproc.o: preproc.f

clean:
	rm -rf *.o *~ preproc.x .libs pile.data frict.data funky.dat
