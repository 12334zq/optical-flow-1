# This is a -*- Makefile -*- 
# Warning: this file contains tabs which cannot be converted to spaces

# TOP_SRCDIR must be relative
TOP_SRCDIR = ../..
SRCDIR = .

#CONFIG=debug

#CC	=	$(CROSS)gcc
CXX	=	$(CROSS)g++
#CPP	=	$(CROSS)gcc -E
RANLIB	=	$(CROSS)ranlib
LD	=	$(CXX)

# CFLAGS_OPT is normally in the config$(CONFIG).mk file,
# but here's a reasonable default value
# add -ftree-vectorizer-verbose=3 to get information about vectorized loops
#CFLAGS_OPT=-Wall -g -O2 -ftree-vectorize -msse3 -mssse3 -ffast-math
#CFLAGS_OPT=-Wall -O3 -DNDEBUG -march=core2 -ftree-vectorize -msse3 -mssse3 -ffast-math
#CFLAGS_OPT=-Wall -g -O2 -march=pentium4 -ftree-vectorize -msse3 -mssse3 -ffast-math -DDISPARITYTAGGER_DISABLE_TEXTURE_RECTANGLE -DDISPARITYTAGGER_DISABLE_PBO
#CFLAGS_OPT=-Wall -g -O2 -march=pentium4 -ftree-vectorize -msse3 -mssse3 -ffast-math -DVIDEOINPUT_WITHOUT_NVSDI
CFLAGS_OPT= -Wall  -Wmissing-declarations
#CFLAGS_OPT= -Wall -g -O2
CFLAGS_OPENMP = -fopenmp

# the flag for .cpp files
#CXXFLAGS = -std=c++11

PNG_CPPFLAGS = $(shell pkg-config libpng --cflags)
PNG_LDFLAGS = $(shell pkg-config libpng --libs-only-L)
PNG_LIBS = $(shell pkg-config libpng --libs-only-l) $(shell pkg-config libpng --libs-only-other)
JPEG_CPPFLAGS =
JPEG_LDFLAGS =
JPEG_LIBS = -ljpeg
TIFF_CPPFLAGS = $(shell pkg-config libtiff-4 --cflags)
TIFF_LDFLAGS = $(shell pkg-config libtiff-4 --libs-only-L)
TIFF_LIBS = $(shell pkg-config libtiff-4 --libs-only-l) $(shell pkg-config libtiff-4 --libs-only-other)

LIBTOOL = libtool

-include config$(CONFIG).mk

CPPFLAGS= \
	$(PNG_CPPFLAGS) \
	$(JPEG_CPPFLAGS) \
	$(TIFF_CPPFLAGS)

LIBS = \
	$(PNG_LDFLAGS) $(PNG_LIBS) \
	$(JPEG_LDFLAGS) $(JPEG_LIBS) \
	$(TIFF_LDFLAGS) $(TIFF_LIBS)


CFLAGS_EXTRA=$(CFLAGS_OPENMP)
CXXFLAGS=$(CFLAGS_OPT) $(CFLAGS_EXTRA)
CFLAGS=$(CFLAGS_OPT) $(CFLAGS_EXTRA)
LDFLAGS=$(CFLAGS_OPT) $(CFLAGS_EXTRA) $(LDFLAGS_EXTRA) $(OPENCL_LDFLAGS)

# for gcov profiling add:
#-fprofile-arcs -ftest-coverage

COMPILE.c=$(CC) -c $(CFLAGS) $(CPPFLAGS)
COMPILE.cpp=$(CXX) -c $(CXXFLAGS) $(CPPFLAGS)

PROGRAMS = horn_schunck_classic horn_schunck_pyramidal brox_spatial brox_temporal tvl1flow tvl1occflow

LIBRARY = libof.a

#VPATH = $(SRCDIR):$(SRCDIR)/$(TOP_SRCDIR)/tools/log:$(SRCDIR)/$(TOP_SRCDIR)/tools/iniParser:$(SRCDIR)/$(TOP_SRCDIR)/tools/files

.PHONY: all lib

all: $(LIBRARY) $(PROGRAMS)

lib: $(LIBRARY)

of_SOURCES_CPP = \
	horn_schunck_classic.cpp \
	horn_schunck_pyramidal.cpp \
	brox_optic_flow_spatial.cpp \
	brox_optic_flow_temporal.cpp \
	brox_spatial_mask.cpp \
	brox_temporal_mask.cpp \
	tvl1flow.cpp \
	tvl1occflow.cpp \
	tvl1occflow_solvers.cpp \
	tvl1occflow_tv_rof_box.cpp \
	bicubic_interpolation.cpp \
	operators.cpp \
	utils.cpp \
	zoom.cpp

#	tvl1occflow.cpp \
#	tvl1occflow_solvers.cpp \
#	tvl1occflow_tv_rof_box.cpp \

of_SOURCES_C =
of_HEADERS = \
	horn_schunck.h \
	brox_optic_flow.h \
	brox_spatial_mask.h \
	tvl1flow.h \
	tvl1occflow.h \
	tvl1occflow_solvers.h \
	tvl1occflow_tv_rof_box.h \
	tvl1occflow_constants.h \
	bicubic_interpolation.h \
	operators.h \
	zoom.h \
	xmalloc.h

of_OBJS =  $(of_SOURCES_CPP:.cpp=.o) $(of_SOURCES_C:.c=.o)

horn_schunck_classic_SOURCES_CPP = \
	horn_schunck_classic_main.cpp \
	iio.cpp
horn_schunck_classic_SOURCES_C =
horn_schunck_classic_HEADERS = \
	iio.h
horn_schunck_classic_OBJS =  $(horn_schunck_classic_SOURCES_CPP:.cpp=.o) $(horn_schunck_classic_SOURCES_C:.c=.o)
horn_schunck_classic_LIBS = $(LIBRARY) $(LIBS)

horn_schunck_pyramidal_SOURCES_CPP = \
	horn_schunck_pyramidal_main.cpp \
	iio.cpp
horn_schunck_pyramidal_SOURCES_C =
horn_schunck_pyramidal_HEADERS = \
	iio.h
horn_schunck_pyramidal_OBJS =  $(horn_schunck_pyramidal_SOURCES_CPP:.cpp=.o) $(horn_schunck_pyramidal_SOURCES_C:.c=.o)
horn_schunck_pyramidal_LIBS = $(LIBRARY) $(LIBS)

brox_spatial_SOURCES_CPP = \
	brox_spatial_main.cpp \
	iio.cpp
brox_spatial_SOURCES_C =
brox_spatial_HEADERS = \
	iio.h
brox_spatial_OBJS =  $(brox_spatial_SOURCES_CPP:.cpp=.o) $(brox_spatial_SOURCES_C:.c=.o)
brox_spatial_LIBS = $(LIBRARY) $(LIBS)

brox_temporal_SOURCES_CPP = \
	brox_temporal_main.cpp \
	iio.cpp
brox_temporal_SOURCES_C =
brox_temporal_HEADERS = \
	iio.h
brox_temporal_OBJS =  $(brox_temporal_SOURCES_CPP:.cpp=.o) $(brox_temporal_SOURCES_C:.c=.o)
brox_temporal_LIBS = $(LIBRARY) $(LIBS)

tvl1flow_SOURCES_CPP = \
	tvl1flow_main.cpp \
	iio.cpp
tvl1flow_SOURCES_C =
tvl1flow_HEADERS = \
	iio.h
tvl1flow_OBJS =  $(tvl1flow_SOURCES_CPP:.cpp=.o) $(tvl1flow_SOURCES_C:.c=.o)
tvl1flow_LIBS = $(LIBRARY) $(LIBS)

tvl1occflow_SOURCES_CPP = \
	tvl1occflow_main.cpp \
	iio.cpp
tvl1occflow_SOURCES_C =
tvl1occflow_HEADERS = \
	iio.h
tvl1occflow_OBJS =  $(tvl1occflow_SOURCES_CPP:.cpp=.o) $(tvl1occflow_SOURCES_C:.c=.o)
tvl1occflow_LIBS = $(LIBRARY) $(LIBS)

SRCS_CPP = \
	$(of_SOURCES_CPP) \
	$(horn_schunck_classic_SOURCES_CPP) \
	$(horn_schunck_pyramidal_SOURCES_CPP) \
	$(brox_spatial_SOURCES_CPP) \
	$(brox_temporal_SOURCES_CPP) \
	$(tvl1flow_SOURCES_CPP)  \
	$(tvl1occflow_SOURCES_CPP) 

SRCS_C = \
	$(of_SOURCES_C) \
	$(horn_schunck_classic_SOURCES_C) \
	$(horn_schunck_pyramidal_SOURCES_C) \
	$(brox_spatial_SOURCES_C) \
	$(brox_temporal_SOURCES_C) \
	$(tvl1flow_SOURCES_C) \
	$(tvl1occflow_SOURCES_C)

HEADERS = \
	$(of_HEADERS) \
	$(horn_schunck_classic_HEADERS) \
	$(horn_schunck_pyramidal_HEADERS) \
	$(brox_spatial_HEADERS) \
	$(brox_temporal_HEADERS) \
	$(tvl1flow_HEADERS) \
	$(tvl1occflow_HEADERS)

horn_schunck_classic: $(horn_schunck_classic_OBJS) $(LIBRARY)
	$(LD) -o $@ $^ $(LDFLAGS) $(horn_schunck_classic_LIBS) $(LDADD)

horn_schunck_pyramidal: $(horn_schunck_pyramidal_OBJS) $(LIBRARY)
	$(LD) -o $@ $^ $(LDFLAGS) $(horn_schunck_pyramidal_LIBS) $(LDADD)

brox_spatial: $(brox_spatial_OBJS) $(LIBRARY)
	$(LD) -o $@ $^ $(LDFLAGS) $(brox_spatial_LIBS) $(LDADD)

brox_temporal: $(brox_temporal_OBJS) $(LIBRARY)
	$(LD) -o $@ $^ $(LDFLAGS) $(brox_temporal_LIBS) $(LDADD)

tvl1flow: $(tvl1flow_OBJS) $(LIBRARY)
	$(LD) -o $@ $^ $(LDFLAGS) $(tvl1flow_LIBS) $(LDADD)

tvl1occflow: $(tvl1occflow_OBJS) $(LIBRARY)
	$(LD) -o $@ $^ $(LDFLAGS) $(tvl1occflow_LIBS) $(LDADD)


$(LIBRARY): $(of_OBJS)
	$(AR) rvu $@ $^
	$(RANLIB) $@

.SUFFIXES: .c .o .cpp

## gcc-only version:
%.o : %.c
	$(COMPILE.c) -MD -o $@ $<
	@cp $*.d $*.P; \
	    sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	        -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
	    rm -f $*.d

%.o : %.cpp
	$(COMPILE.cpp) -MD -o $@ $<
	@cp $*.d $*.P; \
	    sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	        -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
	    rm -f $*.d

## general version:
# %.o : %.c
# 	@$(MAKEDEPEND); \
# 	    cp $*.d $*.P; \
# 	    sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
# 		-e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
# 	    rm -f $*.d
# 	$(COMPILE.c) -o $@ $<

# %.o : %.cpp
# 	@$(MAKEDEPEND); \
# 	    cp $*.d $*.P; \
# 	    sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
# 		-e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
# 	    rm -f $*.d
# 	$(COMPILE.cpp) -o $@ $<

.PHONY: clean distclean
clean:
	-rm -f $(PROGRAMS) $(LIBRARY) *.o  *~
	-rm -rf *.dSYM

distclean: clean
	-rm -f $(SRCS_CPP:.cpp=.P) $(SRCS_C:.c=.P)

count:
	 wc -l $(SRCS_CPP) $(HEADERS)

-include $(SRCS_CPP:.cpp=.P) $(SRCS_C:.c=.P)
