AUTOMAKE_OPTIONS = foreign
INCLUDES= \
 -I m4 \
 -I`root-config --incdir`

libscalar_la_LDFLAGS =  `root-config --libs`
lib_LTLIBRARIES = \
	libscalar.la
libscalar_la_SOURCES = Scalar.cc Scalar.h
