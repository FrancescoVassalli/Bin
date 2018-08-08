AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I m4
INCLUDES= \
	-I`root-config --incdir`
libfrootutils_la_LDFLAGS =\
	`root-config --libs`
lib_LTLIBRARIES = \
	libfrootutils.la
include_HEADERS = \
	Scalar.h\
	Photon.h

libfrootutils_la_SOURCES = \
	Scalar.cc\
	Photon_Dict.C\
	Scalar_Dict.C

BUILT_SOURCES = \
	testexternals.C

noinst_PROGRAMS = \
	testexternals

testexternals_LDADD = \
	libfrootutils.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

##############################################
# please add new classes in alphabetical order
# Rule for generating CINT dictionaries from class headers.
%_Dict.C: %.h %LinkDef.h
	rootcint -f $@ -c $(CINTFLAGS) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict* testexternals.C
