##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I${G4_MAIN}/include \
	-isystem`root-config --incdir`

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

lib_LTLIBRARIES = \
   libg4eval_hp.la

libg4eval_hp_la_LIBADD = \
  -lCLHEP \
  -lffarawobjects\
  -lfun4all \
  -lg4detectors_io \
  -lintt_io \
  -lmbd_io \
  -lmicromegas_io \
  -lmvtx_io \
  -lphg4hit \
  -lphhepmc_io \
  -ltpc \
  -ltpc_io \
  -ltrackbase_historic_io \
  -ltrackbase_historic \
  -ltrack_io \
  -lgsl \
  -lgslcblas

pkginclude_HEADERS = \
  EventCounter_hp.h \
  MicromegasClusterEvaluator_hp.h \
  MicromegasEvaluator_hp.h \
  MicromegasGeometryContainer.h \
  MicromegasTrackEvaluator_hp.h \
  PHTrackPruner.h \
  SimEvaluator_hp.h \
  TrackingEvaluator_hp.h

ROOTDICTS = \
  MicromegasClusterEvaluator_hp_Dict.cc \
  MicromegasEvaluator_hp_Dict.cc \
  MicromegasGeometryContainer_Dict.cc \
  MicromegasTrackEvaluator_hp_Dict.cc \
  SimEvaluator_hp_Dict.cc \
  TrackingEvaluator_hp_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  MicromegasClusterEvaluator_hp_Dict_rdict.pcm \
  MicromegasEvaluator_hp_Dict_rdict.pcm \
  MicromegasGeometryContainer_Dict_rdict.pcm \
  MicromegasTrackEvaluator_hp_Dict_rdict.pcm \
  SimEvaluator_hp_Dict_rdict.pcm \
  TrackingEvaluator_hp_Dict_rdict.pcm

libg4eval_hp_la_SOURCES = \
  $(ROOTDICTS) \
  EventCounter_hp.cc \
  MicromegasClusterEvaluator_hp.cc \
  MicromegasEvaluator_hp.cc \
  MicromegasGeometryContainer.cc \
  MicromegasTrackEvaluator_hp.cc \
  PHTrackPruner.cc \
  SimEvaluator_hp.cc \
  TrackingEvaluator_hp.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_g4eval_hp

testexternals_g4eval_hp_SOURCES = testexternals.cc
testexternals_g4eval_hp_LDADD = libg4eval_hp.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
