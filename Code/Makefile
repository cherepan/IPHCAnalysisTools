include make.inc
#################################################
### Defining Compiling variables              ###
#################################################
#-Wno-parentheses
CXX           = c++ -Wno-deprecated -Wall -std=c++0x
#-m32  -Wno-deprecated
LD            = g++ -Wno-deprecated -Wall
#-m32  -Wno-deprecated 
INSTALL	      = /usr/bin/install
SHELL = /bin/sh

####################################
###                              ###
####################################

SHAREDLIBFLAGS =-L./CommonUtils/lib -lTauNtuple -L./DataFormats/lib -lDataFormats
SHAREDCXXFLAGS =-I./CommonUtils -I./DataFormats
SUBDIRS = CommonUtils/ 
TAUSPINERFLAG = 

override CXXFLAGS += 
ALLCXXFLAGS = $(CXXFLAGS) -I$(ROOTSYS)/include -c -I/include $(SHAREDCXXFLAGS) -c -I./ -c
LDFLAGS       =  
ROOTLIBS      =  -L$(ROOTSYS)/lib -L/usr/lib/ -L/lib/i686/  -lCore -lCint -lHist -lGraf  -lGraf3d -lGpad -lTree -lRint -lReflexDict -lReflex -lPostscript -lMatrix -lPhysics -lGui -LObj -lMinuit2 -lThread -lTreePlayer -lRIO -lMathCore -rdynamic -Wl,--rpath $(ROOTSYS)/lib

# -ldl -lpthread

DEFS  	      = -DSTANDALONE=1
LIBS          = $(ROOTLIBS) $(SHAREDLIBFLAGS)

TARGETS       =	Ntuple_Controller \
		NtupleReader \
		Parameters \
		Selection_Base \
		Selection_Factory \
		Selection \
		Tables \
		Plots \
		DoubleEventRemoval \
		PDG_Var \
		Tools \
		SkimConfig \
		HistoConfig \
		Example \
		TauSolver \
		TauSpinExample \
		TauSpinerInterface \
		ReferenceScaleFactors \
		rochcor2012jan22 \
		Objects \
		PUReweight \
		ScaleFactor \
		tauTrigSFreader \
		DataMCCorrections \
		UncertaintyValue

CINTTARGETS = 

# now include conditional options set by flags in make.inc
ifdef USE_TauSpinner
	SUBDIRS += TauSpiner/
	SHAREDLIBFLAGS +=	-L./TauSpiner/lhapdf/lib -lLHAPDF \
				-L./TauSpiner/tauola++/1.1.3/lib -lTauolaCxxInterface -lTauolaFortran -lTauolaTauSpinner \
				-L./TauSpiner/hepmc/install/lib -lHepMC
        SHAREDCXXFLAGS +=	-I./TauSpiner/hepmc/install/include \
				-I./TauSpiner/tauola++/1.1.3/include/Tauola/ \
				-I./TauSpiner/tauola++/1.1.3/include/TauSpinner/ \
				-I./TauSpiner/tauola++/1.1.3/include/ \
				-I./TauSpiner/lhapdf/include/ \
				-I./TauSpiner/pythia8/176/include/
	DEFS += -DUSE_TauSpinner=1
endif

ifdef USE_SVfit
	TARGETS += SVfitProvider DataStorage SVFitStorage
	CINTTARGETS += SVFitObject
	SHAREDLIBFLAGS += -L./CommonUtils/lib -lSVfit
	DEFS += -DUSE_SVfit=1
endif

ifdef USE_cherepanov
	TARGETS += MyTest ZTauHTauH SingleMuSkim ZTauMuTauH TTBar
	DEFS += -DUSE_cherepanov=1
endif


ifdef USE_goe
	TARGETS += MyTest
	DEFS += -DUSE_goe=1
endif

ifdef USE_lebihan
	TARGETS += 
	DEFS += -DUSE_lebihan=1
endif

ifdef USE_gbourgat
	TARGETS += 
	DEFS += -DUSE_lebihan=1
endif



HDRS          = $(TARGETS:=.h) 		
SRCS          = $(TARGETS:=.cxx)  Analysis.cxx
OBJS          = $(SRCS:.cxx=.o)

PROGRAM       = Analysis.exe

$(PROGRAM):     $(OBJS)
	@echo "Linking $(PROGRAM) ..."
	@echo $(LD) $(LDFLAGS) i386_linux/*.o $(LIBS) -o $(PROGRAM)
	@$(LD) $(LDFLAGS) i386_linux/*.o $(LIBS) -o $(PROGRAM)
	@echo "done"

VPATH = utilities:i386_linux

ifdef USE_cherepanov
	vpath %.cxx cherepanov
	vpath %.h cherepanov
endif
ifdef USE_goe
	vpath %.cxx goe
	vpath %.h goe
endif
ifdef USE_lebihan
	vpath %.cxx lebihan
	vpath %.h lebihan
endif

ifdef USE_gbourgat
	vpath %.cxx gbourgat
	vpath %.h gbourgat
endif

vpath %.cxx CommonFiles
vpath %.h CommonFiles

#	@$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o ../$(PROGRAM)

##########################################################
###                                                    ###
##########################################################

$(OBJS): %.o : %.cxx
	$(CXX) $(ALLCXXFLAGS) $(DEFS) $< -o i386_linux/$@ 

.PHONY: clean cleanall cleandf all dataformats install sharedlib 

install: dataformats Analysis.exe


dataformats: 
	@cd DataFormats; gmake install; cd ../

clean:
	@rm i386_linux/*.o
	@rm Analysis.exe

cleandf:
	@cd DataFormats; gmake clean; cd ../

cleanall:
	@for subdir in $(SUBDIRS); do \
		echo "Making all in $$subdir"; \
		cd $$subdir && $(MAKE) -f Makefile.standalone clean && cd ../../; \
	done
	@cd DataFormats; gmake clean; cd ../
	@rm i386_linux/*.o
	@rm Analysis.exe

all: sharedlib dataformats install

sharedlib:
	@for subdir in $(SUBDIRS); do \
		echo "Making all in $$subdir"; \
		if [ $$subdir == "TauSpiner/" ]; then \
			echo with flag; \
			cd $$subdir && $(MAKE) -f Makefile.standalone $(TAUSPINERFLAG) && cd ../ && pwd; \
		else \
			echo without flag ; \
			cd $$subdir echo $(CXXFLAGS) && $(MAKE) -f Makefile.standalone CXXFLAGS+='$(CXXFLAGS)' ; \
			pwd; \
			if [ $$subdir == "CommonUtils/" ]; then \
				echo $(USE_SVfit); \
				if [ $(USE_SVfit) -eq 1 ]; then \
					$(MAKE) -f Makefile.SVfitStandalone CXXFLAGS+='$(CXXFLAGS)' ; \
				fi \
			fi \
			&& cd ../ && pwd; \
		fi \
	done
	$(install)



