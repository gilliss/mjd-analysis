include $(MGDODIR)/buildTools/config.mk

# Give the list of applications, which must be the stems of cc files with 'main'.  There
# can be more than one.
APPS = 2nuBB_Systematics 2nuBB_RateByDet 2nuBB_RateByDS Submit AutoGenBatch
# The next two lines are important
SHLIB =
ARCHIVE =
TAMDIR ?= $(ROOTSYS)
# Include the correct flags,
INCLUDEFLAGS = $(CLHEP_INCLUDE_FLAGS) -I$(MGDODIR)/Base -I$(MGDODIR)/Root -I$(MGDODIR)/Transforms
INCLUDEFLAGS += -I$(MGDODIR)/Majorana -I$(MGDODIR)/MJDB $(ROOT_INCLUDE_FLAGS) -I$(TAMDIR)/inc -I$(TAMDIR)/include -I$(MGDODIR)/Tabree
INCLUDEFLAGS += -I$(GATDIR)/BaseClasses -I$(GATDIR)/MGTEventProcessing -I$(GATDIR)/MGOutputMCRunProcessing -I$(GATDIR)/Analysis -I$(GATDIR)/MJDAnalysis -I$(GATDIR)/DCProcs
LIBFLAGS = -L$(MGDODIR)/lib -lMGDORoot -lMGDOBase -lMGDOTransforms -lMGDOMajorana -lMGDOGerdaTransforms -lMGDOMJDB -lMGDOTabree
LIBFLAGS += -L$(GATDIR)/lib -lGATBaseClasses -lGATMGTEventProcessing -lGATMGOutputMCRunProcessing -lGATAnalysis -lGATMJDAnalysis -lGATDCProcs $(ROOT_LIB_FLAGS) -lSpectrum -lTreePlayer -L$(TAMDIR)/lib -lTAM

include $(MGDODIR)/buildTools/BasicMakefile
