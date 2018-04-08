VERSION=0.1

# Check for CXX11
CXX11_AVAIL ?= 1

CPP0X=$(shell if [ $(CCVERSION) \> 4.3.0 ] ; then echo " -std=c++0x" ; fi)

ifeq ($(CXX11_AVAIL), 0)
CPP0X=
endif

CURRENT_PATH := $(dir $(lastword $(MAKEFILE_LIST)))

LIB_PATH_GENERAL = ../libStatGen

ACTUAL_PATH = $(if $(filter-out /%,$(LIB_PATH_GENERAL)), \
        $(addprefix $(CURRENT_PATH), $(LIB_PATH_GENERAL)), \
        $(LIB_PATH_GENERAL))
INCLUDE_MAKE = $(ACTUAL_PATH)/Makefiles/$(PARENT_MAKE)

########################
# Include the base Makefile
include $(INCLUDE_MAKE)

# .phony: cloneLib

########################
$(INCLUDE_MAKE):
	@echo Unable to locate: $(INCLUDE_MAKE)
	@echo To change the location, set LIB_PATH_GENERAL to the appropriate path to libStatGen.  Or specify \"make LIB_PATH_GENERAL=yourPath\" or \"make LIB_PATH_GENERAL=yourPath\"
	@echo Use make cloneLib if you have git and want to clone the current libStatGen at that location.


cloneLib:
	@if test -d $(LIB_PATH_GENERAL); \
	then echo $(LIB_PATH_GENERAL) already exists; \
	else git clone git://github.com/statgen/libStatGen.git $(LIB_PATH_GENERAL); fi
	@echo Call make to compile libStatGen and this tool.
