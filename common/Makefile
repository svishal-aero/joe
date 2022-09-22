# Build makefile for common files
# Author: Vincent Terrapon
# Version: 1.0
# Date: 05/2009

SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .c .cpp .o

THIS_DIR := $(shell pwd)
JOE_HOME := $(THIS_DIR)/..

# Includes machine specific compilation options and dependencies of include files <file.h>
include $(JOE_HOME)/Makefile.in
include $(JOE_HOME)/common/Makefile.dep

COMMON_OBJ = $(JOE_HOME)/common/obj/MpiStuff.o     \
             $(JOE_HOME)/common/obj/Gp.o           \
             $(JOE_HOME)/common/obj/MiscUtils.o    \
             $(JOE_HOME)/common/obj/Ugp.o          \
             $(JOE_HOME)/common/obj/Param.o        \
             $(JOE_HOME)/common/obj/UgpWithTools.o \
             $(JOE_HOME)/common/obj/MshFilter.o    \
             $(JOE_HOME)/common/obj/UgpWithCv2.o   \
             $(JOE_HOME)/common/obj/tc_vec3d.o

.PHONY: default common help clean

default: common

help:
	@echo ' '
	@echo 'You must specifiy a target to make by typing: make <arg> '
	@echo 'where <arg> is one of the following options: ' 
	@echo 'all :    compiles common files used by joe, ray and charles'
	@echo 'clean :  clean common files'
	@echo ' '

common: $(JOE_HOME)/common/lib/libJoeCommon.so #$(COMMON_OBJ)

# Object files <file.o> dependencies
# Dependencies of the include files <file.h> are explicitly given in the above included files Makefile.dep
$(JOE_HOME)/common/obj/MpiStuff.o:         $(JOE_HOME)/common/src/MpiStuff.cpp $(MpiStuff_H)
	$(CC) $(CFLAGS) -c $< -o $@

$(JOE_HOME)/common/obj/Gp.o:               $(JOE_HOME)/common/src/Gp.cpp $(Gp_H)
	$(CC) $(CFLAGS) -c $< -o $@

$(JOE_HOME)/common/obj/MiscUtils.o:        $(JOE_HOME)/common/src/MiscUtils.cpp $(MiscUtils_H)
	$(CC) $(CFLAGS) -c $< -o $@

$(JOE_HOME)/common/obj/Ugp.o:              $(JOE_HOME)/common/src/Ugp.cpp $(Ugp_H)
	$(CC) $(CFLAGS) -c $< -o $@

$(JOE_HOME)/common/obj/Param.o:            $(JOE_HOME)/common/src/Param.cpp $(Param_H)
	$(CC) $(CFLAGS) -c $< -o $@

$(JOE_HOME)/common/obj/UgpWithTools.o:     $(JOE_HOME)/common/src/UgpWithTools.cpp $(UgpWithTools_H)
	$(CC) $(CFLAGS) -c $< -o $@

$(JOE_HOME)/common/obj/CdpFilter.o:        $(JOE_HOME)/common/src/CdpFilter.cpp $(CdpFilter_H)
	$(CC) $(CFLAGS) -c $< -o $@

$(JOE_HOME)/common/obj/MshFilter.o:        $(JOE_HOME)/common/src/MshFilter.cpp $(MshFilter_H)
	$(CC) $(CFLAGS) -c $< -o $@

$(JOE_HOME)/common/obj/UgpWithCv2.o:       $(JOE_HOME)/common/src/UgpWithCv2.cpp $(UgpWithCv2_H)
	$(CC) $(CFLAGS) -c $< -o $@

$(JOE_HOME)/common/obj/tc_vec3d.o:         $(JOE_HOME)/common/src/tc_vec3d.cpp $(tc_vec3d_H)
	$(CC) $(CFLAGS) -c $< -o $@

$(JOE_HOME)/common/obj/Logging.o:          $(JOE_HOME)/common/src/Logging.cpp $(CTR_TOMMIETOO_LOGGING_H)
	$(CC) $(CFLAGS) -c $< -o $@

$(JOE_HOME)/common/lib/libJoeCommon.so: $(JOE_HOME)/common/obj/MpiStuff.o     \
                                        $(JOE_HOME)/common/obj/Gp.o           \
                                        $(JOE_HOME)/common/obj/MiscUtils.o    \
                                        $(JOE_HOME)/common/obj/Ugp.o          \
                                        $(JOE_HOME)/common/obj/Param.o        \
                                        $(JOE_HOME)/common/obj/UgpWithTools.o \
                                        $(JOE_HOME)/common/obj/MshFilter.o    \
                                        $(JOE_HOME)/common/obj/UgpWithCv2.o   \
                                        $(JOE_HOME)/common/obj/tc_vec3d.o     \
                                        $(JOE_HOME)/common/obj/Logging.o
	$(CC) $(CFLAGS) --shared $^ -o $@ $(CLIBS)

# Cleaning
clean:
	rm -f ./obj/*.o
	rm -f ./lib/*.so
