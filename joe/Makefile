# Build makefile for joe
# Author: Vincent Terrapon
# Version: 1.0
# Date: 05/2009

SHELL := /bin/bash

.SUFFIXES:
.SUFFIXES: .c .cpp .o

THIS_DIR := $(shell pwd)
JOE_HOME := $(THIS_DIR)/..

# Includes machine specific compilation options and dependencies of include files <file.h>
include $(JOE_HOME)/Makefile.in
include $(JOE_HOME)/common/Makefile.dep
include $(JOE_HOME)/joe/Makefile.dep

JOE_OBJS = $(JOE_HOME)/joe/obj/UgpWithCvCompFlow.o   \
           $(JOE_HOME)/joe/obj/Scalars.o             \
					 $(JOE_HOME)/joe/obj/JoeWithModels.o       \
					 $(JOE_HOME)/joe/obj/UgpWithCvCompFlowAD.o \
					 $(JOE_HOME)/joe/obj/ScalarsAD.o           \
					 $(JOE_HOME)/joe/obj/JoeWithModelsAD.o

#JOE_LIBS = $(CLIBS) -lm

#OBJS = $(OBJS_AUX)

.PHONY: default joe clean

default: joe

joe: $(JOE_HOME)/joe/lib/libJoe.so

$(JOE_HOME)/joe/lib/libJoe.so: $(JOE_OBJS)
	$(CC) $(CFLAGS) --shared $^ -o $@ $(CLIBS) -L../common/lib -lJoeCommon

$(JOE_HOME)/joe/obj/UgpWithCvCompFlow.o: $(JOE_HOME)/joe/src/UgpWithCvCompFlow.cpp $(UgpWithCvCompFlow_H)
	$(CC) $(CFLAGS) -I$(JOE_HOME)/common/include -c $< -o $@

$(JOE_HOME)/joe/obj/UgpWithCvCompFlowAD.o: $(JOE_HOME)/joe/src/UgpWithCvCompFlowAD.cpp $(UgpWithCvCompFlowAD_H)
	$(CC) $(CFLAGS) -I$(JOE_HOME)/common/include -c $< -o $@

$(JOE_HOME)/joe/obj/Scalars.o: $(JOE_HOME)/joe/src/Scalars.cpp $(UgpWithCvCompFlow_H)
	$(CC) $(CFLAGS) -I$(JOE_HOME)/common/include -c $< -o $@

$(JOE_HOME)/joe/obj/ScalarsAD.o: $(JOE_HOME)/joe/src/ScalarsAD.cpp $(UgpWithCvCompFlowAD_H)
	$(CC) $(CFLAGS) -I$(JOE_HOME)/common/include -c $< -o $@

$(JOE_HOME)/joe/obj/JoeWithModels.o: $(JOE_HOME)/joe/src/JoeWithModels.cpp $(JoeWithModels_H)
	$(CC) $(CFLAGS) -I$(JOE_HOME)/common/include -c $< -o $@

$(JOE_HOME)/joe/obj/JoeWithModelsAD.o: $(JOE_HOME)/joe/src/JoeWithModelsAD.cpp $(JoeWithModelsAD_H)
	$(CC) $(CFLAGS) -I$(JOE_HOME)/common/include -c $< -o $@

clean:
	rm -f $(JOE_HOME)/joe/obj/*.o
	rm -f $(JOE_HOME)/joe/lib/*.so
