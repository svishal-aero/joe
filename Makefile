SHELL := /bin/bash

JOE_HOME := $(shell pwd)

include $(JOE_HOME)/Makefile.in

.PHONY: default help all clean_common clean_joe clean

default: all

all: libJoeCommon libJoe

libJoeCommon:
	$(MAKE) common -C $(JOE_HOME)/common

libJoe:
	$(MAKE) joe -C $(JOE_HOME)/joe

clean_common:
	$(MAKE) clean -i -C $(JOE_HOME)/common

clean_joe:
	$(MAKE) clean -i -C $(JOE_HOME)/joe

clean:
	$(MAKE) clean -i -C $(JOE_HOME)/common
	$(MAKE) clean -i -C ./joe
