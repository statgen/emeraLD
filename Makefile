SUBDIRS = src

PARENT_MAKE := Makefile.tool
include Makefile.inc

$(shell mkdir -p obj/tabix_util)
