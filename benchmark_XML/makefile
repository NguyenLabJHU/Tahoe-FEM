# $Id: makefile,v 1.5 2005-08-05 07:18:07 paklein Exp $
#
# type "make clean" to remove test results
#

# list of subdirectories
SUB_DIR = \
level.0 \
level.1 \
level.2

all:
	@echo "make clean to remove test results"
	
clean:
	@ for n in $(SUB_DIR); do \
		if test -d $$n; then \
			echo " cd $$n"; \
			cd $$n; \
			make clean; \
			cd ../; \
		elif test "$$n" != "NONE"; then \
			echo "*** WARNING: missing directory: $(CURR_DIR)/$$n"; \
			fi; \
	done    
	-@ rm tahoe.xsd tahoe.dtd

init: dummy

build: dummy

dummy:
	@echo "dummy target"
