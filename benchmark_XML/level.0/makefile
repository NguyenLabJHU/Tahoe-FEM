# $Id: makefile,v 1.2 2004-03-06 17:24:28 paklein Exp $
#
# type "make clean" to remove test results
#

# list of subdirectories
SUB_DIR = \
	2D.elastostatic \
	diffusion

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
