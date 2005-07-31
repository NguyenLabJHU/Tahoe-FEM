# $Id: makefile,v 1.1.1.1 2004-03-04 06:42:09 paklein Exp $
#
# type "make clean" to remove test results
#

# list of subdirectories
SUB_DIR = level.0

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
