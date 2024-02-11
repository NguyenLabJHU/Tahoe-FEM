# $Id: makefile,v 1.1 2004-09-02 20:44:05 rdorgan Exp $
#
# type "make clean" to remove test results
#

# list of subdirectories
SUB_DIR = \
c.0

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
