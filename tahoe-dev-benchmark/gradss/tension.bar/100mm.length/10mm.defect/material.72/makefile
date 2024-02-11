# $Id: makefile,v 1.1 2004-09-02 20:42:13 rdorgan Exp $
#
# type "make clean" to remove test results
#

# list of subdirectories
SUB_DIR = \
H16L \
L4L \
L6G \
L7G \
R8L

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
