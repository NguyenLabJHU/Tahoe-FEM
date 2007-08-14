# $Id: makefile,v 1.4 2007-08-14 16:16:43 d-farrell2 Exp $
#
# type "make clean" to remove test results
#

# list of sub-directories
# NOTE: if list is empty, set SUB_DIR = NONE
SUB_DIR = \
interpolate \
interpolate.meshfree.1 \
interpolate.meshfree.2 \
projection \
meshfree \
meshfree.1 \
meshfree.2 \
meshfree.3 \
BSM

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
	-@rm *.io*.ps* *.out *.geo *.run *.log

