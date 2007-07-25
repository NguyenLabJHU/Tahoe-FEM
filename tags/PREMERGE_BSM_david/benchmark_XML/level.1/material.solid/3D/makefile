# $Id: makefile,v 1.11 2005-01-21 18:13:25 paklein Exp $
#
# type "make clean" to remove test results
#

# list of sub-directories
# NOTE: if list is empty, set SUB_DIR = NONE
SUB_DIR = \
material.01 \
material.02 \
material.03 \
material.04 \
material.05 \
material.06 \
material.07 \
material.08 \
material.09 \
material.10 \
material.11 \
material.14 \
material.15 \
material.16 \
material.17 \
material.18 \
material.20 \
material.21 \
material.22 \
material.24 \
material.25 \
material.50 \
material.80 \
material.81 \
material.90 \
material.110

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
