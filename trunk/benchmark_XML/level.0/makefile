# $Id: makefile,v 1.17 2005-11-06 00:38:41 paklein Exp $
#
# type "make clean" to remove test results
#

# list of subdirectories
SUB_DIR = \
matrix_check \
2D.spring \
3D.spring \
2D.elastodynamic \
2D.elastostatic \
3D.elastostatic \
3D.elastodynamic \
diffusion \
coupled \
enhanced.strain \
integrator \
2D.particle \
3D.particle \
adhesion \
axisymmetric \
meshfree \
inputoutput \
bridging

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
