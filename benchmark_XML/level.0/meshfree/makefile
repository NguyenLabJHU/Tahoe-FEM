# $Id: makefile,v 1.4 2005-03-01 08:26:53 paklein Exp $
#
# type "make clean" to remove test results
#

# list of sub-directories
# NOTE: if list is empty, set SUB_DIR = NONE
SUB_DIR = NONE 

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
	-@rm *.io*.ps* *.out *.geo *.run *.echo.xml *.valid.xml RKPM.2D.expdyn.rs*of* RKPM.2D.contact.1.rs*of* SCNI.cell.2D.J2.1.rs*of*
