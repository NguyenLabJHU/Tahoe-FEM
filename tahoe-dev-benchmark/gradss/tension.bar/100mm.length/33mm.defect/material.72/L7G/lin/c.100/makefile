# $Id: makefile,v 1.1 2004-09-02 20:44:30 rdorgan Exp $
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
	-@ls -1 | grep -e io0 -e gp1\.ps | xargs rm;
	-@rm *.io*.ps* *.out *.geo *.run *.gp1.* *.echo.xml *.valid.xml *.case
