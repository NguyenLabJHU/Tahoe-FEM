# $Id: makefile,v 1.2 2004-07-12 16:25:54 paklein Exp $
#
# type "make clean" to remove test results
#

# list of sub-directories
# NOTE: if list is empty, set SUB_DIR = NONE
SUB_DIR = \
3D.particle

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
	-@rm *.io*.ps* *.out *.geo *.run *.log console* *io.ID *io.map explicit.dyn*.rs* *.echo.xml *.valid.xml
	-@cd geometry; rm *.n*.p*.geom* *.n*.part*
