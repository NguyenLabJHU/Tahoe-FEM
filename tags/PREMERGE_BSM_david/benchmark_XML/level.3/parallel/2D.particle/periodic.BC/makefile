# $Id: makefile,v 1.1 2005-05-25 00:54:53 paklein Exp $
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
	-@rm *.io*.ps* *.out *.geo *.run md.0.*rs1000* md.0.*rs500* console* *log *io.map *io.ID ./geometry/hex_lattice.n*.* ./geometry/four-atoms.n*.* ./geometry/five-atoms.n*.* *.echo.xml *.valid.xml
