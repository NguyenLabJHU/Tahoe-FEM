#
# type "make clean" to remove test results
#

# list of sub-directories
# NOTE: if list is empty, set SUB_DIR = NONE
SUB_DIR = \
cohesive.00 \
cohesive.01 \
cohesive.07

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
