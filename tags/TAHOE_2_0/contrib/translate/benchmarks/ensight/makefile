#
# type "make clean" to remove test results
#

# list of subdirectories
SUB_DIR = NONE

all:
	@$(ECHO) "make clean to remove test results"
	@$(ECHO) "make trans to run benchmarks"

clean:
	@ for n in $(SUB_DIR); do \
		if test -d $$n; then \
			$(ECHO) " cd $$n"; \
			cd $$n; \
			$(MAKE) clean "MAKE=$(MAKE)" "ECHO=$(ECHO)"; \
			cd ../; \
		elif test "$$n" != "NONE"; then \
			$(ECHO) "*** WARNING: missing directory: $(CURR_DIR)/$$n"; \
		fi; \
	done    
	-@rm trans* 

trans:
	-@$(EXEC) run.batch > trans.log
