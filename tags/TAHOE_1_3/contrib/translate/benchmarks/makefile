#
# type "make clean" to remove test results
#

# executable
EXEC = $(TAHOE)/contrib/translate/translate
MAKE = make
ECHO = echo

# if you have Exodus libraries
EXODUS = exodus
# else use
#EXODUS = NONE



# list of subdirectories
SUB_DIR = abaqus ensight $(EXODUS) patran tahoe

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
	@ for n in $(SUB_DIR); do \
		if test -d $$n; then \
			$(ECHO) " cd $$n"; \
			cd $$n; \
			$(MAKE) trans "EXEC=$(EXEC)" "MAKE=$(MAKE)" "ECHO=$(ECHO)"; \
			cd ../; \
		elif test "$$n" != "NONE"; then \
			$(ECHO) "*** WARNING: missing directory: $(CURR_DIR)/$$n"; \
		fi; \
	done
