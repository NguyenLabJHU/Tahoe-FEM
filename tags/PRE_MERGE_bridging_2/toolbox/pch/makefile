# $Id: makefile,v 1.3 2002-11-25 07:11:43 paklein Exp $
#
# NOTE: calling makefile must define:
#
#     COMP_CC - the c++ compiler
#   CFLAGS_CC - c++ compiler flags
#     INC_DIR - path to the include directory
# MAKE_DEPEND - makedepend (usu. in X11)
#

.SUFFIXES: .pfe .cpp .d

PCH_FILE = cc1plus_ppc.pfe

all: $(PCH_FILE)

clean:
	rm $(PCH_FILE) *.d *.d.bak

force:

.cpp.d:
	@echo " rebuilding precompiled header dependencies"; \
	echo "# DO NOT DELETE" > $@; \
	$(MAKE_DEPEND) -o .pfe -f$@ $<

.cpp.pfe:
	@echo " rebuilding precompiled header file"; \
	$(COMP_CC) --dump-pch . $(CFLAGS_CC) -c $<; \
	echo " pch: `du -sk $(PCH_FILE)`"

# dependencies
include cc1plus_ppc.d
