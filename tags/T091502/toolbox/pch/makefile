# $Id: makefile,v 1.1 2002-09-12 16:41:40 paklein Exp $
#
# NOTE: calling makefile must define:
#
#   COMP_CC - the c++ compiler
# CFLAGS_CC - c++ compiler flags
#   INC_DIR - path to the include directory
#
PCH_FILE = cc1plus_ppc.pfe
PCH_INCL = \
	$(INC_DIR)/Environment.h \
	$(INC_DIR)/ios_fwd_decl.h \
	$(INC_DIR)/iArrayT.h \
	$(INC_DIR)/dArrayT.h \
	$(INC_DIR)/StringT.h

$(PCH_FILE): pch.cpp $(PCH_INCL)
	@echo " rebuilding precompiled header file"
	$(COMP_CC) --dump-pch . $(CFLAGS_CC) -c pch.cpp

clean:
	rm $(PCH_FILE)