# $Id: makefile,v 1.1 2002-10-23 18:00:15 paklein Exp $
#
# NOTE: calling makefile must define:
#
#   COMP_CC - the c++ compiler
# CFLAGS_CC - c++ compiler flags
#   INC_DIR - path to the include directory
#    TB_INC - path to the include directory for the toolbox
#
PCH_FILE = cc1plus_ppc.pfe
PCH_INCL = \
	$(TB_INC)/Environment.h \
	$(TB_INC)/ios_fwd_decl.h \
	$(TB_INC)/ArrayT.h \
	$(TB_INC)/nArrayT.h \
	$(TB_INC)/iArrayT.h \
	$(TB_INC)/dArrayT.h \
	$(TB_INC)/StringT.h

$(PCH_FILE): pch.cpp $(PCH_INCL)
	@echo " rebuilding precompiled header file"
	$(COMP_CC) --dump-pch . $(CFLAGS_CC) -c pch.cpp

clean:
	rm $(PCH_FILE)