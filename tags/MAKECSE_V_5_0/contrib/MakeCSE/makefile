# INSTALLATION INSTRUCTIONS
#
# (1) set HOMEDIR to the path to this makefile
# (2) select the system dependent options by
#     commenting or uncommenting the flags below
# (3) for the first build type: "make headers depend [platform]"
#     where platform is "gnu", "dec", "sgi", "sun", or "alaska"
# (4) for all subsequent builds type: "make [platform]"
#
#
# top level makefile - calls recursive make targets and links

# default target name
TARGET = MakeCSE
LIBRARY  = $(HOMEDIR)/lib/lib$(TARGET).a

# path to current directory
HOMEDIR = /u/jwfoulk/sawimme_fracture/MakeCSEv5

# directories:
OBJDIR = $(HOMEDIR)/obj
INCDIR = $(HOMEDIR)/inc
SRCDIR = $(HOMEDIR)/src

############################### build options #####################

##### SEACAS/ACCESS software
AX_CFLAGS = -I$(ACCESS)/inc -DACCESS
AX_LFLAGS = -L$(ACCESS)/lib -lexoIIv2c -lnetcdf

# grouped
OPT_CFLAGS = $(AX_CFLAGS) 
OPT_LFLAGS = $(AX_LFLAGS) 
####################################################################

# paths for file dependencies
DEPENDPATH = \
	-I$(INCDIR) \
	-I/usr/include \
	-I/usr/include/cxx \
	-I/usr/include/CC \
	-I$(ACCESS)/inc \
	-I/opt/SUNWspro/SC5.0/include/CC \
	-I/export/home/opt/SUNWspro/SC5.0/include/CC \
	-I/usr/lib/compaq/cxx-6.3.6.9/alpha-linux/include/cxx \
	-I/usr/include/g++ 

# Sun spends too much time searching standard directories
# DEPENDPATH = -Y -I$(INCDIR)

# object list from external file - defines OBJ
include object_list

# DEC
DEC_COMP_CC     = cxx
DEC_COMP_C      = cxx
DEC_AR          = ar
DEC_LINK        = $(DEC_COMP_CC)
DEC_TR          = $(OBJDIR)/cxx_repository
DEC_MY_CFLAGS   = -c -ptr $(DEC_TR) -msg_disable 161,111 -D__DEC__ -I$(INCDIR)
DEC_CFLAGS      = $(OPT_CFLAGS) $(DEC_MY_CFLAGS)
DEC_MORE_CFLAGS = -O3 -DNDEBUG 
#DEC_MORE_CFLAGS = -g 
DEC_LFLAGS      = $(OPT_LFLAGS) -lm -ptr $(DEC_TR)
DEC_ARFLAGS     = -ru
DEC_MAKEOPTS    =

# CA-CPLANT (DEC)
CPLANT_COMP_CC     = cxx
CPLANT_COMP_C      = cxx
CPLANT_AR          = ar
CPLANT_LINK        = $(CPLANT_COMP_CC)
CPLANT_TR          = $(OBJDIR)/cxx_repository
CPLANT_MY_CFLAGS   = -c -ptr $(DEC_TR) -I$(INCDIR) -I/usr/local/mpich/include -msg_display_number -msg_disable 161 -D__CPLANT__ -D__CA_CPLANT__
CPLANT_CFLAGS      = $(OPT_CFLAGS) $(CPLANT_MY_CFLAGS)
CPLANT_MORE_CFLAGS = -O3 -DNDEBUG 
#CPLANT_MORE_CFLAGS = -g 
CPLANT_LFLAGS      = $(OPT_LFLAGS) -lm -ptr $(DEC_TR) -L/usr/local/lam/lib
CPLANT_ARFLAGS     = -ru
CPLANT_MAKEOPTS    = -s

# alaska (DEC)
ALASKA_COMP_CC     = c++
ALASKA_COMP_C      = cc
ALASKA_AR          = ar
ALASKA_LINK        = $(ALASKA_COMP_CC)
ALASKA_TR          = $(DEC_TR)
ALASKA_MY_CFLAGS   = $(DEC_MY_CFLAGS) -I/usr/local/cplant/current/include -D__CPLANT__ -D__ALASKA__
ALASKA_CFLAGS      = $(OPT_CFLAGS) $(ALASKA_MY_CFLAGS)
ALASKA_MORE_CFLAGS = $(DEC_MORE_CFLAGS)
#ALASKA_MORE_CFLAGS = -g
ALASKA_LFLAGS      = $(DEC_LFLAGS) -L/usr/local/cplant/current/lib
ALASKA_ARFLAGS     = $(DEC_ARFLAGS)
ALASKA_MAKEOPTS    = $(DEC_MAKEOPTS)

# GNU
#GNU_COMP_CC     = /usr/community/g++/gcc-2.8.1/g++
GNU_COMP_CC     = g++
GNU_COMP_C      = gcc
GNU_AR          = ar
GNU_LINK        = $(GNU_COMP_CC)
GNU_CFLAGS      = -c -fguiding-decls -I$(INCDIR) $(OPT_CFLAGS)
GNU_MORE_CFLAGS = -O3 -DNDEBUG -Dpascal=""
GNU_LFLAGS      = $(OPT_LFLAGS)
GNU_ARFLAGS     = -ru
GNU_MAKEOPTS    = 

# SGI Irix 6.2 (or newer)
SGI_COMP_CC     = CC
SGI_COMP_C      = cc
SGI_AR          = ar
SGI_LINK        = $(SGI_COMP_CC)
SGI_MY_CFLAGS   = -c -mips4 -ptused -woff 1110,1155,1174,1188,1682 -LANG:exceptions=ON -D__SGI__ -I$(INCDIR)
SGI_CFLAGS      = $(OPT_CFLAGS) $(SGI_MY_CFLAGS)
SGI_MORE_CFLAGS = -O3 -OPT:Olimit=4000 -DNDEBUG 
#SGI_MORE_CFLAGS = -g
SGI_LFLAGS      = $(OPT_LFLAGS) -lm -w
SGI_ARFLAGS     = -ru
SGI_MAKEOPTS    = -sP
# parallel builds
# env. var. PARALLEL controls concurrency (default = 2)

# Sun
SUN_COMP_CC     = CC
SUN_COMP_C      = cc
SUN_AR          = ar
SUN_LINK        = $(SUN_COMP_CC)
#SUN_MY_CFLAGS   = -c -ptr$(HOMEDIR)/obj -I$(INCDIR) -features=rtti
SUN_MY_CFLAGS   = -c -instances=static -I$(INCDIR) -features=rtti
SUN_CFLAGS      = $(OPT_CFLAGS) $(SUN_MY_CFLAGS)
#SUN_MORE_CFLAGS = -g
SUN_MORE_CFLAGS = -fast -DNDEBUG 
#SUN_LFLAGS      = $(OPT_LFLAGS) -ptr$(HOMEDIR)/obj
SUN_LFLAGS      = $(OPT_LFLAGS)
SUN_ARFLAGS     = -ru
SUN_MAKEOPTS    = 

# utility targets
default:
	@ echo " no target specified:"
	@ echo "    [alaska] : build for Alaska CPLANT"
	@ echo "    [cplant] : build for CA-CPLANT"
	@ echo "       [dec] : build for DEC"
	@ echo "       [gnu] : build using GNU compiler"
	@ echo "       [sgi] : build for SGI"
	@ echo "       [sun] : build for Sun" 
	@ echo "     [clean] : remove object files and build temporaries"
	@ echo " [clean_all] : complete clean-up"
	@ echo "    [depend] : update dependencies"
	@ echo "   [headers] : copy all headers files to $(INCDIR)"

# create work folders/files
install_dirs:
	@ echo "creating working directories and files"
	@ mkdir $(INCDIR); \
	echo "DO NOT DELETE" > $(INCDIR)/inc.junk
	@ mkdir $(OBJDIR); \
	echo "DO NOT DELETE" > $(OBJDIR)/obj.junk
	@ mkdir $(OBJDIR)/cxx_repository; \
	echo "DO NOT DELETE" > $(OBJDIR)/cxx_repository/repo.junk

clean:
	- rm $(TARGET)
	- cd $(SRCDIR); make clean
	- cd $(OBJDIR); rm *.o
	- cd $(OBJDIR)/cxx_repository; rm -r *.o *.DB

clean_all: clean
	- cd $(INCDIR); rm *.h

depend:
	@ echo "cd $(SRCDIR)"; \
	cd $(SRCDIR); \
	make depend "DEPENDPATH = $(DEPENDPATH)"

headers:
	@ echo "cd $(SRCDIR)"; \
	cd $(SRCDIR); \
	make headers "INCDIR = $(INCDIR)"

touch:
	@ echo "cd $(SRCDIR)"; \
	cd $(SRCDIR); \
	make touch

# build targets
dec: dec_dir_check update_headers dec_objects
	@ echo " "; make $(TARGET) \
	"ARCH    = $(DEC_AR)" \
	"LINK    = $(DEC_LINK)" \
	"ARFLAGS = $(DEC_ARFLAGS)" \
	"LFLAGS  = $(DEC_LFLAGS)"

cplant: cplant_dir_check update_headers cplant_objects
	@ echo " "; make $(TARGET) \
	"ARCH    = $(CPLANT_AR)" \
	"LINK    = $(CPLANT_LINK)" \
	"ARFLAGS = $(CPLANT_ARFLAGS)" \
	"LFLAGS  = $(CPLANT_LFLAGS)"

alaska: alaska_dir_check update_headers alaska_objects
	@ echo " "; make $(TARGET) \
	"ARCH    = $(ALASKA_AR)" \
	"LINK    = $(ALASKA_LINK)" \
	"ARFLAGS = $(ALASKA_ARFLAGS)" \
	"LFLAGS  = $(ALASKA_LFLAGS)"

gnu: gnu_dir_check update_headers gnu_objects
	@ echo " "; make $(TARGET) \
	"ARCH    = $(GNU_AR)" \
	"LINK    = $(GNU_LINK)" \
	"ARFLAGS = $(GNU_ARFLAGS)" \
	"LFLAGS  = $(GNU_LFLAGS)"

sgi: sgi_dir_check sgi_update_headers sgi_objects
	@ echo " "; make $(TARGET) \
	"ARCH    = $(SGI_AR)" \
	"LINK    = $(SGI_LINK)" \
	"ARFLAGS = $(SGI_ARFLAGS)" \
	"LFLAGS  = $(SGI_LFLAGS)"

sun: sun_dir_check update_headers sun_objects
	@ echo " "; make $(TARGET) \
	"ARCH    = $(SUN_AR)" \
	"LINK    = $(SUN_LINK)" \
	"ARFLAGS = $(SUN_ARFLAGS)" \
	"LFLAGS  = $(SUN_LFLAGS)"

# dependent targets (DO NOT make DIRECTLY)
$(TARGET): $(LIBRARY)
	@ echo " linking..."	
	@ $(LINK) $(OBJDIR)/main.o -o $(TARGET) -L$(HOMEDIR)/lib -l$(TARGET) $(LFLAGS)
	@ echo "    ...wrote $(TARGET)"; echo " "

$(LIBRARY): $(OBJDIR)/obj.junk
	@ echo " building archive..."
	@ cd $(OBJDIR); $(ARCH) $(ARFLAGS) $(LIBRARY) $(OBJ)

objects:
	@ echo "cd $(SRCDIR)"; \
	cd $(SRCDIR); \
	make objects \
	"TARGET  = $(HOMEDIR)/$(TARGET) " \
	"INCDIR  = $(INCDIR)" \
	"OBJDIR  = $(OBJDIR)" \
	"TR      = $(TR)" \
	"COMP_CC = $(COMP_CC)" \
	"COMP_C  = $(COMP_C)" \
	"CFLAGS  = $(CFLAGS)" \
	"MORE_CFLAGS = $(MORE_CFLAGS)"

update_headers:
	@ cd $(SRCDIR); \
	make update_headers \
	"INCDIR = $(INCDIR)"
	@ touch $(INCDIR)/inc.junk

sgi_update_headers:
	@ cd $(SRCDIR); \
	make update_headers \
	"INCDIR   = $(INCDIR)" \
	"MAKEOPTS = $(SGI_MAKEOPTS)"
	@ touch $(INCDIR)/inc.junk

# GNU compilation instructions
gnu_objects:
	@ make objects \
	"TARGET      = $(TARGET)" \
	"COMP_CC     = $(GNU_COMP_CC)" \
	"COMP_C      = $(GNU_COMP_C)" \
	"CFLAGS      = $(GNU_CFLAGS)" \
	"MORE_CFLAGS = $(GNU_MORE_CFLAGS) " \
	"MAKEOPTS    = $(GNU_MAKEOPTS)"

sgi_objects:
	@ make $(SGI_MAKEOPTS) objects \
	"COMP_CC     = $(SGI_COMP_CC)" \
	"COMP_C      = $(SGI_COMP_C)" \
	"CFLAGS      = $(SGI_CFLAGS)" \
	"MORE_CFLAGS = $(SGI_MORE_CFLAGS) " \
	"MAKEOPTS    = $(SGI_MAKEOPTS)"

dec_objects:
	@ make objects \
	"TR          = $(TR)" \
	"COMP_CC     = $(DEC_COMP_CC)" \
	"COMP_C      = $(DEC_COMP_C)" \
	"CFLAGS      = $(DEC_CFLAGS)" \
	"MORE_CFLAGS = $(DEC_MORE_CFLAGS) " \
	"MAKEOPTS    = $(DEC_MAKEOPTS)"

cplant_objects:
	@ make objects \
	"TR          = $(TR)" \
	"COMP_CC     = $(CPLANT_COMP_CC)" \
	"COMP_C      = $(CPLANT_COMP_C)" \
	"CFLAGS      = $(CPLANT_CFLAGS)" \
	"MORE_CFLAGS = $(CPLANT_MORE_CFLAGS) " \
	"MAKEOPTS    = $(CPLANT_MAKEOPTS)"

alaska_objects:
	@ make objects \
	"TR          = $(TR)" \
	"COMP_CC     = $(ALASKA_COMP_CC)" \
	"COMP_C      = $(ALASKA_COMP_C)" \
	"CFLAGS      = $(ALASKA_CFLAGS)" \
	"MORE_CFLAGS = $(ALASKA_MORE_CFLAGS) " \
	"MAKEOPTS    = $(ALASKA_MAKEOPTS)"

sun_objects:
	@ make objects \
	"TR          = $(TR)" \
	"COMP_CC     = $(SUN_COMP_CC)" \
	"COMP_C      = $(SUN_COMP_C)" \
	"CFLAGS      = $(SUN_CFLAGS)" \
	"MORE_CFLAGS = $(SUN_MORE_CFLAGS) " \
	"MAKEOPTS    = $(SUN_MAKEOPTS)"

# make sure target directories exist (must be better way than to "cp" ?)
sgi_dir_check: dir_check

gnu_dir_check: dir_check

dec_dir_check: dir_check
	@ echo " template repository    : $(DEC_TR)"
	@ cp makefile $(DEC_TR); rm $(DEC_TR)/makefile

alaska_dir_check: dir_check
	@ echo " template repository    : $(ALASKA_TR)"
	@ cp makefile $(ALASKA_TR); rm $(ALASKA_TR)/makefile

cplant_dir_check: dir_check
	@ echo " template repository    : $(CPLANT_TR)"
	@ cp makefile $(CPLANT_TR); rm $(CPLANT_TR)/makefile

sun_dir_check: dir_check

dir_check: inc_dir_check obj_dir_check

inc_dir_check:
	@ echo " include file directory : $(INCDIR)"
	@ cp makefile $(INCDIR); rm $(INCDIR)/makefile

obj_dir_check:
	@ echo " object file directory  : $(OBJDIR)"
	@ cp makefile $(OBJDIR); rm $(OBJDIR)/makefile
