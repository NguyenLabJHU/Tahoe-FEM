.SUFFIXES: .cpp .h .java .class .o

JAVA_OBJ = test.class app.class

CPP_OBJ = testClass.o testImpl.o

COMP_CC = c++
CC_FLAGS = -c -I /System/Library/Frameworks/JavaVM.framework/Headers

LD       = c++
LD_FLAGS = -dynamiclib -v
DSO_EXT  = dylib

all: $(JAVA_OBJ) java_native cpp_dso

java_native:
	javah test

# not needed by JNI on Mac OS X
#javah -stubs test

cpp_dso: $(CPP_OBJ)
	$(LD) $(LD_FLAGS) -o testClass.$(DSO_EXT) $(CPP_OBJ)

run:
	java app

clean:
	rm *.class *.$(DSO_EXT) *.o

.java.class:
	javac $<

.cpp.o:
	$(COMP_CC) $(CC_FLAGS) $<
