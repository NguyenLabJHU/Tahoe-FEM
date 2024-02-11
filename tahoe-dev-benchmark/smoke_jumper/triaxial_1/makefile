# Make file to clean directory of Tahoe files

ifdef COMSPEC
  SUFFIX = .exe
  RM = del
else
  SUFFIX =
  RM = rm -f
endif

clean:
	$(RM) *.io0.*
	$(RM) *.io1.*
	$(RM) history.xml
	$(RM) runtime.properties
	$(RM) *echo*
	$(RM) *valid*
	$(RM) *.rs*
	$(RM) *.out