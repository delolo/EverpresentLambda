# MAKEFILE: A makefile is a text file that is referenced by the make command that 
# describes the building of targets, and contains information such as source-level 
# dependencies and build-order dependencies.
# http://help.eclipse.org/luna/index.jsp?topic=%2Forg.eclipse.cdt.doc.user%2Fconcepts%2Fcdt_c_makefile.htm


# Here are some Make variables to be used in the commands defined below.
# Make Macros do not need to be defined before their Macro Expansions,
# but they normally should be defined before they appear in any Rules.
# Consequently Make Macros often appear first in a Makefile.
INCRANDOM = -I include/tools
INCTOOLS  = -I include 

# Set default make recipe. This means that the command $ make will run $ make "main".
default: main

# Here are the compile recipes:
test: src/test.cpp
	g++ -Wall $(INCTOOLS) -o test.out src/test.cpp include/Tools.cpp
	@echo Successfully compiled to "test.out".

main: src/main.cpp
	g++ -Wall $(INCRANDOM) -o main.out src/main.cpp include/tools/random.cpp
	@echo Successfully compiled to "main.out".


# Here is the clean-up recipe. Typically it just deletes the binaries.
clean:
	rm -f *.out
	@echo Clean Done.