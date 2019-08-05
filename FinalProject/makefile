# Start of the makefile
# Defining variables
cpp_comp = mpiicpc 

# Makefile
all: structure

structure: main.cpp
	mpiicpc -std=c++14 -I ./eigen-eigen-323c052e1731  main.cpp -o main_exe

.PHONY: clean

# Cleaning everything
clean:
	rm *_exe
# End of the makefile
