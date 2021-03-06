# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

# Default target executed when no arguments are given to make.
default_target: all
.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/benjc/Documents/treenome

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/benjc/Documents/treenome

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake --regenerate-during-build -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache
.PHONY : rebuild_cache/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/usr/bin/ccmake -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/benjc/Documents/treenome/CMakeFiles /home/benjc/Documents/treenome//CMakeFiles/progress.marks
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/benjc/Documents/treenome/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean
.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named treenome

# Build rule for target.
treenome: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 treenome
.PHONY : treenome

# fast build rule for target.
treenome/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/treenome.dir/build.make CMakeFiles/treenome.dir/build
.PHONY : treenome/fast

src/GTree.o: src/GTree.cpp.o
.PHONY : src/GTree.o

# target to build an object file
src/GTree.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/treenome.dir/build.make CMakeFiles/treenome.dir/src/GTree.cpp.o
.PHONY : src/GTree.cpp.o

src/GTree.i: src/GTree.cpp.i
.PHONY : src/GTree.i

# target to preprocess a source file
src/GTree.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/treenome.dir/build.make CMakeFiles/treenome.dir/src/GTree.cpp.i
.PHONY : src/GTree.cpp.i

src/GTree.s: src/GTree.cpp.s
.PHONY : src/GTree.s

# target to generate assembly for a file
src/GTree.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/treenome.dir/build.make CMakeFiles/treenome.dir/src/GTree.cpp.s
.PHONY : src/GTree.cpp.s

src/InputFile.o: src/InputFile.cpp.o
.PHONY : src/InputFile.o

# target to build an object file
src/InputFile.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/treenome.dir/build.make CMakeFiles/treenome.dir/src/InputFile.cpp.o
.PHONY : src/InputFile.cpp.o

src/InputFile.i: src/InputFile.cpp.i
.PHONY : src/InputFile.i

# target to preprocess a source file
src/InputFile.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/treenome.dir/build.make CMakeFiles/treenome.dir/src/InputFile.cpp.i
.PHONY : src/InputFile.cpp.i

src/InputFile.s: src/InputFile.cpp.s
.PHONY : src/InputFile.s

# target to generate assembly for a file
src/InputFile.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/treenome.dir/build.make CMakeFiles/treenome.dir/src/InputFile.cpp.s
.PHONY : src/InputFile.cpp.s

src/Node.o: src/Node.cpp.o
.PHONY : src/Node.o

# target to build an object file
src/Node.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/treenome.dir/build.make CMakeFiles/treenome.dir/src/Node.cpp.o
.PHONY : src/Node.cpp.o

src/Node.i: src/Node.cpp.i
.PHONY : src/Node.i

# target to preprocess a source file
src/Node.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/treenome.dir/build.make CMakeFiles/treenome.dir/src/Node.cpp.i
.PHONY : src/Node.cpp.i

src/Node.s: src/Node.cpp.s
.PHONY : src/Node.s

# target to generate assembly for a file
src/Node.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/treenome.dir/build.make CMakeFiles/treenome.dir/src/Node.cpp.s
.PHONY : src/Node.cpp.s

src/SeqRead.o: src/SeqRead.cpp.o
.PHONY : src/SeqRead.o

# target to build an object file
src/SeqRead.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/treenome.dir/build.make CMakeFiles/treenome.dir/src/SeqRead.cpp.o
.PHONY : src/SeqRead.cpp.o

src/SeqRead.i: src/SeqRead.cpp.i
.PHONY : src/SeqRead.i

# target to preprocess a source file
src/SeqRead.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/treenome.dir/build.make CMakeFiles/treenome.dir/src/SeqRead.cpp.i
.PHONY : src/SeqRead.cpp.i

src/SeqRead.s: src/SeqRead.cpp.s
.PHONY : src/SeqRead.s

# target to generate assembly for a file
src/SeqRead.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/treenome.dir/build.make CMakeFiles/treenome.dir/src/SeqRead.cpp.s
.PHONY : src/SeqRead.cpp.s

src/TreeTop.o: src/TreeTop.cpp.o
.PHONY : src/TreeTop.o

# target to build an object file
src/TreeTop.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/treenome.dir/build.make CMakeFiles/treenome.dir/src/TreeTop.cpp.o
.PHONY : src/TreeTop.cpp.o

src/TreeTop.i: src/TreeTop.cpp.i
.PHONY : src/TreeTop.i

# target to preprocess a source file
src/TreeTop.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/treenome.dir/build.make CMakeFiles/treenome.dir/src/TreeTop.cpp.i
.PHONY : src/TreeTop.cpp.i

src/TreeTop.s: src/TreeTop.cpp.s
.PHONY : src/TreeTop.s

# target to generate assembly for a file
src/TreeTop.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/treenome.dir/build.make CMakeFiles/treenome.dir/src/TreeTop.cpp.s
.PHONY : src/TreeTop.cpp.s

src/cli.o: src/cli.cpp.o
.PHONY : src/cli.o

# target to build an object file
src/cli.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/treenome.dir/build.make CMakeFiles/treenome.dir/src/cli.cpp.o
.PHONY : src/cli.cpp.o

src/cli.i: src/cli.cpp.i
.PHONY : src/cli.i

# target to preprocess a source file
src/cli.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/treenome.dir/build.make CMakeFiles/treenome.dir/src/cli.cpp.i
.PHONY : src/cli.cpp.i

src/cli.s: src/cli.cpp.s
.PHONY : src/cli.s

# target to generate assembly for a file
src/cli.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/treenome.dir/build.make CMakeFiles/treenome.dir/src/cli.cpp.s
.PHONY : src/cli.cpp.s

src/main.o: src/main.cpp.o
.PHONY : src/main.o

# target to build an object file
src/main.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/treenome.dir/build.make CMakeFiles/treenome.dir/src/main.cpp.o
.PHONY : src/main.cpp.o

src/main.i: src/main.cpp.i
.PHONY : src/main.i

# target to preprocess a source file
src/main.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/treenome.dir/build.make CMakeFiles/treenome.dir/src/main.cpp.i
.PHONY : src/main.cpp.i

src/main.s: src/main.cpp.s
.PHONY : src/main.s

# target to generate assembly for a file
src/main.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/treenome.dir/build.make CMakeFiles/treenome.dir/src/main.cpp.s
.PHONY : src/main.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... rebuild_cache"
	@echo "... treenome"
	@echo "... src/GTree.o"
	@echo "... src/GTree.i"
	@echo "... src/GTree.s"
	@echo "... src/InputFile.o"
	@echo "... src/InputFile.i"
	@echo "... src/InputFile.s"
	@echo "... src/Node.o"
	@echo "... src/Node.i"
	@echo "... src/Node.s"
	@echo "... src/SeqRead.o"
	@echo "... src/SeqRead.i"
	@echo "... src/SeqRead.s"
	@echo "... src/TreeTop.o"
	@echo "... src/TreeTop.i"
	@echo "... src/TreeTop.s"
	@echo "... src/cli.o"
	@echo "... src/cli.i"
	@echo "... src/cli.s"
	@echo "... src/main.o"
	@echo "... src/main.i"
	@echo "... src/main.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

