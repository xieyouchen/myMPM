# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.29

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.29.2/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.29.2/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/xieyouchen/usually/Projects/LiuLong/Project/myMPM

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/xieyouchen/usually/Projects/LiuLong/Project/myMPM/build

# Include any dependencies generated for this target.
include CMakeFiles/myMPM.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/myMPM.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/myMPM.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/myMPM.dir/flags.make

CMakeFiles/myMPM.dir/main.cpp.o: CMakeFiles/myMPM.dir/flags.make
CMakeFiles/myMPM.dir/main.cpp.o: /Users/xieyouchen/usually/Projects/LiuLong/Project/myMPM/main.cpp
CMakeFiles/myMPM.dir/main.cpp.o: CMakeFiles/myMPM.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/xieyouchen/usually/Projects/LiuLong/Project/myMPM/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/myMPM.dir/main.cpp.o"
	/Applications/Xcode-14.2.0.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/myMPM.dir/main.cpp.o -MF CMakeFiles/myMPM.dir/main.cpp.o.d -o CMakeFiles/myMPM.dir/main.cpp.o -c /Users/xieyouchen/usually/Projects/LiuLong/Project/myMPM/main.cpp

CMakeFiles/myMPM.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/myMPM.dir/main.cpp.i"
	/Applications/Xcode-14.2.0.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/xieyouchen/usually/Projects/LiuLong/Project/myMPM/main.cpp > CMakeFiles/myMPM.dir/main.cpp.i

CMakeFiles/myMPM.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/myMPM.dir/main.cpp.s"
	/Applications/Xcode-14.2.0.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/xieyouchen/usually/Projects/LiuLong/Project/myMPM/main.cpp -o CMakeFiles/myMPM.dir/main.cpp.s

# Object files for target myMPM
myMPM_OBJECTS = \
"CMakeFiles/myMPM.dir/main.cpp.o"

# External object files for target myMPM
myMPM_EXTERNAL_OBJECTS =

myMPM: CMakeFiles/myMPM.dir/main.cpp.o
myMPM: CMakeFiles/myMPM.dir/build.make
myMPM: CMakeFiles/myMPM.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/xieyouchen/usually/Projects/LiuLong/Project/myMPM/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable myMPM"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/myMPM.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/myMPM.dir/build: myMPM
.PHONY : CMakeFiles/myMPM.dir/build

CMakeFiles/myMPM.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/myMPM.dir/cmake_clean.cmake
.PHONY : CMakeFiles/myMPM.dir/clean

CMakeFiles/myMPM.dir/depend:
	cd /Users/xieyouchen/usually/Projects/LiuLong/Project/myMPM/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/xieyouchen/usually/Projects/LiuLong/Project/myMPM /Users/xieyouchen/usually/Projects/LiuLong/Project/myMPM /Users/xieyouchen/usually/Projects/LiuLong/Project/myMPM/build /Users/xieyouchen/usually/Projects/LiuLong/Project/myMPM/build /Users/xieyouchen/usually/Projects/LiuLong/Project/myMPM/build/CMakeFiles/myMPM.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/myMPM.dir/depend

