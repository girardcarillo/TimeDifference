# - Basic CMake setup
# Check version meets ou requirements
# Declare project, which will configure compiler for us
cmake_minimum_required(VERSION 3.3)
project(TimeDifference)

# Build a dynamic library from our sources
add_library(TimeDifference SHARED TimeDifference.h TimeDifference.cc)

# Link it to the FalaiseModule library
# This ensures the correct compiler flags, include paths
# and linker flags are applied to our dynamic library.
target_link_libraries(TimeDifference PUBLIC Falaise Falaise_ParticleIdentification)

# Install it:
install(TARGETS TimeDifference DESTINATION ${CMAKE_INSTALL_LIBDIR}/Falaise/modules)
