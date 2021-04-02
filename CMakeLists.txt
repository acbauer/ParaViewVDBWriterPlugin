cmake_minimum_required(VERSION 3.10)
project(VDBWriter)

include(GNUInstallDirs)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
set(BUILD_SHARED_LIBS ON)

find_package(ParaView REQUIRED)


#paraview_add_plugin(VDBWriter
#  REQUIRED_ON_SERVER
#  VERSION "1.0"
#  MODULES VDBWriter::vtkVDBWriter
#  MODULE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/IOVDB/vtk.module")


# Scan the plugin file in order to set up internal data structures for building
# plugins.
paraview_plugin_scan(
  # The `paraview.plugin` file describing the plugin.
  PLUGIN_FILES "${CMAKE_CURRENT_SOURCE_DIR}/Plugin/paraview.plugin"
  # A result variable for the (enabled) plugins found during the scan.
  PROVIDES_PLUGINS plugins
  # Enable plugins during this scan by default.
  ENABLE_BY_DEFAULT ON)

# Build the plugins discovered during the scan.
paraview_plugin_build(
  PLUGINS ${plugins})
