paraview_add_plugin(VDBWriter
  REQUIRED_ON_SERVER
  VERSION "1.0"
  SERVER_MANAGER_XML "${CMAKE_CURRENT_SOURCE_DIR}/../IOVDB/Resources/writers_pv_iovdb.xml"
  MODULE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/../IOVDB/vtk.module"
  MODULES VDBWriters)

#  MODULES VDBWriter::vtkVDBWriter
#  MODULE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/IOVDB/vtk.module")
