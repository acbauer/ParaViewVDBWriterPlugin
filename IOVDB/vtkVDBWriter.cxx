/*=========================================================================

  Program:   ParaView
  Module:    vtkVDBWriter.cxx

  Copyright (c) Kitware, Inc.
  All rights reserved.
  See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkVDBWriter.h"

#include "vtkCellData.h"
#include "vtkCommunicator.h"
#include "vtkDiscretizableColorTransferFunction.h"
#include "vtkFloatArray.h"
#include "vtkImageData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMath.h"
#include "vtkMultiProcessController.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPointSet.h"
#include "vtkScalarsToColors.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkUnsignedCharArray.h"

#include "vtksys/FStream.hxx"
#include <vtksys/SystemTools.hxx>

#include <algorithm>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include <openvdb/openvdb.h>
#include <openvdb/points/PointConversion.h>
#include <openvdb/points/PointCount.h>
#include <openvdb/tools/Dense.h>

namespace
{
std::string GetVDBGridName(const char* arrayName, int component, int numberOfComponents)
{
  std::string vdbName = arrayName;
  if (numberOfComponents != 1 and numberOfComponents != 3)
  {
    vdbName = vdbName + "_" + std::to_string(component);
  }
  return vdbName;
}

}

vtkStandardNewMacro(vtkVDBWriter);
vtkCxxSetObjectMacro(vtkVDBWriter, Controller, vtkMultiProcessController);
vtkCxxSetObjectMacro(vtkVDBWriter, LookupTable, vtkScalarsToColors);
//-----------------------------------------------------------------------------
vtkVDBWriter::vtkVDBWriter()
{
  // openvdb::initialize() can be called multiple times
  openvdb::initialize();
  this->FileName = nullptr;
  this->WriteAllTimeSteps = false;
  this->Controller = nullptr;
  this->CurrentTimeIndex = 0;
  this->NumberOfTimeSteps = 1;
  this->SetController(vtkMultiProcessController::GetGlobalController());


  //this->ArrayName = nullptr;
  this->Component = 0;
  this->LookupTable = nullptr;
  this->EnableColoring = false;
  this->EnableAlpha = false;
  //this->SetArrayName("vtkVDBWriterColors");
}

//-----------------------------------------------------------------------------
vtkVDBWriter::~vtkVDBWriter()
{
  this->SetFileName(nullptr);
  this->SetController(nullptr);
  //this->SetArrayName(nullptr);
}

//-----------------------------------------------------------------------------
int vtkVDBWriter::FillInputPortInformation(int vtkNotUsed(port), vtkInformation* info)
{
  info->Append(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
  info->Append(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  return 1;
}

//----------------------------------------------------------------------------
int vtkVDBWriter::ProcessRequest(
  vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  if (request->Has(vtkStreamingDemandDrivenPipeline::REQUEST_UPDATE_EXTENT()))
  {
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(),
       (this->Controller ? this->Controller->GetNumberOfProcesses() : 1));
     inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(),
       (this->Controller ? this->Controller->GetLocalProcessId() : 0));
     inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(), 0);

     if (this->WriteAllTimeSteps)
     {
       double* inTimes =
         inputVector[0]->GetInformationObject(0)->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
       if (inTimes && this->WriteAllTimeSteps)
       {
         double timeReq = inTimes[this->CurrentTimeIndex];
         inputVector[0]->GetInformationObject(0)->Set(
           vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), timeReq);
       }
     }
     return 1;
  }
  else if (request->Has(vtkStreamingDemandDrivenPipeline::REQUEST_INFORMATION()))
  {
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    if (inInfo->Has(vtkStreamingDemandDrivenPipeline::TIME_STEPS()))
    {
      // reset the CurrentTimeIndex in case we're writing out all of the time steps
      this->CurrentTimeIndex = 0;
      this->NumberOfTimeSteps = inInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    }
    else
    {
      this->NumberOfTimeSteps = 1;
    }
  }
  else if (request->Has(vtkStreamingDemandDrivenPipeline::REQUEST_DATA()))
  {
    if (this->WriteAllTimeSteps && this->CurrentTimeIndex == 0)
    {
      // Tell the pipeline to start looping.
      request->Set(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1);
    }
  }

  int retVal = this->Superclass::ProcessRequest(request, inputVector, outputVector);

  if (request->Has(vtkStreamingDemandDrivenPipeline::REQUEST_DATA()))
  {
    if (this->WriteAllTimeSteps && this->CurrentTimeIndex == this->NumberOfTimeSteps)
    {
      // Tell the pipeline to stop looping.
      request->Remove(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING());
      this->CurrentTimeIndex = 0;
    }
  }

  return retVal;
}



//-----------------------------------------------------------------------------
void vtkVDBWriter::WriteData()
{
  // double time = vtkMath::Nan();
  // if (this->AddTime && input && input->GetInformation()->Has(vtkDataObject::DATA_TIME_STEP()))
  // {
  //   time = input->GetInformation()->Get(vtkDataObject::DATA_TIME_STEP());
  // }


  if (vtkImageData* imageData = vtkImageData::SafeDownCast(this->GetInput()))
  {
    this->WriteImageData(imageData);
  }
  else if (vtkPointSet* pointSet = vtkPointSet::SafeDownCast(this->GetInput()))
  {
    this->WritePointSet(pointSet);
  }
  this->CurrentTimeIndex++;


}

//-----------------------------------------------------------------------------
void vtkVDBWriter::WriteImageData(vtkImageData* imageData)
{



  std::vector<openvdb::GridBase::Ptr> grids;

  double dx(0), dy(0), dz(0);
  imageData->GetSpacing(dx, dy, dz);

  openvdb::math::Mat4d mat(dx, 0., 0., 0.,
                           0., dy, 0., 0.,
                           0., 0., dz, 0.,
                           0., 0., 0., 1.);

  openvdb::math::Transform::Ptr linearTransform =
    openvdb::math::Transform::createLinearTransform(mat);

  int extent[6], wholeExtent[6];
  imageData->GetExtent(extent);
  this->Controller->AllReduce(extent, wholeExtent, 6, vtkCommunicator::MAX_OP);
  imageData->GetExtent(extent);
  int pointExtent[6] = {extent[0], extent[1], extent[2], extent[3], extent[4], extent[5]};

  // since we don't want duplicate data in parallel for the pointdata
  // we chop of the upper points if they are not on the boundary and let
  // other processes handle that data
  for (int i=0;i<3;i++)
  {
    if (extent[2*i+1] != wholeExtent[2*i+1])
    {
      pointExtent[2*i+1] = extent[2*i+1]-1;
    }
  }

  openvdb::Coord denseDim(pointExtent[1]-pointExtent[0]+1, pointExtent[3]-pointExtent[2]+1, pointExtent[5]-pointExtent[4]+1);
  openvdb::Coord denseMin(pointExtent[0], pointExtent[2], pointExtent[4]);

  cerr << this->Controller->GetLocalProcessId() << " pointextents " << pointExtent[0] << " " << pointExtent[1] << " " << pointExtent[2] << " " << pointExtent[3] << " " << pointExtent[4] << " " << pointExtent[5] << endl;

  // compute colors, if any
  vtkNew<vtkPointData> pointData;
  pointData->ShallowCopy(imageData->GetPointData());
  vtkNew<vtkCellData> cellData;
  cellData->ShallowCopy(imageData->GetCellData());
  if (this->LookupTable && this->EnableColoring)
  {
    int fieldAssociation = 0;
    vtkAbstractArray* scalars = this->GetInputAbstractArrayToProcess(0, imageData, fieldAssociation);
    vtkDiscretizableColorTransferFunction* dctf =
      vtkDiscretizableColorTransferFunction::SafeDownCast(this->LookupTable);
    bool enableOpacityMapping = dctf->GetEnableOpacityMapping();
    dctf->SetEnableOpacityMapping(this->EnableAlpha);
    vtkUnsignedCharArray* rgba =
      this->LookupTable->MapScalars(scalars, VTK_COLOR_MODE_MAP_SCALARS, -1);
    vtkIdType numPoints = imageData->GetNumberOfPoints();
    if (rgba && rgba->GetNumberOfTuples() == numPoints)
    {
      this->SetRGBA(numPoints, rgba, pointData);
      // if (pointColors)
      // {
      //   pointColors->SetName("color");
      //   pointData->AddArray(pointColors);
      // }
    }
    fieldAssociation = 1;
    scalars = this->GetInputAbstractArrayToProcess(0, imageData, fieldAssociation);
    dctf = vtkDiscretizableColorTransferFunction::SafeDownCast(this->LookupTable);
    rgba = this->LookupTable->MapScalars(scalars, VTK_COLOR_MODE_MAP_SCALARS, -1);
    vtkIdType numCells = imageData->GetNumberOfCells();
    if (rgba && rgba->GetNumberOfTuples() == numCells)
    {
      this->SetRGBA(numCells, rgba, cellData);
      // if (cellColors)
      // {
      //   cellColors->SetName("color");
      //   cellData->AddArray(cellColors);
      // }
    }
    dctf->SetEnableOpacityMapping(enableOpacityMapping);
  }



  for (int array=0;array<pointData->GetNumberOfArrays();array++)
  {
    bool useDenseGrid = false;
    vtkDataArray* data = pointData->GetArray(array);
    const char* arrayName = data->GetName();
    const int numberOfComponents = data->GetNumberOfComponents();
    for (int component=0;component<numberOfComponents;component++)
    {
      if (numberOfComponents == 3 && component > 0)
      {
        continue;
      }
      // Vec3SGrid is sparse but for image data want Vec3DGrid for dense storage for
      // smaller size
      openvdb::Vec3DGrid::Ptr vecGrid = openvdb::Vec3DGrid::create();
      // see https://www.openvdb.org/documentation/doxygen/namespaceopenvdb_1_1v8__0.html#ae93f92d10730a52ed3b207d5811f6a6e
      if (strcmp(arrayName, "color"))
      {
        vecGrid->setVectorType(openvdb::VEC_CONTRAVARIANT_RELATIVE);
      }
      else
      {
        vecGrid->setVectorType(openvdb::VEC_INVARIANT);
      }
      vecGrid->setGridClass(openvdb::GRID_FOG_VOLUME);


      openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create();
      // enum GridClass {
      //     GRID_UNKNOWN = 0,
      //     GRID_LEVEL_SET,
      //     GRID_FOG_VOLUME,
      //     GRID_STAGGERED
      // };
      grid->setGridClass(openvdb::GRID_FOG_VOLUME);
      //auto denseGrid = new openvdb::tools::Dense(denseDim, denseMin);
      // openvdb::tools::Dense<float>::Ptr pp;
      // auto grid222 = openvdb::FloatGrid::createGrid<openvdb::tools::Dense<float> >();
      // const openvdb::CoordBBox bbox(openvdb::Coord(pointExtent[0], pointExtent[2], pointExtent[4]),
      //                               openvdb::Coord(pointExtent[1]+1, pointExtent[3]+1, pointExtent[5]+1));
      // openvdb::tools::Dense<float>* d = new openvdb::tools::Dense<float>(bbox);
      // //openvdb::tools::Dense<float>* d = new openvdb::tools::Dense<float>(bbox);
      // //openvdb::tools::Dense<float>::Ptr denseGrid(d);
      // auto denseGrid = std::make_shared<openvdb::tools::Dense<float>* >(d);

      std::string vdbName = GetVDBGridName(arrayName, component, numberOfComponents);
      grid->setName(vdbName);
      vecGrid->setName(vdbName);
      //denseGrid->setName(vdbName);
      openvdb::FloatGrid::Accessor accessor = grid->getAccessor();
      openvdb::Vec3DGrid::Accessor vecAccessor = vecGrid->getAccessor();
      //auto denseAccessor = denseGrid->getAccessor();
      openvdb::Coord ijk;

      int &i = ijk[0], &j = ijk[1], &k = ijk[2];
      int vtkijk[3];
      for(k = pointExtent[4]; k <= pointExtent[5]; ++k)
      {
        vtkijk[2] = k;
        for(j = pointExtent[2]; j <= pointExtent[3]; ++j)
        {
          vtkijk[1] = j;
          for(i = pointExtent[0]; i <= pointExtent[1]; ++i)
          {
            vtkijk[0] = i;
            if (numberOfComponents == 3)
            {
              double tuple[3];
              data->GetTuple(imageData->ComputePointId(vtkijk), tuple);
              openvdb::Vec3f ftuple(tuple[0], tuple[1], tuple[2]);
              vecAccessor.setValue(ijk, ftuple);
            }
            else
            {
              if (useDenseGrid)
              {
                //denseAccessor.setValue(ijk, data->GetComponent(imageData->ComputePointId(vtkijk), component));
              }
              else
              {
                accessor.setValue(ijk, data->GetComponent(imageData->ComputePointId(vtkijk), component));
              }
            }
          }
        }
      grid->setTransform(linearTransform);
      //denseGrid->setTransform(linearTransform);
      vecGrid->setTransform(linearTransform);
      }
      if (numberOfComponents == 3)
      {
        grids.push_back(vecGrid);
      }
      else
      {
        if (useDenseGrid)
        {
          //grids.push_back(denseGrid);
        }
        else
        {
          grids.push_back(grid);
        }
      }
    }
  }

  for (int array=0;array<cellData->GetNumberOfArrays();array++)
  {
    vtkDataArray* data = cellData->GetArray(array);
    const char* arrayName = data->GetName();
    openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create();
// enum GridClass {
//     GRID_UNKNOWN = 0,
//     GRID_LEVEL_SET,
//     GRID_FOG_VOLUME,
//     GRID_STAGGERED
// };
    grid->setGridClass(openvdb::GRID_FOG_VOLUME);
    grid->setName(arrayName);
    openvdb::FloatGrid::Accessor accessor = grid->getAccessor();
    openvdb::Coord ijk;

    vtkIdType counter = 0;
    int &i = ijk[0], &j = ijk[1], &k = ijk[2];
    // for cell data we don't have to worry about ghost cells. hopefully they're gone in parallel. -- check this ACB
    for(k = extent[4]; k < extent[5]; ++k)
    {
      for(j = extent[2]; j < extent[3]; ++j)
      {
        for(i = extent[0]; i < extent[1]; ++i)
        {
          //cerr << this->Controller->GetLocalProcessId() << " is writing " << i << " " << j << " " << k << " " << data->GetTuple1(counter) << endl;
          accessor.setValue(ijk, data->GetTuple1(counter));
          counter++;
        }
      }


      grid->setTransform(linearTransform);

    }
    grids.push_back(grid);
  }
  // if no point data or no cell data then we just add in a value of 1 for the voxels for the cells
  if (grids.empty())
  {
    openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create();
// enum GridClass {
//     GRID_UNKNOWN = 0,
//     GRID_LEVEL_SET,
//     GRID_FOG_VOLUME,
//     GRID_STAGGERED
// };
    grid->setGridClass(openvdb::GRID_FOG_VOLUME);
    grid->setName("blank"); // ACB find a better name
    openvdb::FloatGrid::Accessor accessor = grid->getAccessor();
    openvdb::Coord ijk;

    int &i = ijk[0], &j = ijk[1], &k = ijk[2];
    // for cell data we don't have to worry ghost cells. hopefully they're gone in parallel. -- check this ACB
    for(k = extent[4]; k < extent[5]; ++k)
    {
      for(j = extent[2]; j < extent[3]; ++j)
      {
        for(i = extent[0]; i < extent[1]; ++i)
        {
          accessor.setValue(ijk, 1.);
        }
      }

      grid->setTransform(linearTransform);

    }
    grids.push_back(grid);


  }

  //openvdb::io::File(output.data()).write({grid});
  //openvdb::io::File(this->FileName).write({grid});
  std::ostringstream oss;
  oss << std::setw(5) << std::setfill('0') << this->CurrentTimeIndex;
  if (this->Controller->GetNumberOfProcesses() == 1)
  {
    if (this->WriteAllTimeSteps && this->NumberOfTimeSteps > 1)
    {
      std::string fileNameBase = vtksys::SystemTools::GetFilenameWithoutExtension(this->FileName);
      std::string ext = vtksys::SystemTools::GetFilenameExtension(this->FileName);
      std::string newFileName = fileNameBase+"_"+oss.str()+ext;
      openvdb::io::File(newFileName).write(grids);
    }
    else
    {
      openvdb::io::File(this->FileName).write(grids);
    }
  }
  else
  {
    if (this->WriteAllTimeSteps && this->NumberOfTimeSteps > 1)
    {
      std::string fileNameBase = vtksys::SystemTools::GetFilenameWithoutExtension(this->FileName);
      std::string ext = vtksys::SystemTools::GetFilenameExtension(this->FileName);
      std::string newFileName = fileNameBase+"_"+std::to_string(this->Controller->GetLocalProcessId())+
        "_"+oss.str()+ext;
      openvdb::io::File(newFileName).write(grids);
    }
    else
    {
      std::string fileNameBase = vtksys::SystemTools::GetFilenameWithoutExtension(this->FileName);
      std::string ext = vtksys::SystemTools::GetFilenameExtension(this->FileName);
      std::string newFileName = fileNameBase+"_"+std::to_string(this->Controller->GetLocalProcessId())+ext;
      openvdb::io::File(newFileName).write(grids);
    }
  }
}

//-----------------------------------------------------------------------------
void vtkVDBWriter::WritePointSet(vtkPointSet* pointSet)
{
  // openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create();
  // openvdb::FloatGrid::Accessor accessor = grid->getAccessor();
  // openvdb::Coord xyz(1000, -200000000, 30000000); // ACB ???

  std::vector<openvdb::Vec3R> positions;
  double coord[3];
  for (vtkIdType i=0;i<pointSet->GetNumberOfPoints();i++)
  {
    pointSet->GetPoint(i, coord);
    positions.push_back(openvdb::Vec3R(static_cast<float>(coord[0]),
                                       static_cast<float>(coord[1]),
                                       static_cast<float>(coord[2])));
  }

  // The VDB Point-Partioner is used when bucketing points and requires a
  // specific interface. For convenience, we use the PointAttributeVector
  // wrapper around an stl vector wrapper here, however it is also possible to
  // write one for a custom data structure in order to match the interface
  // required.
  openvdb::points::PointAttributeVector<openvdb::Vec3R> positionsWrapper(positions);
  // This method computes a voxel-size to match the number of
  // points / voxel requested. Although it won't be exact, it typically offers
  // a good balance of memory against performance.
  int pointsPerVoxel = 8;
  float voxelSize =
    openvdb::points::computeVoxelSize(positionsWrapper, pointsPerVoxel);
  // Create a transform using this voxel-size.
  openvdb::math::Transform::Ptr transform =
    openvdb::math::Transform::createLinearTransform(voxelSize);
  // Create a PointIndexGrid. This can be done automatically on creation of
  // the grid, however as this index grid is required for the position and
  // radius attributes, we create one we can use for both attribute creation.
  openvdb::tools::PointIndexGrid::Ptr pointIndexGrid =
    openvdb::tools::createPointIndexGrid<openvdb::tools::PointIndexGrid>(
      positionsWrapper, *transform);
  // Create a PointDataGrid containing these four points and using the point
  // index grid. This requires the positions wrapper.
  openvdb::points::PointDataGrid::Ptr grid =
    openvdb::points::createPointDataGrid<openvdb::points::NullCodec,
                                         openvdb::points::PointDataGrid>(*pointIndexGrid, positionsWrapper, *transform);

  // Set the name of the grid
  grid->setName("Points");



  for (int array=0;array<pointSet->GetPointData()->GetNumberOfArrays();array++)
  {
    std::vector<float> values;
    vtkDataArray* data = pointSet->GetPointData()->GetArray(array);
    const char* arrayName = data->GetName();
    for (vtkIdType i=0;i<pointSet->GetNumberOfPoints();i++)
    {
      values.push_back(static_cast<float>(data->GetTuple1(i)));
    }

    // Append a data->GetName() attribute to the grid to hold the radius.
    // This attribute storage uses a unit range codec to reduce the memory
    // storage requirements down from 4-bytes to just 1-byte per value. This is
    // only possible because accuracy of the radius is not that important to us
    // and the values are always within unit range (0.0 => 1.0).
    // Note that this attribute type is not registered by default so needs to be
    // explicitly registered.
    // using Codec = openvdb::points::FixedPointCodec</*1-byte=*/false,
    //         openvdb::points::UnitRange>;
    // openvdb::points::TypedAttributeArray<float, Codec>::registerType();
    openvdb::NamePair radiusAttribute =
      openvdb::points::TypedAttributeArray<float, openvdb::points::NullCodec>::attributeType();
    openvdb::points::appendAttribute(grid->tree(), data->GetName(), radiusAttribute);
    // Create a wrapper around the values vector.
    openvdb::points::PointAttributeVector<float> valuesWrapper(values);
    // Populate the data->GetName() attribute on the points
    openvdb::points::populateAttribute<openvdb::points::PointDataTree,
        openvdb::tools::PointIndexTree, openvdb::points::PointAttributeVector<float>>(
          grid->tree(), pointIndexGrid->tree(), data->GetName(), valuesWrapper);



  }
  if (this->Controller->GetNumberOfProcesses() == 1)
  {
    openvdb::io::File(this->FileName).write({grid});
  }
  else
  {
    std::string fileNameBase = vtksys::SystemTools::GetFilenameWithoutExtension(this->FileName);
    std::string ext = vtksys::SystemTools::GetFilenameExtension(this->FileName);
    std::string newFileName = fileNameBase+"_"+std::to_string(this->Controller->GetLocalProcessId())+ext;
    openvdb::io::File(newFileName).write({grid});
  }
}

//-----------------------------------------------------------------------------
void vtkVDBWriter::SetRGBA(
  vtkIdType num, vtkUnsignedCharArray* rgbaArray, vtkDataSetAttributes* attributes)
{
  vtkIdType i;

  int numComp = rgbaArray->GetNumberOfComponents();
  if (numComp == 3)
  { // have unsigned char array of three components, copy it
    rgbaArray->SetName("color");
    attributes->AddArray(rgbaArray);
    if (this->EnableAlpha)
    {
      vtkWarningMacro("No alpha channel to set even though requested");
    }
    return;
  }
  else if (numComp == 4)
  {
    // have unsigned char array of four components (RGBA), copy it without the `A`.
    vtkSmartPointer<vtkFloatArray> colors = vtkSmartPointer<vtkFloatArray>::New();
    colors->SetNumberOfComponents(3);
    colors->SetNumberOfTuples(num);
    colors->SetName("color");
    float* c = colors->WritePointer(0, 3 * num);
    vtkSmartPointer<vtkFloatArray> alpha = vtkSmartPointer<vtkFloatArray>::New();
    alpha->SetNumberOfComponents(1);
    alpha->SetNumberOfTuples(num);
    alpha->SetName("alpha");
    float* a = alpha->WritePointer(0, num);
    const unsigned char* rgba = rgbaArray->GetPointer(0);
    for (i = 0; i < num; i++)
    {
      *c++ = (*rgba++)/255.;
      *c++ = (*rgba++)/255.;
      *c++ = (*rgba++)/255.;
      *a++ = (*rgba++)/255.;
    }
    attributes->AddArray(colors);
    if (this->EnableAlpha)
    {
      attributes->AddArray(alpha);
    }
    return;
  }
  // else if (this->LookupTable != nullptr)
  // { // use the data array mapped through lookup table
  //   vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
  //   colors->SetNumberOfComponents(3);
  //   colors->SetNumberOfTuples(num);
  //   c = colors->WritePointer(0, 3 * num);
  //   for (i = 0; i < num; i++)
  //   {
  //     double* tuple = da->GetTuple(i);
  //     const unsigned char* rgb = this->LookupTable->MapValue(tuple[this->Component]);
  //     *c++ = rgb[0];
  //     *c++ = rgb[1];
  //     *c++ = rgb[2];
  //   }
  //   return colors;
  // }
  // no lookup table
  return;
}

//-----------------------------------------------------------------------------
void vtkVDBWriter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "FileName: " << (this->FileName ? this->FileName : "none") << endl;
  os << indent << "WriteAllTimeSteps: " << this->WriteAllTimeSteps << endl;
  if (this->Controller)
  {
    os << indent << "Controller: " << this->Controller << endl;
  }
  else
  {
    os << indent << "Controller: (none)" << endl;
  }
}
