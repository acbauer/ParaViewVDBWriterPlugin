/*=========================================================================

  Program:   ParaView
  Module:    vtkVDBWriter.h

  Copyright (c) Kitware, Inc.
  All rights reserved.
  See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class   vtkVDBWriter
 * @brief   VDB writer for vtkImageData or vtkPointSet
 * Writes a vtkImageData or vtkPointSEt as a VDB file.
*/

#ifndef vtkVDBWriter_h
#define vtkVDBWriter_h

#include "vtkVDBWritersModule.h"
#include "vtkSmartPointer.h" // For protected ivars
#include "vtkWriter.h"
#include <string>

class vtkDataSetAttributes;
class vtkImageData;
class vtkMultiProcessController;
class vtkPointSet;
class vtkScalarsToColors;
class vtkUnsignedCharArray;
class vtkVDBWriterInternals;

class VTKVDBWRITERS_EXPORT vtkVDBWriter : public vtkWriter
{
public:
  static vtkVDBWriter* New();
  vtkTypeMacro(vtkVDBWriter, vtkWriter);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  //@{
  /**
   * Get/Set the filename for the file.
   */
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);
  //@}

  //@{
  /**
   * Get/Set whether or not to save all time steps or
   * just the current time step. Default is false
   * (save only the current time step).
   */
  vtkSetMacro(WriteAllTimeSteps, bool);
  vtkGetMacro(WriteAllTimeSteps, bool);
  //@}

  //@{
  /**
   * A lookup table can be specified in order to convert data arrays to
   * RGBA colors.
   */
  virtual void SetLookupTable(vtkScalarsToColors*);
  vtkGetObjectMacro(LookupTable, vtkScalarsToColors);
  //@}

  //@{
  /**
   * Enable coloring channel output based on LookupTable. The output
   * channel will be named 'color'.
   */
  vtkSetMacro(EnableColoring, bool);
  vtkGetMacro(EnableColoring, bool);
  //@}

  //@{
  /**
   * Enable alpha channel output based on LookupTable. The output
   * channel will be name 'alpha'.
   */
  vtkSetMacro(EnableAlpha, bool);
  vtkGetMacro(EnableAlpha, bool);
  //@}

  //@{
  /**
   * Get/Set the controller to use. By default,
   * `vtkMultiProcessController::GetGlobalController` will be used.
   */
  void SetController(vtkMultiProcessController*);
  vtkGetObjectMacro(Controller, vtkMultiProcessController);
  //@}

protected:
  vtkVDBWriter();
  ~vtkVDBWriter() override;

  void WriteData() override;

  void WriteImageData(vtkImageData* imageData);
  void WritePointSet(vtkPointSet* pointSet);

  // see algorithm for more info.
  // This writer takes in vtkTable, vtkDataSet or vtkCompositeDataSet.
  int FillInputPortInformation(int port, vtkInformation* info) override;

  // see algorithm for more info. needed here so we can request pieces.
  int ProcessRequest(vtkInformation* request, vtkInformationVector** inputVector,
    vtkInformationVector* outputVector) override;

  void SetRGBA(vtkIdType num, vtkUnsignedCharArray* rgba, vtkDataSetAttributes* attributes);

  char* FileName;

  bool WriteAllTimeSteps;

  // For outputting the Lookup Table in the VDB file.
  // Copying what's done in vtkPLYWriter
  vtkScalarsToColors* LookupTable;
  bool EnableColoring;
  bool EnableAlpha;


private:
  vtkVDBWriter(const vtkVDBWriter&) = delete;
  void operator=(const vtkVDBWriter&) = delete;

  vtkMultiProcessController* Controller;
  // Internal variable to keep track of the current time step index
  // when writing out all the time steps.
  vtkIdType CurrentTimeIndex;
  // Internal variable to keep track of the number of time steps.
  vtkIdType NumberOfTimeSteps;

  vtkVDBWriterInternals* Internals;
  friend class vtkVDBWriterInternals;

  //class VDBFile;
};

#endif
