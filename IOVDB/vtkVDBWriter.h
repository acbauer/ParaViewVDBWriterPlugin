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
 * Writes a vtkImageData or vtkPointSet as a VDB file.
*/

#ifndef vtkVDBWriter_h
#define vtkVDBWriter_h

#include "vtkVDBWritersModule.h" //needed for exports
#include "vtkWriter.h"

class vtkImageData;
class vtkPointSet;
class vtkMultiProcessController;

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

  char* FileName;

private:
  vtkVDBWriter(const vtkVDBWriter&) = delete;
  void operator=(const vtkVDBWriter&) = delete;

  vtkMultiProcessController* Controller;

  //class VDBFile;
};

#endif
