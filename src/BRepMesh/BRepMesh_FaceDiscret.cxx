// Created on: 2016-04-19
// Copyright (c) 2016 OPEN CASCADE SAS
// Created by: Oleg AGASHIN
//
// This file is part of Open CASCADE Technology software library.
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License version 2.1 as published
// by the Free Software Foundation, with special exception defined in the file
// OCCT_LGPL_EXCEPTION.txt. Consult the file LICENSE_LGPL_21.txt included in OCCT
// distribution for complete text of the license and disclaimer of any warranty.
//
// Alternatively, this file may be used under the terms of Open CASCADE
// commercial license or contractual agreement.

#include <BRepMesh_FaceDiscret.hxx>
#include <IMeshData_Model.hxx>
#include <IMeshData_Face.hxx>
#include <IMeshData_Wire.hxx>
#include <IMeshData_Edge.hxx>
#include <IMeshTools_MeshAlgo.hxx>
#include <OSD_Parallel.hxx>

#include <map>
#include <string>

IMPLEMENT_STANDARD_RTTIEXT(BRepMesh_FaceDiscret, IMeshTools_ModelAlgo)

//=======================================================================
// Function: Constructor
// Purpose : 
//=======================================================================
BRepMesh_FaceDiscret::BRepMesh_FaceDiscret(
  const Handle(IMeshTools_MeshAlgoFactory)& theAlgoFactory)
  : myAlgoFactory(theAlgoFactory)
{
}

//=======================================================================
// Function: Destructor
// Purpose : 
//=======================================================================
BRepMesh_FaceDiscret::~BRepMesh_FaceDiscret()
{
}

//! Auxiliary functor for parallel processing of Faces.
class BRepMesh_FaceDiscret::FaceListFunctor
{
public:
  FaceListFunctor (BRepMesh_FaceDiscret* theAlgo,
                   const Message_ProgressRange& theRange)
  : myAlgo (theAlgo),
    myScope (theRange, "Face Discret", theAlgo->myModel->FacesNb())
  {
    myRanges.reserve (theAlgo->myModel->FacesNb());
    for (Standard_Integer aFaceIter = 0; aFaceIter < theAlgo->myModel->FacesNb(); ++aFaceIter)
    {
      myRanges.push_back (myScope.Next());
    }
  }

  void operator() (const Standard_Integer theFaceIndex) const
  {
    if (!myScope.More())
    {
      return;
    }
    Message_ProgressScope aFaceScope(myRanges[theFaceIndex], NULL, 1);
    myAlgo->process(theFaceIndex, aFaceScope.Next());
  }

private:
  mutable BRepMesh_FaceDiscret* myAlgo;
  Message_ProgressScope myScope;
  std::vector<Message_ProgressRange> myRanges;
};

//=======================================================================
// Function: Perform
// Purpose : 
//=======================================================================
Standard_Boolean BRepMesh_FaceDiscret::performInternal(
  const Handle(IMeshData_Model)& theModel,
  const IMeshTools_Parameters&   theParameters,
  const Message_ProgressRange&   theRange)
{
  myModel      = theModel;
  myParameters = theParameters;
  if (myModel.IsNull())
  {
    return Standard_False;
  }

  FaceListFunctor aFunctor(this, theRange);
  OSD_Parallel::For(0, myModel->FacesNb(), aFunctor, Standard_True);
  // OSD_Parallel::For(0, myModel->FacesNb(), aFunctor, !(myParameters.InParallel && myModel->FacesNb() > 1));
  if (!theRange.More())
  {
    return Standard_False;
  }

  myModel.Nullify(); // Do not hold link to model.
  return Standard_True;
}

std::string GetFaceTypeString(const GeomAbs_SurfaceType &facetype)
{
  static std::map<GeomAbs_SurfaceType, std::string> facetype2string;
  if (facetype2string.size() == 0)
  {
    facetype2string[GeomAbs_Plane] = "GeomAbs_Plane";
    facetype2string[GeomAbs_Cylinder] = "GeomAbs_Cylinder";
    facetype2string[GeomAbs_Cone] = "GeomAbs_Cone";
    facetype2string[GeomAbs_Sphere] = "GeomAbs_Sphere";
    facetype2string[GeomAbs_Torus] = "GeomAbs_Torus";
    facetype2string[GeomAbs_BezierSurface] = "GeomAbs_BezierSurface";
    facetype2string[GeomAbs_BSplineSurface] = "GeomAbs_BSplineSurface";
    facetype2string[GeomAbs_SurfaceOfRevolution] = "GeomAbs_SurfaceOfRevolution";
    facetype2string[GeomAbs_SurfaceOfExtrusion] = "GeomAbs_SurfaceOfExtrusion";
    facetype2string[GeomAbs_OffsetSurface] = "GeomAbs_OffsetSurface";
    facetype2string[GeomAbs_OtherSurface] = "GeomAbs_OtherSurface";
  }
  return facetype2string[facetype];
}

//=======================================================================
// Function: process
// Purpose : 
//=======================================================================
void BRepMesh_FaceDiscret::process(const Standard_Integer theFaceIndex, 
                                   const Message_ProgressRange& theRange) const
{
  const IMeshData::IFaceHandle& aDFace = myModel->GetFace(theFaceIndex);
  if (aDFace->IsSet(IMeshData_Failure) ||
      aDFace->IsSet(IMeshData_Reused))
  {
    return;
  }

  try
  {
    OCC_CATCH_SIGNALS

    Handle(IMeshTools_MeshAlgo) aMeshingAlgo = 
      myAlgoFactory->GetAlgo(aDFace->GetSurface()->GetType(), myParameters);
  
    if (aMeshingAlgo.IsNull())
    {
      aDFace->SetStatus(IMeshData_Failure);
      return;
    }
  
    if (!theRange.More())
    {
      aDFace->SetStatus (IMeshData_UserBreak);
      return;
    }
    GeomAbs_SurfaceType surfacetype;
    const TopoDS_Face &oceface = aDFace->GetFace();
    if (!BRep_Tool::Surface(oceface).IsNull())
    { // Get the surface type
      BRepAdaptor_Surface First(oceface, Standard_False);
      surfacetype = First.GetType();
      std::cout << "Begin " << GetFaceTypeString(surfacetype) << std::endl;
    }
    aMeshingAlgo->Perform(aDFace, myParameters, theRange);
    std::cout << "End " << GetFaceTypeString(surfacetype) << std::endl;
  }
  catch (Standard_Failure const&)
  {
    aDFace->SetStatus (IMeshData_Failure);
  }
}
