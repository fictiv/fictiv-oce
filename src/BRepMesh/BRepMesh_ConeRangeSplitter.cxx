// Created on: 2016-07-07
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

#include <BRepMesh_ConeRangeSplitter.hxx>
#include <GCPnts_TangentialDeflection.hxx>
#include <set>
#include <algorithm>
#include <iostream>

//=======================================================================
// Function: GetSplitSteps
// Purpose : 
//=======================================================================
std::pair<Standard_Real, Standard_Real> BRepMesh_ConeRangeSplitter::GetSplitSteps(
  const IMeshTools_Parameters&                   theParameters,    
  std::pair<Standard_Integer, Standard_Integer>& theStepsNb) const
{
  const std::pair<Standard_Real, Standard_Real>& aRangeU = GetRangeU();
  const std::pair<Standard_Real, Standard_Real>& aRangeV = GetRangeV();

  gp_Cone aCone = GetDFace()->GetSurface()->Cone();
  Standard_Real aRefR = aCone.RefRadius();
  Standard_Real aSAng = aCone.SemiAngle();
  Standard_Real aRadius = Max(Abs(aRefR + aRangeV.first  * Sin(aSAng)),
                              Abs(aRefR + aRangeV.second * Sin(aSAng)));

  Standard_Real Dv, Du = GCPnts_TangentialDeflection::ArcAngularStep(
    aRadius, GetDFace()->GetDeflection(),
    theParameters.Angle, theParameters.MinSize);

  const Standard_Real aDiffU = aRangeU.second - aRangeU.first;
  const Standard_Real aDiffV = aRangeV.second - aRangeV.first;
  Standard_Integer nbU = (Standard_Integer) (aDiffU / Du);
  Standard_Integer nbV = (Standard_Integer) (nbU * (aDiffV) / (aDiffU * aRadius));
  Du = aDiffU / (nbU + 1);
  Dv = aDiffV / (nbV + 1);

  theStepsNb.first  = nbU;
  theStepsNb.second = nbV;
  return std::make_pair (Du, Dv);
}

//=======================================================================
// Function: GenerateSurfaceNodes
// Purpose : 
//=======================================================================
Handle(IMeshData::ListOfPnt2d) BRepMesh_ConeRangeSplitter::GenerateSurfaceNodes(
  const IMeshTools_Parameters& theParameters) const
{
  const std::pair<Standard_Real, Standard_Real>& aRangeU = GetRangeU();
  const std::pair<Standard_Real, Standard_Real>& aRangeV = GetRangeV();

  std::pair<Standard_Integer, Standard_Integer> aStepsNb;
  std::pair<Standard_Real, Standard_Real> aSteps = GetSplitSteps (theParameters, aStepsNb);

  const Handle(NCollection_IncAllocator) aTmpAlloc =
    new NCollection_IncAllocator(IMeshData::MEMORY_BLOCK_SIZE_HUGE);
  Handle(IMeshData::ListOfPnt2d) aNodes = new IMeshData::ListOfPnt2d(aTmpAlloc);


  const Standard_Real aPasMaxV = aRangeV.second - aSteps.second*0.5;
  const Standard_Real aPasMaxU = aRangeU.second - aSteps.first *0.5;
  std::vector<Standard_Real> vpnts;
  std::vector<Standard_Real> upnts;
  for(auto&& pnt : boundaryPnts){
    upnts.push_back(pnt.X());
    vpnts.push_back(pnt.Y());
  }
  std::sort(upnts.begin(), upnts.end());
  std::sort(vpnts.begin(), vpnts.end());
  if(vpnts.size()<1 || upnts.size()<1){
    return aNodes;
  }

  //for (Standard_Real aPasV = aRangeV.first + aSteps.second; aPasV < aPasMaxV; aPasV += aSteps.second){
  Standard_Real aPasVprev = vpnts[0];
  for (unsigned int i=0; i<vpnts.size(); i++){
    Standard_Real aPasV = vpnts[i];
    if(i>0){
      Standard_Real deltaV = aPasV-aPasVprev;
      if(deltaV==0){
        continue;
      }
      else if(false && deltaV>aSteps.second){ //I'm not sure how this would happen given that the boundary sampling should follow the same logic
        //Insert an intermediate point
        aPasV = aPasVprev + aSteps.second;
        i--;  //We'll try this boundary point again
      }
    }
    //for (Standard_Real aPasU = aRangeU.first + aSteps.first; aPasU < aPasMaxU; aPasU += aSteps.first){
    Standard_Real aPasUprev = upnts[0];
    for(unsigned int j=0; j<upnts.size(); j++){
      Standard_Real aPasU = upnts[j];
      Standard_Real deltaU = aPasU-aPasUprev;
      if(j>0){
        if(deltaU==0){
          continue;
        }
        else if(deltaU>aSteps.first){ //I'm not sure how this would happen given that the boundary sampling should follow the same logic
          //Insert an intermediate point
          aPasU = aPasUprev + aSteps.first;
          j--;  //We'll try this boundary point again
        }
      }

      aNodes->Append(gp_Pnt2d(aPasU, aPasV));

      aPasUprev = aPasU;
    }
    aPasVprev = aPasV;
  }

  return aNodes;
}

void BRepMesh_ConeRangeSplitter::AddPoint(const gp_Pnt2d &thePoint)
{
  std::cout<<"Using BRepMesh_ConeRangeSplitter::AddPoint"<<std::endl;
  BRepMesh_DefaultRangeSplitter::AddPoint(thePoint);
  boundaryPnts.push_back(thePoint);
}

