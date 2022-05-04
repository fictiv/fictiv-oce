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

#include <BRepMesh_NURBSRangeSplitter.hxx>

#include <algorithm>
#include <BRepMesh_GeomTool.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <GeomLib.hxx>
#include <IMeshData_Edge.hxx>
#include <IMeshData_Wire.hxx>
#include <NCollection_Handle.hxx>
#include <Geom_BezierSurface.hxx>
#include <Geom_BSplineSurface.hxx>
#include <GCPnts_AbscissaPoint.hxx>

#include <iostream>

namespace
{
  class AnalyticalFilter
  {
  public:
    //! Constructor.
    AnalyticalFilter(
      const IMeshData::IFaceHandle&             theDFace,
      const GeomAbs_IsoType                     theIsoType,
      const Handle(IMeshData::SequenceOfReal)&  theParams,
      const Handle(IMeshData::SequenceOfReal)&  theControlParams,
      const Handle(IMeshData::MapOfReal)&       theParamsForbiddenToRemove,
      const Handle(IMeshData::MapOfReal)&       theControlParamsForbiddenToRemove)
      : myDFace(theDFace),
        mySurface(myDFace->GetSurface()->Surface().Surface()),
        myIsoU(theIsoType == GeomAbs_IsoU),
        myParams(theParams),
        myControlParams(theControlParams),
        myParamsForbiddenToRemove(theParamsForbiddenToRemove),
        myControlParamsForbiddenToRemove(theControlParamsForbiddenToRemove),
        myAllocator(new NCollection_IncAllocator(IMeshData::MEMORY_BLOCK_SIZE_HUGE)),
        myControlParamsToRemove(new IMeshData::MapOfReal(1, myAllocator)),
        myCurrParam(0.0),
        myCurrControlParam(0.0),
        myPrevControlParam(0.0)
    {
    }

    //! Returns map of parameters supposed to be removed.
    const Handle(IMeshData::MapOfReal)& GetControlParametersToRemove(
      const IMeshTools_Parameters& theParameters)
    {
      myParameters = theParameters;

      Standard_Integer aStartIndex, aEndIndex;
      if (myIsoU)
      {
        aStartIndex = 1;
        aEndIndex = myParams->Length();
      }
      else
      {
        aStartIndex = 2;
        aEndIndex = myParams->Length() - 1;
      }

      for (Standard_Integer i = aStartIndex; i <= aEndIndex; ++i)
      {
        myCurrParam = myParams->Value(i);
        myIso = new GeomAdaptor_Curve(myIsoU ? mySurface->UIso(myCurrParam) : mySurface->VIso(myCurrParam));

        myPrevControlParam = myControlParams->Value(1);
        myIso->D1(myPrevControlParam, myPrevControlPnt, myPrevControlVec);
        for (Standard_Integer j = 2; j <= myControlParams->Length();)
        {
          j += checkControlPointAndMoveOn(j);
        }
      }

      return myControlParamsToRemove;
    }

  private:

    //! Checks the given control point for deviation.
    //! Returns number of steps to be used to move point iterator.
    Standard_Integer checkControlPointAndMoveOn(const Standard_Integer theIndex)
    {
      Standard_Integer aMoveSteps = 0;
      myCurrControlParam = myControlParams->Value(theIndex);
      myIso->D1(myCurrControlParam, myCurrControlPnt, myCurrControlVec);

      const Standard_Real aMidParam = 0.5 * (myPrevControlParam + myCurrControlParam);
      const gp_Pnt aMidPnt = myIso->Value(aMidParam);

      const Standard_Real aSqDist = BRepMesh_GeomTool::SquareDeflectionOfSegment(
        myPrevControlPnt, myCurrControlPnt, aMidPnt);

      Standard_Real anAngle = 0.0;
      
      if ((myPrevControlVec.SquareMagnitude() > Precision::SquareConfusion()) &&
          (myCurrControlVec.SquareMagnitude() > Precision::SquareConfusion()))
      {
        anAngle = myPrevControlVec.Angle(myCurrControlVec);
      }

      const Standard_Real aSqMaxDeflection = myDFace->GetDeflection() *
        myDFace->GetDeflection();

      if (((aSqDist > aSqMaxDeflection) || (anAngle > myParameters.AngleInterior)) &&
          aSqDist > myParameters.MinSize * myParameters.MinSize)
      {
        // insertion 
        myControlParams->InsertBefore(theIndex, aMidParam);
      }
      else
      {
        // Here we should leave at least 3 parameters as far as
        // we must have at least one parameter related to surface
        // internals in order to prevent movement of triangle body
        // outside the surface in case of highly curved ones, e.g.
        // BSpline springs.
        if (((aSqDist < aSqMaxDeflection) || (anAngle < myParameters.AngleInterior)) &&
            myControlParams->Length() > 3 && theIndex < myControlParams->Length())
        {
          // Remove too dense points
          const Standard_Real aTmpParam = myControlParams->Value(theIndex + 1);
          if (checkParameterForDeflectionAndUpdateCache(aTmpParam))
          {
            ++aMoveSteps;
          }
        }

        myPrevControlParam = myCurrControlParam;
        myPrevControlPnt   = myCurrControlPnt;
        myPrevControlVec   = myCurrControlVec;

        ++aMoveSteps;
      }

      return aMoveSteps;
    }

    //! Checks whether the given param suits specified deflection. Updates cache.
    Standard_Boolean checkParameterForDeflectionAndUpdateCache(const Standard_Real theParam)
    {
      gp_Pnt aTmpPnt;
      gp_Vec aTmpVec;
      myIso->D1(theParam, aTmpPnt, aTmpVec);

      const Standard_Real aTmpMidParam = 0.5 * (myPrevControlParam + theParam);
      const gp_Pnt        aTmpMidPnt = myIso->Value(aTmpMidParam);

      // Lets check next parameter.
      // If it also fits deflection, we can remove previous parameter.
      const Standard_Real aSqDist = BRepMesh_GeomTool::SquareDeflectionOfSegment(
        myPrevControlPnt, aTmpPnt, aTmpMidPnt);

      if (aSqDist < myDFace->GetDeflection() * myDFace->GetDeflection())
      {
        // Lets check parameters for angular deflection.
        if (myPrevControlVec.SquareMagnitude() < gp::Resolution() ||
            aTmpVec.SquareMagnitude()          < gp::Resolution() ||
            myPrevControlVec.Angle(aTmpVec)    < myParameters.AngleInterior)
        {
          // For current Iso line we can remove this parameter.
          myControlParamsToRemove->Add(myCurrControlParam);
          myCurrControlParam = theParam;
          myCurrControlPnt   = aTmpPnt;
          myCurrControlVec   = aTmpVec;
          return Standard_True;
        }
        else
        {
          // We have found a place on the surface refusing 
          // removement of this parameter.
          myParamsForbiddenToRemove       ->Add(myCurrParam);
          myControlParamsForbiddenToRemove->Add(myCurrControlParam);
        }
      }

      return Standard_False;
    }

  private:

    IMeshData::IFaceHandle                myDFace;
    Handle(Geom_Surface)                  mySurface;
    Standard_Boolean                      myIsoU;
    Handle(IMeshData::SequenceOfReal)     myParams;
    Handle(IMeshData::SequenceOfReal)     myControlParams;

    Handle(IMeshData::MapOfReal)          myParamsForbiddenToRemove;
    Handle(IMeshData::MapOfReal)          myControlParamsForbiddenToRemove;

    Handle(NCollection_IncAllocator)      myAllocator;
    Handle(IMeshData::MapOfReal)          myControlParamsToRemove;


    IMeshTools_Parameters                 myParameters;
    NCollection_Handle<GeomAdaptor_Curve> myIso;

    Standard_Real                         myCurrParam;

    Standard_Real                         myCurrControlParam;
    gp_Pnt                                myCurrControlPnt;
    gp_Vec                                myCurrControlVec;

    Standard_Real                         myPrevControlParam;
    gp_Pnt                                myPrevControlPnt;
    gp_Vec                                myPrevControlVec;
  };

  //! Adds param to map if it fits specified range.
  Standard_Boolean addParam(
    const Standard_Real&                           theParam,
    const std::pair<Standard_Real, Standard_Real>& theRange,
    IMeshData::IMapOfReal&                         theParams)
  {
    if (theParam < theRange.first ||
        theParam > theRange.second)
    {
      return Standard_False;
    }

    theParams.Add(theParam);
    return Standard_True;
  }

  //! Initializes parameters map using CN intervals.
  Standard_Boolean initParamsFromIntervals(
    const TColStd_Array1OfReal&                    theIntervals,
    const std::pair<Standard_Real, Standard_Real>& theRange,
    const Standard_Boolean                         isSplitIntervals,
    IMeshData::IMapOfReal&                         theParams)
  {
    Standard_Boolean isAdded = Standard_False;

    for (Standard_Integer i = theIntervals.Lower(); i <= theIntervals.Upper(); ++i)
    {
      const Standard_Real aStartParam = theIntervals.Value(i);
      if (addParam(aStartParam, theRange, theParams))
      {
        isAdded = Standard_True;
      }

      if (isSplitIntervals && i < theIntervals.Upper())
      {
        const Standard_Real aMidParam = (aStartParam + theIntervals.Value(i + 1)) / 2.;
        if (addParam(aMidParam, theRange, theParams))
        {
          isAdded = Standard_True;
        }
      }
    }

    return isAdded;
  }

  //! Checks whether intervals should be split.
  //! Returns true in case if it is impossible to compute normal 
  //! directly on intervals, false is returned elsewhere.
  Standard_Boolean toSplitIntervals (const Handle (Geom_Surface)&  theSurf,
                                     const TColStd_Array1OfReal  (&theIntervals)[2])
  {
    Standard_Integer aIntervalU = theIntervals[0].Lower ();
    for (; aIntervalU <= theIntervals[0].Upper (); ++aIntervalU)
    {
      const Standard_Real aParamU = theIntervals[0].Value(aIntervalU);
      if (Precision::IsInfinite (aParamU))
        continue;

      Standard_Integer aIntervalV = theIntervals[1].Lower ();
      for (; aIntervalV <= theIntervals[1].Upper (); ++aIntervalV)
      {
        gp_Dir aNorm;
        const Standard_Real aParamV = theIntervals[1].Value(aIntervalV);
        if (Precision::IsInfinite (aParamV))
          continue;

        if (GeomLib::NormEstim (theSurf, gp_Pnt2d (aParamU, aParamV), Precision::Confusion (), aNorm) != 0)
        {
          return Standard_True;
        }
        // TODO: do not split intervals if there is no normal in the middle of interval.
      }
    }

    return Standard_False;
  }
}

//=======================================================================
// Function: AdjustRange
// Purpose : 
//=======================================================================
void BRepMesh_NURBSRangeSplitter::AdjustRange()
{
  BRepMesh_DefaultRangeSplitter::AdjustRange();
  mySurfaceType = GetSurface()->GetType();

  if (mySurfaceType == GeomAbs_BezierSurface)
  {
    const std::pair<Standard_Real, Standard_Real>& aRangeU = GetRangeU();
    const std::pair<Standard_Real, Standard_Real>& aRangeV = GetRangeV();

    myIsValid = !(aRangeU.first  < -0.5 ||
                  aRangeU.second >  1.5 ||
                  aRangeV.first  < -0.5 ||
                  aRangeV.second >  1.5);
  }
}

//=======================================================================
// Function: AddPoint
// Purpose :
//=======================================================================
void BRepMesh_NURBSRangeSplitter::AddPoint(const gp_Pnt2d &thePoint)
{
  BRepMesh_DefaultRangeSplitter::AddPoint(thePoint);
  GetParametersU().Add(thePoint.X());
  GetParametersV().Add(thePoint.Y());
}

//=======================================================================
// Function: GenerateSurfaceNodes
// Purpose :
//=======================================================================
Handle(IMeshData::ListOfPnt2d) BRepMesh_NURBSRangeSplitter::GenerateSurfaceNodes(
    const IMeshTools_Parameters &theParameters) const
{
  std::cout << "BRepMesh_NURBSRangeSplitter::GenerateSurfaceNodes" << std::endl;
  if (!initParameters())
  {
    return Handle(IMeshData::ListOfPnt2d)();
  }

  const std::pair<Standard_Real, Standard_Real> &aRangeU = GetRangeU();
  const std::pair<Standard_Real, Standard_Real> &aRangeV = GetRangeV();
  const std::pair<Standard_Real, Standard_Real> &aDeltaOld = GetDelta();

  const Standard_Real aDelta[2] = {
      aDeltaOld.first,
      aDeltaOld.second};

  const Standard_Real aDefFace = GetDFace()->GetDeflection();
  const Handle(BRepAdaptor_Surface) &gFace = GetSurface();
  Handle(Geom_Surface) aSurface = gFace->Surface().Surface();

  const Handle(NCollection_IncAllocator) aTmpAlloc =
      new NCollection_IncAllocator(IMeshData::MEMORY_BLOCK_SIZE_HUGE);

  /*
  //V7
  const Handle(IMeshData::SequenceOfReal) aParams[2] = {
      computeGrainAndFilterParameters(GetParametersU(), gFace->UResolution(aDefFace),
                                      (aRangeU.second - aRangeU.first), aDelta.first, theParameters, aTmpAlloc),

      computeGrainAndFilterParameters(GetParametersV(), gFace->VResolution(aDefFace),
                                      (aRangeV.second - aRang
                                      eV.first), aDelta.second, theParameters, aTmpAlloc)};
  */
  // V0.18 of equivalent of computeGrainAndFilterParameters
  const Standard_Real aRange[2][2] = {
      {aRangeU.second, aRangeU.first},
      {aRangeV.second, aRangeV.first}};
  Standard_Real myMinSize = theParameters.MinSize;
  Handle(IMeshData::SequenceOfReal) aParams[2];
  for (Standard_Integer i = 0; i < 2; ++i)
  {
    Standard_Boolean isU = (i == 0);
    Standard_Real aRes = isU ? gFace->UResolution(aDefFace) : gFace->VResolution(aDefFace);

    // Sort and filter sequence of parameters
    Standard_Real aMinDiff = Precision::PConfusion();
    if (aDelta[i] < 1.)
      aMinDiff /= aDelta[i];

    aMinDiff = Max(myMinSize, aMinDiff);

    Standard_Real aRangeDiff = aRange[i][0] - aRange[i][1];
    Standard_Real aDiffMaxLim = 0.1 * aRangeDiff;
    Standard_Real aDiffMinLim = Max(0.005 * aRangeDiff, 2. * aRes);
    Standard_Real aDiff = Max(myMinSize, Min(aDiffMaxLim, aDiffMinLim));
    aParams[i] = filterParameters(isU ? GetParametersU() : GetParametersV(), aMinDiff, aDiff, aTmpAlloc);
  }
  // End V0.18 of equivalent of computeGrainAndFilterParameters

  /*
  std::cout << "Raw UV params" << std::endl;
  for (int i = 1; i <= GetParametersU().Extent(); i++)
  {
    for (int j = 1; j <= GetParametersV().Extent(); j++)
    {
      std::cout << GetParametersU()(i) << ", " << GetParametersV()(j) << std::endl;
    }
  }
  std::cout << "End Raw UV params" << std::endl;
  */

  Handle(IMeshData::ListOfPnt2d) aNodes = new IMeshData::ListOfPnt2d(
      new NCollection_IncAllocator(IMeshData::MEMORY_BLOCK_SIZE_HUGE));
  /*
  //V7
  // check intermediate isolines
  Handle(IMeshData::MapOfReal) aFixedParams[2] = {
      new IMeshData::MapOfReal(1, aTmpAlloc),
      new IMeshData::MapOfReal(1, aTmpAlloc)};

  const Handle(IMeshData::MapOfReal) aParamsToRemove[2] = {
      AnalyticalFilter(GetDFace(), GeomAbs_IsoV, aParams[1], aParams[0],
                       aFixedParams[1], aFixedParams[0])
          .GetControlParametersToRemove(theParameters),

      AnalyticalFilter(GetDFace(), GeomAbs_IsoU, aParams[0], aParams[1],
                       aFixedParams[0], aFixedParams[1])
          .GetControlParametersToRemove(theParameters),
  };

  aParamsToRemove[0]->Subtract(*aFixedParams[0]);
  aParamsToRemove[1]->Subtract(*aFixedParams[1]);

  // insert nodes of the regular grid

  // insert nodes of the regular grid
  for (Standard_Integer i = 1; i <= aParams[0]->Length(); ++i)
  {
    const Standard_Real aParam1 = aParams[0]->Value(i);
    if (aParamsToRemove[0]->Contains(aParam1))
    {
      // continue;
    }

    for (Standard_Integer j = 1; j <= aParams[1]->Length(); ++j)
    {
      const Standard_Real aParam2 = aParams[1]->Value(j);
      if (aParamsToRemove[1]->Contains(aParam2))
      {
        // continue;
      }
      std::cout << aParam1 << ", " << aParam2 << std::endl;
      aNodes->Append(gp_Pnt2d(aParam1, aParam2));
    }
  }
  //*/
  // check intermediate isolines
  GeomAbs_SurfaceType aSurfType = gFace->GetType();
  Handle_Geom_Surface aBSpline;
  if (aSurfType == GeomAbs_BezierSurface)
  {
    aBSpline = gFace->Bezier();
  }
  else if (aSurfType == GeomAbs_BSplineSurface)
  {
    aBSpline = gFace->BSpline();
  }

  // precision for compare square distances
  const Standard_Real aPrecision = Precision::Confusion();
  const Standard_Real aSqPrecision = aPrecision * aPrecision;
  for (Standard_Integer k = 0; k < 2; ++k)
  {
    Handle(IMeshData::SequenceOfReal) aParams1 = aParams[k];
    Handle(IMeshData::SequenceOfReal) aParams2 = aParams[(k + 1) % 2];
    const Standard_Boolean isU = (k == 0);
    Standard_Integer aStartIndex, aEndIndex;
    if (isU)
    {
      aStartIndex = 1;
      aEndIndex = aParams1->Length();
    }
    else
    {
      aStartIndex = 2;
      aEndIndex = aParams1->Length() - 1;
    }

    for (Standard_Integer i = aStartIndex; i <= aEndIndex; ++i)
    {
      const Standard_Real aParam1 = aParams1->Value(i);
      GeomAdaptor_Curve aIso(isU ? aBSpline->UIso(aParam1) : aBSpline->VIso(aParam1));

      Standard_Real aPrevParam2 = aParams2->Value(1);
      gp_Pnt aPrevPnt2 = aIso.Value(aPrevParam2);
      for (Standard_Integer j = 2; j <= aParams2->Length();)
      {
        Standard_Real aParam2 = aParams2->Value(j);
        gp_Pnt aPnt2 = aIso.Value(aParam2);
        Standard_Real aMidParam = 0.5 * (aPrevParam2 + aParam2);
        gp_Pnt aMidPnt = aIso.Value(aMidParam);

        // 23.03.2010 skl for OCC21645 - change precision for comparison
        Standard_Real aDist;
        if (aPrevPnt2.SquareDistance(aPnt2) > aSqPrecision)
        {
          gp_Lin aLin(aPrevPnt2, gp_Dir(gp_Vec(aPrevPnt2, aPnt2)));
          aDist = aLin.Distance(aMidPnt);
        }
        else
          aDist = aPrevPnt2.Distance(aMidPnt);

        if (aDist > aDefFace && aDist > myMinSize)
        {
          // insertion
          aParams2->InsertBefore(j, aMidParam);
        }
        else
        {
          // put regular grig for normals
          gp_Pnt2d aStPnt1, aStPnt2;
          if (isU)
          {
            aStPnt1 = gp_Pnt2d(aParam1, aPrevParam2);
            aStPnt2 = gp_Pnt2d(aParam1, aMidParam);
          }
          else
          {
            aStPnt1 = gp_Pnt2d(aPrevParam2, aParam1);
            aStPnt2 = gp_Pnt2d(aMidParam, aParam1);
          }

          gp_Dir N1(0, 0, 1), N2(0, 0, 1);
          Standard_Integer aSt1 = GeomLib::NormEstim(aBSpline, aStPnt1, aPrecision, N1);
          Standard_Integer aSt2 = GeomLib::NormEstim(aBSpline, aStPnt2, aPrecision, N2);

          Standard_Real aAngle = N2.Angle(N1);
          if (aSt1 < 1 && aSt2 < 1 && aAngle > .5 * theParameters.Angle)
          {
            Standard_Real aLen = GCPnts_AbscissaPoint::Length(aIso,
                                                              aPrevParam2, aMidParam, aDefFace);

            if (aLen > myMinSize)
            {
              // insertion
              aParams2->InsertBefore(j, aMidParam);
              continue;
            }
          }

          aPrevParam2 = aParam2;
          aPrevPnt2 = aPnt2;

          ++j;
        }
      }
    }
  }

  // insert nodes of the regular grid
  std::cout << "Begin node insertion" << std::endl;
  /*
  for (Standard_Integer i = 1; i <= aParams[0]->Length(); ++i)
  {
    const Standard_Real aParam1 = aParams[0]->Value(i);
    for (Standard_Integer j = 1; j <= aParams[1]->Length(); ++j)
    {
      gp_Pnt2d aPnt2d(aParam1, aParams[1]->Value(j));
      std::cout << 0 << ", " << aPnt2d.X() << " , " << aPnt2d.Y() << std::endl;
      aNodes->Append(aPnt2d);
    }
  }
  //*/
  // std::cout << "Begin raw grid" << std::endl;
  for (int i = 1; i <= GetParametersU().Extent(); i++)
  {
    for (int j = 1; j <= GetParametersV().Extent(); j++)
    {
       aNodes->Append(gp_Pnt2d(GetParametersU()(i), GetParametersV()(j)));
      std::cout << 1 << ", " << GetParametersU()(i) << " , " << GetParametersV()(j) << std::endl;
    }
  }
  std::cout << "End node insertion" << std::endl;
  // std::cout << "End raw grid" << std::endl;
  //*/
  std::cout << "End BRepMesh_NURBSRangeSplitter::GenerateSurfaceNodes" << std::endl;
  return aNodes;
}

//=======================================================================
// Function: initParameters
// Purpose :
//=======================================================================
Standard_Boolean BRepMesh_NURBSRangeSplitter::initParameters() const
{
  const Handle(BRepAdaptor_Surface) &aSurface = GetSurface();

  const GeomAbs_Shape aContinuity = GeomAbs_CN;
  const std::pair<Standard_Integer, Standard_Integer> aIntervalsNb(
      aSurface->NbUIntervals(aContinuity),
      aSurface->NbVIntervals(aContinuity));

  TColStd_Array1OfReal aIntervals[2] = {
      TColStd_Array1OfReal(1, aIntervalsNb.first + 1),
      TColStd_Array1OfReal(1, aIntervalsNb.second + 1)};

  aSurface->UIntervals(aIntervals[0], aContinuity);
  aSurface->VIntervals(aIntervals[1], aContinuity);

  const Standard_Boolean isSplitIntervals = toSplitIntervals(aSurface->Surface().Surface(), aIntervals);

  if (!initParamsFromIntervals(aIntervals[0], GetRangeU(), isSplitIntervals,
                               const_cast<IMeshData::IMapOfReal &>(GetParametersU())))
  {
    // if (!grabParamsOfEdges (Edge_Frontier, Param_U))
    {
      return Standard_False;
    }
  }

  if (!initParamsFromIntervals(aIntervals[1], GetRangeV(), isSplitIntervals,
                               const_cast<IMeshData::IMapOfReal &>(GetParametersV())))
  {
    // if (!grabParamsOfEdges (Edge_Frontier, Param_V))
    {
      return Standard_False;
    }
  }

  return grabParamsOfEdges(Edge_Internal, Param_U | Param_V);
}

//=======================================================================
//function : grabParamsOfInternalEdges
//purpose  : 
//=======================================================================
Standard_Boolean BRepMesh_NURBSRangeSplitter::grabParamsOfEdges (
  const EdgeType         theEdgeType,
  const Standard_Integer theParamDimensionFlag) const
{
  if ((theParamDimensionFlag & (Param_U | Param_V)) == 0)
  {
    return Standard_False;
  }

  const IMeshData::IFaceHandle& aDFace = GetDFace ();
  for (Standard_Integer aWireIt = 0; aWireIt < aDFace->WiresNb (); ++aWireIt)
  {
    const IMeshData::IWireHandle& aDWire = aDFace->GetWire (aWireIt);
    for (Standard_Integer aEdgeIt = 0; aEdgeIt < aDWire->EdgesNb (); ++aEdgeIt)
    {
      const IMeshData::IEdgePtr& aDEdge = aDWire->GetEdge (aEdgeIt);
      for (Standard_Integer aPCurveIt = 0; aPCurveIt < aDEdge->PCurvesNb (); ++aPCurveIt)
      {
        const IMeshData::IPCurveHandle& aDPCurve = aDEdge->GetPCurve (aPCurveIt);
        if (aDPCurve->GetFace () == aDFace)
        {
          if (theEdgeType == Edge_Internal && !aDPCurve->IsInternal ())
          {
            continue;
          }

          for (Standard_Integer aPointIt = 0; aPointIt < aDPCurve->ParametersNb (); ++aPointIt)
          {
            const gp_Pnt2d& aPnt2d = aDPCurve->GetPoint (aPointIt);
            if (theParamDimensionFlag & Param_U)
            {
              const_cast<IMeshData::IMapOfReal&>(GetParametersU ()).Add (aPnt2d.X ());
            }

            if (theParamDimensionFlag & Param_V)
            {
              const_cast<IMeshData::IMapOfReal&>(GetParametersV ()).Add (aPnt2d.Y ());
            }
          }
        }
      }
    }
  }

  return Standard_True;
}

//=======================================================================
//function : computeGrainAndFilterParameters
//purpose  : 
//=======================================================================
Handle(IMeshData::SequenceOfReal) BRepMesh_NURBSRangeSplitter::computeGrainAndFilterParameters(
  const IMeshData::IMapOfReal&            theSourceParams,
  const Standard_Real                     theTol2d,
  const Standard_Real                     theRangeDiff,
  const Standard_Real                     theDelta,
  const IMeshTools_Parameters&            theParameters,
  const Handle(NCollection_IncAllocator)& theAllocator) const
{
  // computeGrainAndFilterParameters(GetParametersU(), gFace->UResolution(aDefFace),
  //         (aRangeU.second - aRangeU.first), aDelta.first, theParameters, aTmpAlloc),

  // Sort and filter sequence of parameters
  Standard_Real aMinDiff = Precision::PConfusion();
  if (theDelta < 1.)
  {
    aMinDiff /= theDelta;
  }

  aMinDiff = Max(theParameters.MinSize, aMinDiff);

  const Standard_Real aDiffMaxLim = 0.1 * theRangeDiff;
  const Standard_Real aDiffMinLim = Max(0.005 * theRangeDiff,
                                        2. * theTol2d);
  const Standard_Real aDiff = Max(theParameters.MinSize,
                                  Min(aDiffMaxLim, aDiffMinLim));
  return filterParameters(theSourceParams, aMinDiff, aDiff, theAllocator);
}

//=======================================================================
//function : filterParameters
//purpose  : 
//=======================================================================
Handle(IMeshData::SequenceOfReal) BRepMesh_NURBSRangeSplitter::filterParameters(
  const IMeshData::IMapOfReal&            theParams,
  const Standard_Real                     theMinDist,
  const Standard_Real                     theFilterDist,
  const Handle(NCollection_IncAllocator)& theAllocator) const
{
  Handle(IMeshData::SequenceOfReal) aResult = new IMeshData::SequenceOfReal(theAllocator);

  // Sort sequence of parameters
  const Standard_Integer anInitLen = theParams.Extent();

  if (anInitLen < 1)
  {
    return aResult;
  }

  TColStd_Array1OfReal aParamArray(1, anInitLen);
  Standard_Integer j;
  for (j = 1; j <= anInitLen; j++)
    aParamArray(j) = theParams(j);

  std::sort(aParamArray.begin(), aParamArray.end());

  Standard_Integer aParamLength = 1;
  /*
  //V7
  // mandatory pre-filtering using the first (minimal) filter value
  for (j = 2; j <= anInitLen; j++)
  {
    if ((aParamArray(j) - aParamArray(aParamLength)) > theMinDist)
    {
      if (++aParamLength < j)
        aParamArray(aParamLength) = aParamArray(j);
    }
  }

  //perform filtering on series
  Standard_Real aLastAdded, aLastCandidate;
  Standard_Boolean isCandidateDefined = Standard_False;
  aLastAdded = aParamArray(1);
  aLastCandidate = aLastAdded;
  aResult->Append(aLastAdded);

  for (j = 2; j < aParamLength; j++)
  {
    Standard_Real aVal = aParamArray(j);
    if (aVal - aLastAdded > theFilterDist)
    {
      //adds the parameter
      if (isCandidateDefined)
      {
        aLastAdded = aLastCandidate;
        isCandidateDefined = Standard_False;
        j--;
      }
      else
      {
        aLastAdded = aVal;
      }
      aResult->Append(aLastAdded);
      continue;
    }

    aLastCandidate = aVal;
    isCandidateDefined = Standard_True;
  }
  aResult->Append(aParamArray(aParamLength));
  */

  // V0.18
  //  mandadory pre-filtering using the first (minimal) filter value
  Handle(IMeshData::SequenceOfReal) aParamTmp = new IMeshData::SequenceOfReal(theAllocator);
  Standard_Real aP1, aP2;
  aP1 = aParamArray(1);
  aParamTmp->Append(aP1);
  for (j = 2; j <= anInitLen; j++)
  {
    aP2 = aParamArray(j);
    if ((aP2 - aP1) > theMinDist)
    {
      aParamTmp->Append(aP2);
      aP1 = aP2;
      aParamLength++;
    }
  }

  // add last point if required
  if (aParamArray(anInitLen) - theParams(aParamLength) > theMinDist)
  {
    aParamTmp->Append(aParamArray(anInitLen));
    aParamLength++;
  }

  // perform filtering on series
  Standard_Real aLastAdded, aLastCandidate;
  Standard_Boolean isCandidateDefined = Standard_False;
  aLastAdded = aParamTmp->First();
  aLastCandidate = aLastAdded;
  aResult->Append(aParamTmp->First());

  for (j = 2; j < aParamTmp->Length(); j++)
  {
    Standard_Real aVal = aParamTmp->Value(j);
    if (aVal - aLastAdded > theFilterDist)
    {
      // adds the parameter
      if (isCandidateDefined)
      {
        aLastAdded = aLastCandidate;
        isCandidateDefined = Standard_False;
        j--;
      }
      else
      {
        aLastAdded = aVal;
      }
      aResult->Append(aLastAdded);
      continue;
    }

    aLastCandidate = aVal;
    isCandidateDefined = Standard_True;
  }
  aResult->Append(aParamTmp->Last());

  if (aResult->Length() == 2)
  {
    Standard_Real dist = aResult->Last() - aResult->First();
    Standard_Integer nbint = (Standard_Integer)((dist / theFilterDist) + 0.5);

    if (nbint > 1)
    {
      // Five points more is maximum
      if (nbint > 5)
      {
        nbint = 5;
      }

      Standard_Integer i;
      Standard_Real dU = dist / nbint;
      for (i = 1; i < nbint; i++)
      {
        aResult->InsertAfter(i, aResult->First() + i * dU);
      }
    }
  }

  // Common to V7 and V0.18
  return aResult;
}
