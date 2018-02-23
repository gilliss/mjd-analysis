//---------------------------------------------------------------------------//
//bb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nu//
//                                                                           //
//                         MAJORANA Simulation                               //
//                                                                           //
//      This code implementation is the intellectual property of the         //
//      MAJORANA Collaboration. It is based on Geant4, an intellectual       //
//      property of the RD44 GEANT4 collaboration.                           //
//                                                                           //
//                        *********************                              //
//                                                                           //
//    Neither the authors of this software system, nor their employing       //
//    institutes, nor the agencies providing financial support for this      //
//    work  make  any representation or  warranty, express or implied,       //
//    regarding this software system or assume any liability for its use.    //
//    By copying, distributing or modifying the Program (or any work based   //
//    on on the Program) you indicate your acceptance of this statement,     //
//    and all its terms.                                                     //
//                                                                           //
//bb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nu//
//---------------------------------------------------------------------------//
/**                                                            
 *      
 * CLASS DECLARATION:  MGStartTimeFinder.hh
 *
 * DESCRIPTION: 
 *
 * A class estimating the start time (t0) of a waveform. Based off of MGDCRSlopeCalculator and that class' interaction with GATWFCalcParamProcessor in process_mjd_data.
 *
 * AUTHOR: T. Gilliss
 * CONTACT: gilliss@unc.edu
 * FIRST SUBMISSION: 3/31/2016
 */

#ifndef _MGStartTimeFinder_HH
#define _MGStartTimeFinder_HH

#include "MGVWaveformParameterCalculator.hh"

class MGStartTimeFinder : public MGVWaveformParameterCalculator
{
  public:
    MGStartTimeFinder(const std::string& aName= "MGStartTimeFinder");
    virtual ~MGStartTimeFinder() {}
    virtual void CalculateParameters(const MGWaveform& wf);
};

#endif
