//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: B4cEventAction.hh 75215 2013-10-29 16:07:06Z gcosmo $
// 
/// \file B4cEventAction.hh
/// \brief Definition of the B4cEventAction class

#ifndef B4cEventAction_h
#define B4cEventAction_h 1

#define TIMEBINS 1024

#include <fstream>

#include "G4UserEventAction.hh"
#include "B4cCalorHit.hh"
#include "globals.hh"
//ROOT Class
#include "TH1F.h"
#include "TTimeStamp.h"

/// Event action class
///
/// In EndOfEventAction(), it prints the accumulated quantities of the energy 
/// deposit and track lengths of charged particles in Absober and Gap layers 
/// stored in the hits collections.

typedef signed short SSHORT;
typedef unsigned short USHORT;



class B4cEventAction : public G4UserEventAction
{
public:
	double mpe; // the average energy deposited in the octagon scintilator from muon of 30 GeV
    G4double* fEdepTimeShape; 
    SSHORT* fDetectorTimeShape; 
    //TH1F *UnitResponseShape;

  B4cEventAction();
  virtual ~B4cEventAction();

  virtual void  BeginOfEventAction(const G4Event* event);
  virtual void    EndOfEventAction(const G4Event* event);

  void AddEdep(G4double edep) { fEdep += edep; }
  void AddTimeShape (G4double time, G4double E){
	  //int timebins(TIMEBINS);
    if (time<200.) fEdepTimeShape[G4int(time*100./512.)] += E;/// Time shape bins of DRS4, where binwidth ~ 200 ps
															/// The bin width of DRS4 vary roghly from 100 - 300 ps 
    else G4cout<<" time too LARGE "<<time;
  }


private:
  // methods
  B4cCalorHitsCollection* GetHitsCollection(G4int hcID,
                                            const G4Event* event) const;
  void PrintEventStatistics(G4double absoEdep, G4double absoTrackLength,
                            G4double gapEdep, G4double gapTrackLength) const;
  
  // data members                   
  G4int  fAbsHCID;
  G4int  fGapHCID;
  G4double     fEdep;
  
	unsigned int *year, *month, *day, *hour, *min, *sec;
	USHORT yU, moU, dU, hU, miU, sU;
	TTimeStamp ts;


  //B4RunAction* fRunAction;

};
                     
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
