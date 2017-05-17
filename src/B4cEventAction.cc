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
// $Id: B4cEventAction.cc 88427 2015-02-19 08:19:38Z gcosmo $
// 
/// \file B4cEventAction.cc
/// \brief Implementation of the B4cEventAction class



#include "B4cEventAction.hh"
#include "B4cCalorimeterSD.hh"
#include "B4cCalorHit.hh"
#include "B4Analysis.hh"
/// ******
#include "GetUnitHistogram.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

#define FREE_ARG char*

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cEventAction::B4cEventAction()
 : G4UserEventAction(),
   fAbsHCID(-1),
   fGapHCID(-1),
   mpe(1.1)
{
    fEdepTimeShape = new G4double[TIMEBINS]; //(G4double *)malloc((size_t) (TIMEBINS*sizeof(G4double)));
    if (!fEdepTimeShape) G4cout<<"***************** !!!!!!!!!!!!!!  allocation failure in vector()  ***************!!!!!!!!!!!!!!!!!!!!!!"<<G4endl;
	fDetectorTimeShape = new SSHORT[TIMEBINS];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cEventAction::~B4cEventAction()
{

    
delete fEdepTimeShape;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cCalorHitsCollection* 
B4cEventAction::GetHitsCollection(G4int hcID,
                                  const G4Event* event) const
{
  B4cCalorHitsCollection* hitsCollection 
    = static_cast<B4cCalorHitsCollection*>(
        event->GetHCofThisEvent()->GetHC(hcID));
  
  if ( ! hitsCollection ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << hcID; 
    G4Exception("B4cEventAction::GetHitsCollection()",
      "MyCode0003", FatalException, msg);
  }         

  return hitsCollection;
}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cEventAction::PrintEventStatistics(
                              G4double absoEdep, G4double absoTrackLength,
                              G4double gapEdep, G4double gapTrackLength) const
{
  // print event statistics
  G4cout
     << "   Absorber: total energy: " 
     << std::setw(7) << G4BestUnit(absoEdep, "Energy")
     << "       total track length: " 
     << std::setw(7) << G4BestUnit(absoTrackLength, "Length")
     << G4endl
     << "        Gap: total energy: " 
     << std::setw(7) << G4BestUnit(gapEdep, "Energy")
     << "       total track length: " 
     << std::setw(7) << G4BestUnit(gapTrackLength, "Length")
     << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cEventAction::BeginOfEventAction(const G4Event* /*event*/)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cEventAction::EndOfEventAction(const G4Event* event)
{  
	//const float ZeroFloat(0.);
	const USHORT ZeroUSHORT(0);
	const SSHORT ZeroSSHORT(0);
	const unsigned int ZeroUInt(0);
	unsigned int UeventID;
	RootHistograms *InstHistograms = new RootHistograms();
	TH1F* UnitResponseShape=InstHistograms->MakeUnitHistogram();
	int timebins(TIMEBINS);
	
	year = new unsigned int;
	month = new unsigned int;
	day = new unsigned int;
	hour = new unsigned int;
	min = new unsigned int;
	sec = new unsigned int;	
	
	std::ofstream OutDRS4FileEvents;// file will be open once by B4RunAction::BeginOfRunAction to write headers, then closed
							// then it will be open for each event by B4cEventAction::EndOfEventAction
							// to APPEND the event, and then close, just to be reopen for the next event... and so on
	
	OutDRS4FileEvents.open("DRS4.dat", std::ofstream::binary | std::ofstream::app);
  // Get hits collections IDs (only once)
  if ( fAbsHCID == -1 ) {
    fAbsHCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("AbsorberHitsCollection");
    fGapHCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("GapHitsCollection");
  }

  // Get hits collections
  B4cCalorHitsCollection* absoHC = GetHitsCollection(fAbsHCID, event);
  B4cCalorHitsCollection* gapHC = GetHitsCollection(fGapHCID, event);

  // Get hit with total values
  B4cCalorHit* absoHit = (*absoHC)[absoHC->entries()-1];
  B4cCalorHit* gapHit = (*gapHC)[gapHC->entries()-1];
 
  // Print per event (modulo n)
  //
  G4int eventID = event->GetEventID();
  G4int printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
    G4cout << "---> End of event: " << eventID << G4endl;     

    PrintEventStatistics(
      absoHit->GetEdep(), absoHit->GetTrackLength(),
      gapHit->GetEdep(), gapHit->GetTrackLength());
  }  
  
//fRunAction->AddEdep(fEdep);   ///pedja
//fRunAction->TimeShapeTotal(fEdepTimeShape); /// pedja

  // Fill histograms, ntuple
  //

  // get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
 
  // fill histograms
  analysisManager->FillH1(1, absoHit->GetEdep());
  analysisManager->FillH1(2, gapHit->GetEdep());
  analysisManager->FillH1(3, absoHit->GetTrackLength());
  analysisManager->FillH1(4, gapHit->GetTrackLength());

  //G4cout<<"fEdepTimeShape: "<<fEdepTimeShape<<G4endl;

  //G4cout<<"fEdepTimeShape[1]: "<<fEdepTimeShape[1]<<G4endl;


  for(int i = 0; i<TIMEBINS; i++) {
	  // i correspond to one bin in DRS4
    analysisManager->FillH1(5, G4double(i), fEdepTimeShape[i]);
    
    if(fEdepTimeShape[i]!=0) for(int j=176; j<396; j++) 
			// the Unit signal starts at bin 176 (end 396) in UnitResponseShape, total width  220bins (44 ns)
			// mpe is the most probable energy of Unit response (muon)
			// delay on signal, since DRS4Analysis calculates the base line using first 35 ns of the signal
			if ((i+j)<999) fDetectorTimeShape[i+j+25]+=static_cast<SSHORT>(-10*fEdepTimeShape[i]/mpe*UnitResponseShape->GetBinContent(j));
	fEdepTimeShape[i]=0;

	
  }
	UeventID = static_cast<unsigned int>(eventID);
	OutDRS4FileEvents<<"EHDR";
	OutDRS4FileEvents.write(reinterpret_cast<const char *>(&UeventID), sizeof(UeventID));
	ts.GetDate(1,0,year, month, day);
	ts.GetTime(1,0,hour, min, sec);
	yU = static_cast<USHORT>(*year);
	moU= static_cast<USHORT>(*month);
	dU= static_cast<USHORT>(*day);
	hU= static_cast<USHORT>(*hour);
	miU= static_cast<USHORT>(*min);
	sU= static_cast<USHORT>(*sec);
	OutDRS4FileEvents.write(reinterpret_cast<const char *>(&yU), sizeof(yU));
	OutDRS4FileEvents.write(reinterpret_cast<const char *>(&moU), sizeof(moU));
	OutDRS4FileEvents.write(reinterpret_cast<const char *>(&dU), sizeof(dU));
	OutDRS4FileEvents.write(reinterpret_cast<const char *>(&hU), sizeof(hU));
	OutDRS4FileEvents.write(reinterpret_cast<const char *>(&miU), sizeof(miU));
	OutDRS4FileEvents.write(reinterpret_cast<const char *>(&sU), sizeof(sU));
	OutDRS4FileEvents.write(reinterpret_cast<const char *>(&ZeroUSHORT), sizeof(ZeroUSHORT));//miliseconds
	OutDRS4FileEvents.write(reinterpret_cast<const char *>(&ZeroUSHORT), sizeof(ZeroUSHORT));//range
	
	OutDRS4FileEvents<<"B#";
	OutDRS4FileEvents.write(reinterpret_cast<const char *>(&ZeroUSHORT), sizeof(ZeroUSHORT));
	OutDRS4FileEvents<<"T#";
	OutDRS4FileEvents.write(reinterpret_cast<const char *>(&ZeroUSHORT), sizeof(ZeroUSHORT));
	OutDRS4FileEvents<<"C001";	
	OutDRS4FileEvents.write(reinterpret_cast<const char *>(&ZeroUInt), sizeof(ZeroUInt));// scaller
	OutDRS4FileEvents.write(reinterpret_cast<const char*>(fDetectorTimeShape), 2*timebins);// /sizeof(fDetectorTimeShape[0]));
    
    
    //write same for ch2
    OutDRS4FileEvents<<"C002";
	OutDRS4FileEvents.write(reinterpret_cast<const char *>(&ZeroUInt), sizeof(ZeroUInt));// scaller
    OutDRS4FileEvents.write(reinterpret_cast<const char*>(fDetectorTimeShape), 2*timebins);
   
	// write ch3 and ch4, used in a coincidence as the trigger
	// the advance of 10bins (2ns) introduce for ch1 and ch2
	// maybe the amplitude of ch3,4 should be augmented?
    OutDRS4FileEvents<<"C003";
	OutDRS4FileEvents.write(reinterpret_cast<const char *>(&ZeroUInt), sizeof(ZeroUInt));// scaller
	OutDRS4FileEvents.write(reinterpret_cast<const char*>(fDetectorTimeShape+10),  2*(timebins - 10));
	OutDRS4FileEvents.write(reinterpret_cast<const char *>(&ZeroSSHORT), sizeof(ZeroSSHORT)*10);
    OutDRS4FileEvents<<"C004";
	OutDRS4FileEvents.write(reinterpret_cast<const char *>(&ZeroUInt), sizeof(ZeroUInt));// scaller
	OutDRS4FileEvents.write(reinterpret_cast<const char*>(fDetectorTimeShape+10), 2*(timebins - 10));
	OutDRS4FileEvents.write(reinterpret_cast<const char *>(&ZeroSSHORT), sizeof(ZeroSSHORT)*10);

OutDRS4FileEvents.close();
  UnitResponseShape->Delete(); 
  //delete UnitResponseShape;
  delete InstHistograms;
  
  delete year;
  delete month;
  delete day;
  delete hour;
  delete min;
  delete sec;
  
  //for (i=0; i<1024; i++) fDetectorTimeShape[i]=0;
  memset(fDetectorTimeShape, 0, timebins * sizeof(fDetectorTimeShape[0]));

  // fill ntuple
  analysisManager->FillNtupleDColumn(0, absoHit->GetEdep());
  analysisManager->FillNtupleDColumn(1, gapHit->GetEdep());
  analysisManager->FillNtupleDColumn(2, absoHit->GetTrackLength());
  analysisManager->FillNtupleDColumn(3, gapHit->GetTrackLength());
  analysisManager->AddNtupleRow(); 

}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
