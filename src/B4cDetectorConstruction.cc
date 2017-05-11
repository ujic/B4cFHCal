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
// $Id: B4cDetectorConstruction.cc 87359 2014-12-01 16:04:27Z gcosmo $
// 
/// \file B4cDetectorConstruction.cc
/// \brief Implementation of the B4cDetectorConstruction class

#include "B4cDetectorConstruction.hh"
#include "B4cCalorimeterSD.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4RotationMatrix.hh"
#include "G4SubtractionSolid.hh"
#include "G4ThreeVector.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4UserLimits.hh" // pedja

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal 
G4GlobalMagFieldMessenger* B4cDetectorConstruction::fMagFieldMessenger = 0; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cDetectorConstruction::B4cDetectorConstruction()
 : G4VUserDetectorConstruction(),
   fCheckOverlaps(true),
   fNofLayers(-1)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cDetectorConstruction::~B4cDetectorConstruction()
{ 
	
	
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4cDetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();
  
  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cDetectorConstruction::DefineMaterials()
{ 
    G4double a;  // mass of a mole;
    G4double z;  // z=mean number of protons;
    G4double density;

  // Iron material defined using NIST Manager
  G4NistManager* nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Fe");


  a = 1.01*g/mole;
  G4Element* elH  = new G4Element("Hydrogen","H" , z= 1., a);
  a = 12.01*g/mole;
  G4Element* elC  = new G4Element("Carbon"  ,"C" , z= 6., a);

  density = 1.032*g/cm3;
  G4Material* Sci = new G4Material("Scintillator", density, 2); // number of components = 2
  Sci->AddElement(elC, 9); // natoms=9
  Sci->AddElement(elH, 10); // natoms=10

  /*
  new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
         // The argon by NIST Manager is a gas with a different density
*/
  // Vacuum
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4cDetectorConstruction::DefineVolumes()
{
  // Geometry parameters
  fNofLayers = 1;

  G4double absoThickness = 300.*mm;
  G4double gapThickness =  5.*mm;
  G4double calorSizeXY  = 30.*cm;
  G4double gapSizeXY = 8.5*cm; // the real size is 8.5 cm and the corners are cut (to see how much) ...
  G4double gap2abs  = 10*cm;
  G4double cornerDistance = 4.25*cm ; // at the distance (sqrt(2)+1)/2*gapSizeXY is the touching distance

  G4double layerThickness = 2*absoThickness + gapThickness+2*gap2abs;
  G4double calorThickness = fNofLayers * layerThickness;
  G4double worldSizeXY = 1.2 * calorSizeXY;
  G4double worldSizeZ  = 1.2 * calorThickness; 
  
  // Get materials
  G4Material* defaultMaterial = G4Material::GetMaterial("Galactic");
  G4Material* absorberMaterial = G4Material::GetMaterial("G4_Fe");
  G4Material* gapMaterial = G4Material::GetMaterial("Scintillator");
  
  if ( ! defaultMaterial || ! absorberMaterial || ! gapMaterial ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined."; 
    G4Exception("B4DetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }  
   
  //     
  // World
  //
  G4VSolid* worldS 
    = new G4Box("World",           // its name
                 worldSizeXY/2, worldSizeXY/2, worldSizeZ/2); // its size
                         
  G4LogicalVolume* worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 defaultMaterial,  // its material
                 "World");         // its name
                                   
  G4VPhysicalVolume* worldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume                         
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  //                               
  // Calorimeter
  //  
  G4VSolid* calorimeterS
    = new G4Box("Calorimeter",     // its name
                 calorSizeXY/2, calorSizeXY/2, calorThickness/2); // its size
                         
  G4LogicalVolume* calorLV
    = new G4LogicalVolume(
                 calorimeterS,     // its solid
                 defaultMaterial,  // its material
                 "Calorimeter");   // its name
                                   
  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 calorLV,          // its logical volume                         
                 "Calorimeter",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
   
  /*
  //                                 
  // Layer
  //
  G4VSolid* layerS 
    = new G4Box("Layer",           // its name
                 calorSizeXY/2, calorSizeXY/2, layerThickness/2); //its size
                         
  G4LogicalVolume* layerLV
    = new G4LogicalVolume(
                 layerS,           // its solid
                 defaultMaterial,  // its material
                 "Layer");         // its name

  new G4PVReplica(
                 "Layer",          // its name
                 layerLV,          // its logical volume
                 calorLV,          // its mother
                 kZAxis,           // axis of replication
                 fNofLayers,        // number of replica
                 layerThickness);  // witdth of replica
*/
  //                               
  // Absorber
  //
  G4VSolid* absorberS 
    = new G4Box("Abso",            // its name
                 calorSizeXY/2, calorSizeXY/2, absoThickness/2); // its size
                         
  G4LogicalVolume* absorberLV
    = new G4LogicalVolume(
                 absorberS,        // its solid
                 absorberMaterial, // its material
                 "AbsoLV");        // its name
                                   
   new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., -gapThickness/2-gap2abs-absoThickness/2), // its position
                 absorberLV,       // its logical volume                         
                 "Abso",           // its name
                 calorLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  //
  // Absorber2
  //
  G4VSolid* absorber2S
    = new G4Box("Abso2",            // its name
                 calorSizeXY/2, calorSizeXY/2, absoThickness/2); // its size

  G4LogicalVolume* absorber2LV
    = new G4LogicalVolume(
                 absorber2S,        // its solid
                 absorberMaterial, // its material
                 "Abso2LV");          // its name

  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., gapThickness/2+gap2abs+absoThickness/2), // its position
                 absorber2LV,       // its logical volume
                 "Abso2",           // its name
                 calorLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps

  /// MORA SE DODATI ABSORBER 2 JER JERE SE IZ G4PVPLACEMENTA UZIMA E I TRACK ZA HISTOGRAME ...


  //                               
  // Gap
  //
 /* G4VSolid* gapS 
    = new G4Box("Gap",             // its name
                 gapSizeXY/2, gapSizeXY/2, gapThickness/2); // its size
   */ 
   
  G4VSolid* box 
    = new G4Box("box",             // its name
                 gapSizeXY/2, gapSizeXY/2, gapThickness/2); // its size              
  G4VSolid* CornerCut = new G4Box("ConrnerCut", 1.5*cm, 4*cm, gapThickness);
  G4RotationMatrix* rm= new G4RotationMatrix();
	rm->rotateZ(-45.*deg);
  G4VSolid* subtract1 = new G4SubtractionSolid("subtract1", box, CornerCut, rm, G4ThreeVector(cornerDistance,cornerDistance,0.));
  rm->rotateZ(90.*deg);
  G4VSolid* subtract2 = new G4SubtractionSolid("subtract2", subtract1, CornerCut, rm, G4ThreeVector(cornerDistance,-cornerDistance,0.));
  rm->rotateZ(90.*deg);
  G4VSolid* subtract3 = new G4SubtractionSolid("subtract3", subtract2, CornerCut, rm, G4ThreeVector(-cornerDistance,-cornerDistance,0.));  
  rm->rotateZ(90.*deg);
  G4VSolid* gapS	  = new G4SubtractionSolid("Gap", subtract3, CornerCut, rm, G4ThreeVector(-cornerDistance,cornerDistance,0.));
                
  G4LogicalVolume* gapLV = new G4LogicalVolume(
                 gapS,             // its solid
                 gapMaterial,      // its material
                 "GapLV");         // its name
                 

                                   
  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(10.*cm, 0., 0.), // its position
                 gapLV,            // its logical volume                         
                 "Gap",            // its name
                 calorLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

    fScoringVolume = gapLV; //pedja add as in B1
 
  //
  // print parameters
  //
  G4cout
    << G4endl 
    << "------------------------------------------------------------" << G4endl
    << "---> The calorimeter is " << fNofLayers << " layers of: [ "
    << absoThickness/mm << "mm of " << absorberMaterial->GetName() 
    << " + "
    << gapThickness/mm << "mm of " << gapMaterial->GetName() << " ] " << G4endl
    << "------------------------------------------------------------" << G4endl;
  
  //                                        
  // Visualization attributes
  //
  worldLV->SetVisAttributes (G4VisAttributes::Invisible);

  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  calorLV->SetVisAttributes(simpleBoxVisAtt);

  //  User limits
  
  G4double maxStep = 0.2*gapThickness;				//pedja
  //G4double maxLength = 0.2*gapThickness;  
  //G4double maxTime = 0.1*ns;
  //G4double minEkin = 0.05*MeV;   // 
  gapLV->SetUserLimits(new G4UserLimits(maxStep));//,maxLength,maxTime, minEkin, 0.));	//pedja


  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cDetectorConstruction::ConstructSDandField()
{
  // G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  // 
  // Sensitive detectors
  //
  B4cCalorimeterSD* absoSD 
    = new B4cCalorimeterSD("AbsorberSD", "AbsorberHitsCollection", fNofLayers);
  SetSensitiveDetector("AbsoLV",absoSD);

  B4cCalorimeterSD* gapSD 
    = new B4cCalorimeterSD("GapSD", "GapHitsCollection", fNofLayers);
  SetSensitiveDetector("GapLV",gapSD);

  // 
  // Magnetic field
  //
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue = G4ThreeVector();
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);
  
  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
