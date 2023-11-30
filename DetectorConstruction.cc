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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"

#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSPassageCellCurrent.hh"
#include "G4SDParticleWithEnergyFilter.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),
 fPShield(0), fLShield(0), fMaterial(0), fDetectorMessenger(0)
{
  //fBoxSize = 1*m;
  ShThick = 15*cm;
  DefineMaterials();
  SetMaterial("polyethylene");  
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ delete fDetectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // specific element name for thermal neutronHP
  // (see G4ParticleHPThermalScatteringNames.cc)

  G4int ncomponents, natoms;
  G4double massfraction;

  G4double Vdens = 1.e-25*g/cm3;
  G4double Vpres = 1.e-19*pascal;
  G4double Vtemp = 0.1*kelvin;
  
  G4double a, z;

  // vacuum
  G4Material* Vacc = new G4Material("Galactic", z=1, a=1.01*g/mole, Vdens, kStateGas, Vtemp, Vpres);

  // pressurized water
  G4Element* H  = new G4Element("TS_H_of_Water" ,"H" , 1., 1.0079*g/mole);
  G4Element* O  = new G4Element("Oxygen"        ,"O" , 8., 16.00*g/mole);
  G4Material* H2O = 
  new G4Material("Water_ts", 1.000*g/cm3, ncomponents=2,
                         kStateLiquid, 593*kelvin, 150*bar);
  H2O->AddElement(H, natoms=2);
  H2O->AddElement(O, natoms=1);
  H2O->GetIonisation()->SetMeanExcitationEnergy(78.0*eV);
  
  // heavy water
  G4Isotope* H2 = new G4Isotope("H2",1,2);
  G4Element* D  = new G4Element("TS_D_of_Heavy_Water", "D", 1);
  D->AddIsotope(H2, 100*perCent);  
  G4Material* D2O = new G4Material("HeavyWater", 1.11*g/cm3, ncomponents=2,
                        kStateLiquid, 293.15*kelvin, 1*atmosphere);
  D2O->AddElement(D, natoms=2);
  D2O->AddElement(O, natoms=1);
  
  // graphite
  G4Isotope* C12 = new G4Isotope("C12", 6, 12);  
  G4Element* C   = new G4Element("TS_C_of_Graphite","C", ncomponents=1);
  C->AddIsotope(C12, 100.*perCent);
  G4Material* graphite = 
  new G4Material("graphite", 2.27*g/cm3, ncomponents=1,
                         kStateSolid, 293*kelvin, 1*atmosphere);
  graphite->AddElement(C, natoms=1);

  // air
  G4Element* N = new G4Element("Nitrogen", "N", 7., 14.01*g/mole);
  G4Material* Air = new G4Material("air", 1.290*mg/cm3, ncomponents=2, kStateGas, 293*kelvin, 1*atmosphere);
  Air->AddElement(N, massfraction=70.*perCent);
  Air->AddElement(O, massfraction=30.*perCent);

  // iron
  G4Isotope* Fe56 = new G4Isotope("Fe56", 26, 56);
  G4Element* Fe = new G4Element("TS_Iron_Metal", "Fe", ncomponents=1);
  Fe->AddIsotope(Fe56, 100.*perCent);
  G4Material* iron = new G4Material("iron", 7.874*g/cm3, ncomponents=1, kStateSolid, 293*kelvin, 1*atmosphere);
  iron->AddElement(Fe, natoms=1);

  // boron
  G4Isotope* B10 = new G4Isotope("B10", 5, 10);
  G4Isotope* B11 = new G4Isotope("B11", 5, 11);
  G4Element* B = new G4Element("Boron", "B", ncomponents=2);
  B->AddIsotope(B10, 19.9*perCent);
  B->AddIsotope(B11, 80.1*perCent);
  G4Material* boron = new G4Material("boron", 2.46*g/cm3, ncomponents=1, kStateSolid,293*kelvin, 1*atmosphere);
  boron->AddElement(B, natoms=1);

  // polyethilene
  G4Element* Hpe = new G4Element("TS_H_of_Polyethylene", "H", 1, 1.0079*g/mole);
  G4Element* Cpe = new G4Element("Carbon", "C", 6, 12.01*g/mole);
  G4Material* polyethylene = new G4Material("polyethylene", 0.93*g/cm3, ncomponents=2, kStateSolid, 293*kelvin, 1*atmosphere);
  polyethylene->AddElement(Hpe, natoms=4);
  polyethylene->AddElement(Cpe, natoms=2);

  // borated polyethilene
  G4Material* b_polyethylene = new G4Material("b_polyethylene",0.94*g/cm3,ncomponents=4,kStateSolid,293*kelvin,1*atmosphere);
  b_polyethylene->AddElement(Hpe, 11.6*perCent);
  b_polyethylene->AddElement(Cpe, 61.2*perCent);
  b_polyethylene->AddElement(B, 5*perCent);
  b_polyethylene->AddElement(O, 22.2*perCent);
  
 ///G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::MaterialWithSingleIsotope( G4String name,
                           G4String symbol, G4double density, G4int Z, G4int A)
{
 // define a material from an isotope
 //
 G4int ncomponents;
 G4double abundance, massfraction;

 G4Isotope* isotope = new G4Isotope(symbol, Z, A);
 
 G4Element* element  = new G4Element(name, symbol, ncomponents=1);
 element->AddIsotope(isotope, abundance= 100.*perCent);
 
 G4Material* material = new G4Material(name, density, ncomponents=1);
 material->AddElement(element, massfraction=100.*perCent);

 return material;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // Get materials
  auto iron = G4Material::GetMaterial("iron");
  auto polyethylene = G4Material::GetMaterial("polyethylene");
  auto D2O = G4Material::GetMaterial("HeavyWater");
  auto Vacc = G4Material::GetMaterial("Galactic");

  // world
  fBoxSize = 1*m;

  G4Box*
  sBox = new G4Box("world",                             //its name
                   fBoxSize/2,fBoxSize/2,fBoxSize/2);   //its dimensions

  auto fLBox = new G4LogicalVolume(sBox,                     //its shape
                             Vacc,                      //its material
                             "World");                  //its name

  auto fPBox = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(),            //at (0,0,0)
                            fLBox,                      //its logical volume
                            "World",                    //its name
                            0,                          //its mother  volume
                            false,                      //no boolean operation
                            0);                         //copy number

  // shielding
  //ShThick = 15*cm;
  G4double ShSize = 100*cm;
  G4double ShPos = 0*cm;

  G4Box* sShield = new G4Box("shield",
                    ShThick/2,ShSize/2,ShSize/2);

  fLShield = new G4LogicalVolume(sShield,
                                      fMaterial,
                                      "Shield");

  fPShield = new G4PVPlacement(0,
                              G4ThreeVector(ShPos,0.*cm,0.*cm),
                              fLShield,
                              "Shield",
                              fLBox,
                              false,
                              0);

  // scores
  G4double ScThick = 0.5*cm;
  G4double ScSize = 100*cm;
  G4double ScGap = 2.5*cm;
  //G4double ScPos = 45*cm;
  G4double ScPos = (ShThick/2 + ScGap + ScThick/2);

  G4Box* sScore = new G4Box("score",
                            ScThick/2,ScSize/2,ScSize/2);

  auto fLScore = new G4LogicalVolume(sScore,
                                      Vacc,
                                      "Score");

  auto fPScore_r = new G4PVPlacement(0,
                                    G4ThreeVector(ScPos,0.*cm,0.*cm),
                                    fLScore,
                                    "Score_r",
                                    fLBox,
                                    false,
                                    0);

  //auto fPScore_f = new G4PVPlacement(0,
  //                                  G4ThreeVector(-ScPos,0.*cm,0.*cm),
  //                                  fLScore,
  //                                  "Score_f",
  //                                  fLBox,
  //                                  false,
  //                                  0);          

  PrintParameters();
  
  //always return the root volume
  //
  return fPBox;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
  G4cout << "\n The Box is " << G4BestUnit(fBoxSize,"Length")
         << " of " << "Galactic" 
         << "\n \n" << fMaterial << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  
  if (pttoMaterial) { 
    if(fMaterial != pttoMaterial) {
      fMaterial = pttoMaterial;
      if(fLShield) { fLShield->SetMaterial(pttoMaterial); }
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetMaterial : "
           << materialChoice << " not found" << G4endl;
  }              
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSize(G4double value)
{
  ShThick = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
  auto Score_r = new G4MultiFunctionalDetector("score_r");

  G4int idx = 0;

  for (G4int i=-8; i<1; i++) {
    for (G4int j=1; j<10; j++) {

    //generate names
    G4String Fname = "F_" + std::to_string(i+8) + "_" + std::to_string(j);
    G4String Sname = "re" + std::to_string(i+8) + "_" + std::to_string(j);
    
    //define filters
    G4double Felow = j*pow(10,i);
    G4double Fehigh = (j+1)*pow(10,i);
    
    G4SDParticleWithEnergyFilter* Filt;
    Filt = new G4SDParticleWithEnergyFilter(Fname,Felow*MeV,Fehigh*MeV);
    Filt->add("neutron");

    //define scores (current) and assign filters
    G4VPrimitiveScorer* re;
    re = new G4PSPassageCellCurrent(Sname);
    re->SetFilter(Filt);
    Score_r->RegisterPrimitive(re);

    }
  }

  //define filters
  //auto* F0 = new G4SDParticleWithEnergyFilter("F0",0.*eV,1.*eV);
  //F0->add("neutron");
  //auto* F1 = new G4SDParticleWithEnergyFilter("F1",1*eV,500.*eV);
  //F1->add("neutron");
  //auto* F2 = new G4SDParticleWithEnergyFilter("F2",0.5*keV,1.*keV);
  //F2->add("neutron");
  //auto* F3 = new G4SDParticleWithEnergyFilter("F3",1*keV,500.*keV);
  //F3->add("neutron");
  //auto* F4 = new G4SDParticleWithEnergyFilter("F4",0.5*MeV,1.*MeV);
  //F4->add("neutron");
  //auto* F5 = new G4SDParticleWithEnergyFilter("F5",1.*MeV,2.*MeV);
  //F5->add("neutron");
  //auto* FT = new G4SDParticleFilter("FT");
  //FT->add("neutron");

  //define scores (current) and asign filters
  //auto Score_r = new G4MultiFunctionalDetector("score_r");
  //G4VPrimitiveScorer* re0 = new G4PSPassageCellCurrent("re0");
  //re0->SetFilter(F0);
  //Score_r->RegisterPrimitive(re0);
  //G4VPrimitiveScorer* re1 = new G4PSPassageCellCurrent("re1");
  //re1->SetFilter(F1);
  //Score_r->RegisterPrimitive(re1);
  //G4VPrimitiveScorer* re2 = new G4PSPassageCellCurrent("re2");
  //re2->SetFilter(F2);
  //Score_r->RegisterPrimitive(re2);
  //G4VPrimitiveScorer* re3 = new G4PSPassageCellCurrent("re3");
  //re3->SetFilter(F3);
  //Score_r->RegisterPrimitive(re3);
  //G4VPrimitiveScorer* re4 = new G4PSPassageCellCurrent("re4");
  //re4->SetFilter(F4);
  //Score_r->RegisterPrimitive(re4);
  //G4VPrimitiveScorer* re5 = new G4PSPassageCellCurrent("re5");
  //re5->SetFilter(F5);
  //Score_r->RegisterPrimitive(re5);
  //G4VPrimitiveScorer* reT = new G4PSPassageCellCurrent("reT");
  //reT->SetFilter(FT);
  //Score_r->RegisterPrimitive(reT);

  G4SDManager::GetSDMpointer()->AddNewDetector(Score_r);
  SetSensitiveDetector("Score",Score_r);
}
