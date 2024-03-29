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
// --------------------------------------------------------------
//

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1
#endif

#include "G4VUserDetectorConstruction.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RunManager.hh"
#include "DetectorConstructionMessenger.hh"

#include "globals.hh"


#include "G4UserLimits.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//class G4UserLimits;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    
    DetectorConstruction();
    DetectorConstruction(G4bool);
    ~DetectorConstruction();
    
    void DefineMaterials();
    G4VPhysicalVolume* Construct();

private:
    void ConstructSDandField();

private:
    DetectorConstructionMessenger* fMessenger;
    G4VPhysicalVolume* physAbsor;
    G4VPhysicalVolume* worldPhysical; 
    // addition
    G4UserLimits* stepLimit;             // pointer to user step limits

public:
    G4String GetEC() {return fECfileName;}
    void SetEC(G4String aString) {fECfileName = aString;}

    G4String GetECO() {return fECOfileName;}
    void SetECO(G4String aString) {fECOfileName = aString;}
    
    G4String GetMaterial() {return fMaterialName;}
    void SetMaterial(G4String aString) {fMaterialName = aString;}
    
    G4ThreeVector GetSizes() {return fSizes;}
    void SetSizes(G4ThreeVector a3vec) {fSizes = a3vec;}
    
    G4ThreeVector GetBR() {return fBR;}
    void SetBR(G4ThreeVector a3vec) {fBR = a3vec;}

    G4String GetBRFileName() {return fBRFileName;}
    void SetBRFileName(G4String aString) {fBRFileName = aString;}

    G4ThreeVector GetAngles() {return fAngles;}
    void SetAngles(G4ThreeVector a3vec) {fAngles = a3vec;}
	//mano sudas
    void SetAmorph(G4int i) 	{amorfinis = i;}
    G4int GetAmorp() {return amorfinis;}

// mano

     G4VPhysicalVolume* GetAbsorber()   {return physAbsor;};
     const G4VPhysicalVolume* GetWorld()   {return worldPhysical;};
     G4Material* GetMaterialM() 	{return mat;};
     G4double GetLength()		{return fSizes.z();};


private:
    G4String fECfileName;
    G4String fECOfileName;
    G4String fMaterialName;
    G4ThreeVector fSizes;
    G4ThreeVector fBR;
    G4String fBRFileName;
    G4ThreeVector fAngles;
    G4int amorfinis;


private:
    G4String fWorldMaterial;
    G4String fDetectorMaterialName;
    G4ThreeVector fDetectorSizes;
    G4double fDetectorDistance[5];

    //mano
    G4Material* mat;

public:
    void SetWorldMaterial(G4String aString) {fWorldMaterial = aString;}
    G4String GetWorldMaterial() {return fWorldMaterial;}

    void SetDetectorMaterial(G4String aString) {fDetectorMaterialName = aString;}
    G4String GetDetectorMaterial() {return fDetectorMaterialName;}
    
    void SetDetectorSizes(G4ThreeVector a3vec) {fDetectorSizes = a3vec;}
    G4ThreeVector GetDetectorSizes() {return fDetectorSizes;}

    void SetDetectorDistance(G4int aInt,G4double aDouble) {fDetectorDistance[aInt] = aDouble;}
    G4double GetDetectorDistance(G4int aInt) {return fDetectorDistance[aInt];}

    //pridetas amorfizavimas, testui
    void SetAmorphous(bool a) {fCrystalAmorphous = a;}
    
private:
    G4bool fCrystalAmorphous;
};


