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
//

#include "DetectorConstruction.hh"
#include "G4RunManager.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4CrystalExtension.hh"
#include "G4ExtendedMaterial.hh"
#include "G4LogicalCrystalVolume.hh"

#include "G4ChannelingMaterialData.hh"
#include "G4ChannelingOptrMultiParticleChangeCrossSection.hh"

#include "CrystalDetector.hh"
#include "SensitiveDetector.hh"

#include "G4SDManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction():
fECfileName("data/Si220pl"),stepLimit(0),
fECOfileName(""),
fMaterialName("G4_Si"),
fSizes(G4ThreeVector(1.,1.,1.)),
fBR(G4ThreeVector(0.,0.,0.)),
fBRFileName(""),
fAngles(G4ThreeVector(0.,0.,0.)),
fWorldMaterial("G4_Galactic"),
fDetectorMaterialName(""),
fDetectorSizes(G4ThreeVector(50. * CLHEP::mm,50. * CLHEP::mm,1 * CLHEP::mm)),
fDetectorDistance{-20. * CLHEP::cm,-19. * CLHEP::cm,+19. * CLHEP::cm,+20. * CLHEP::cm,+40. * CLHEP::cm},
fCrystalAmorphous(false){
    fMessenger = new DetectorConstructionMessenger(this);  //adition
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction(G4bool aBool):
fECfileName("data/Si220pl"),stepLimit(0),
fECOfileName(""),
fMaterialName("G4_Si"),
fSizes(G4ThreeVector(1.,1.,1.)),
fBR(G4ThreeVector(0.,0.,0.)),
fBRFileName(""),
fAngles(G4ThreeVector(0.,0.,0.)),
fWorldMaterial("G4_Galactic"),
fDetectorMaterialName(""),
fDetectorSizes(G4ThreeVector(50. * CLHEP::mm,50. * CLHEP::mm,1 * CLHEP::mm)),
fDetectorDistance{-20. * CLHEP::cm,-19. * CLHEP::cm,+19. * CLHEP::cm,+20. * CLHEP::cm,+40. * CLHEP::cm},
fCrystalAmorphous(aBool){
    fMessenger = new DetectorConstructionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

DetectorConstruction::~DetectorConstruction()
{
  //delete stepLimit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void DetectorConstruction::DefineMaterials(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct(){
    
    //** World **//
    G4Material* worldMaterial = G4NistManager::Instance()->FindOrBuildMaterial(fWorldMaterial);
    
    G4double worldSizeXY = 2. * CLHEP::meter;
    G4double worldSizeZ = 10. * CLHEP::meter; // buvo 30 metru
    
    G4Box* worldSolid = new G4Box("world.solid",
                                  worldSizeXY/2.,
                                  worldSizeXY/2.,
                                  worldSizeZ/2.);
    
    G4LogicalVolume* worldLogic = new G4LogicalVolume(worldSolid,
                                                      worldMaterial,
                                                      "world.logic");
    
    /*G4PVPlacement* */worldPhysical = new G4PVPlacement(0,
                                                     G4ThreeVector(),
                                                     worldLogic,
                                                     "world.physic",
                                                     0,
                                                     false,
                                                     0);
    
    
    //** Detectors instantiation **//
    G4Box* ssdSolid = new G4Box("ssd.solid",
                                fDetectorSizes.x()/2.,
                                fDetectorSizes.y()/2.,
                                fDetectorSizes.z()/2.);
    
    G4Material* detectorMaterial;
    if(fDetectorMaterialName == ""){
        detectorMaterial = worldMaterial;
    }
    else{
        detectorMaterial = G4NistManager::Instance()->FindOrBuildMaterial(fDetectorMaterialName);
    }
    
    G4LogicalVolume* ssdLogic = new G4LogicalVolume(ssdSolid,
                                                    detectorMaterial,
                                                    "ssd.logic");
    
    for(size_t i1=0;i1<4;i1++){
        new G4PVPlacement(0,
                          G4ThreeVector(0.,0.,fDetectorDistance[i1]),
                          ssdLogic,
                          "ssd.physic",
                          worldLogic,
                          false,
                          i1);
        
    }
    
    //** Crystal solid parameters **//
    G4Box* crystalSolid = new G4Box("crystal.solid",
                                    fSizes.x()/2.,
                                    fSizes.y()/2.,
                                    fSizes.z()/2.);



    
    //** Crystal Definition Start **//
    //G4Material* mat;
    G4cout << "Material Name: " << fMaterialName << G4endl;
    if(fMaterialName == ""){
        mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
    }
    else{
        mat = G4NistManager::Instance()->FindOrBuildMaterial(fMaterialName);
    }



   //step limiter
  G4double maxStep = 100.*nm;
    
// **************************************************************
// bandymas primesti amorfini suda

if ( amorfinis == 1)  
	{
	G4LogicalVolume* crystalLogic = new G4LogicalVolume (crystalSolid, mat, "crystal.logic"); 
   	physAbsor = new G4PVPlacement(0,
                      G4ThreeVector(),
                      crystalLogic,
                      "crystal.physic",
                      worldLogic,
                      false,
                      0,
		      true);


    //testuojam step limiteri

  stepLimit = new G4UserLimits(maxStep);
  crystalLogic->SetUserLimits(stepLimit);
    // end of test

	}
else 
{    
    G4ExtendedMaterial* Crystal = new G4ExtendedMaterial("crystal.material",mat);
    
    Crystal->RegisterExtension(std::unique_ptr<G4CrystalExtension>(new G4CrystalExtension(Crystal)));
    G4CrystalExtension* crystalExtension = (G4CrystalExtension*)Crystal->RetrieveExtension("crystal");
   

    // Ciklas nustatyti gardelei modeliuoti
// varis
	if(fMaterialName == "G4_Cu")
	{
    crystalExtension->SetUnitCell(new G4CrystalUnitCell(3.615 * CLHEP::angstrom, //buvo 5.43, varis 3.615
                                                        3.615 * CLHEP::angstrom, //4.805 aliuminis
                                                        3.615 * CLHEP::angstrom, //13.116 
                                                        CLHEP::halfpi,
                                                        CLHEP::halfpi,
                                                        CLHEP::halfpi, //2.0944 * radian for Al2O3
                                                        225)); //Cu - 225, Si - 227, al2o3 - 167
	}
	if(fMaterialName == "G4_Fe")
	{
// alfa-gelezis
    crystalExtension->SetUnitCell(new G4CrystalUnitCell(2.867 * CLHEP::angstrom, //buvo 5.43, varis 3.615
                                                        2.867 * CLHEP::angstrom, //4.805 aliuminis
                                                        2.867 * CLHEP::angstrom, //13.116 
                                                        CLHEP::halfpi,
                                                        CLHEP::halfpi,
                                                        CLHEP::halfpi, //2.0944 * radian for Al2O3
                                                        229)); //Cu - 225, Si - 227, al2o3 - 167
	}

	if(fMaterialName == "G4_Pb")
	{
// alfa-gelezis
    crystalExtension->SetUnitCell(new G4CrystalUnitCell(4.951 * CLHEP::angstrom, //buvo 5.43, varis 3.615
                                                        4.951 * CLHEP::angstrom, //4.805 aliuminis
                                                        4.951 * CLHEP::angstrom, //13.116 
                                                        CLHEP::halfpi,
                                                        CLHEP::halfpi,
                                                        CLHEP::halfpi, //2.0944 * radian for Al2O3
                                                        225)); //Cu - 225, Si - 227, al2o3 - 167
	}

	if(fMaterialName == "G4_W")
	{
// alfa-gelezis
    crystalExtension->SetUnitCell(new G4CrystalUnitCell(3.165 * CLHEP::angstrom, //buvo 5.43, varis 3.615
                                                        3.165 * CLHEP::angstrom, //4.805 aliuminis
                                                        3.165 * CLHEP::angstrom, //13.116 
                                                        CLHEP::halfpi,
                                                        CLHEP::halfpi,
                                                        CLHEP::halfpi, //2.0944 * radian for Al2O3
                                                        229)); //Cu - 225, Si - 227, al2o3 - 167
	}

	if(fMaterialName == "G4_GRAPHITE")
	{
    crystalExtension->SetUnitCell(new G4CrystalUnitCell(2.46 * CLHEP::angstrom, //buvo 5.43, varis 3.615
                                                        2.46 * CLHEP::angstrom, //4.805 aliuminis
                                                        3.35 * CLHEP::angstrom, //13.116 
                                                        CLHEP::halfpi,
                                                        CLHEP::halfpi,
                                                        (CLHEP::pi)/3, //2.0944 * radian for Al2O3
                                                        194)); //Cu - 225, Si - 227, al2o3 - 167
	}
	// jei nei vienas, tada Silicis
	else 
	{
     crystalExtension->SetUnitCell(new G4CrystalUnitCell(5.43 * CLHEP::angstrom,
                                                        5.43 * CLHEP::angstrom,
                                                        5.43 * CLHEP::angstrom,
                                                        CLHEP::halfpi,
                                                        CLHEP::halfpi,
                                                        CLHEP::halfpi,
                                                        227));
	}
	
    Crystal->RegisterExtension(std::unique_ptr<G4ChannelingMaterialData>(new G4ChannelingMaterialData("channeling")));
    G4ChannelingMaterialData* crystalChannelingData = (G4ChannelingMaterialData*)Crystal->RetrieveExtension("channeling");
    crystalChannelingData->SetFilename(fECfileName);

    if(fBRFileName != ""){
        crystalChannelingData->SetBR(fBRFileName);
//        for(int i=0;i<100;i++){
//            G4double z = fSizes.z()*i/100;
//            G4ThreeVector pos = G4ThreeVector(z,z,z);
//            G4cout << i << " " << z << " " << crystalChannelingData->GetBR(pos).x()/CLHEP::meter << G4endl;
//        }
    }
    else if(fBR!=G4ThreeVector()){
        crystalChannelingData->SetBR(fBR.x());
    }

   // while(!getchar());

    G4LogicalCrystalVolume* crystalLogic = new G4LogicalCrystalVolume(crystalSolid,
                                                                      Crystal,
                                                                      "crystal.logic");
    crystalLogic->SetVerbose(1);
    //** Crystal Definition End **//



    G4RotationMatrix* rot = new G4RotationMatrix;
    if(fAngles.x()!=0.){
        rot->rotateX(fAngles.x());
    }
    if(fAngles.y()!=0.){
        rot->rotateY(fAngles.y());
    }
    if(fAngles.z()!=0.){
        rot->rotateZ(fAngles.z());
    }
    
   physAbsor = new G4PVPlacement(rot,
                      G4ThreeVector(),
                      crystalLogic,
                      "crystal.physic",
                      worldLogic,
                      false,
                      0,
		      true);


    //testuojam step limiteri
  //G4double maxStep = 0.2*nm;
  stepLimit = new G4UserLimits(maxStep);
  crystalLogic->SetUserLimits(stepLimit);
    // end of test


} 	// uzdarymas nuo mano else funkcijos    
    

   G4cout<< "Vardas: " << physAbsor->GetName() <<G4endl;
   G4cout<< " Ar amorfinis? " << fCrystalAmorphous << G4endl;




    
    return worldPhysical;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField(){
    if(fCrystalAmorphous == true){
        G4LogicalVolume* crystalLogic = G4LogicalVolumeStore::GetInstance()->GetVolume("crystal.logic");
        G4VSensitiveDetector* crystaldetector = new CrystalDetector("/crystaldetector");
        G4SDManager::GetSDMpointer()->AddNewDetector(crystaldetector);
        crystalLogic->SetSensitiveDetector(crystaldetector);
    }
    else{
        G4LogicalVolume* crystalLogic = G4LogicalVolumeStore::GetInstance()->GetVolume("crystal.logic");
        //G4ChannelingOptrMultiParticleChangeCrossSection* testMany =
        //    new G4ChannelingOptrMultiParticleChangeCrossSection();
        //testMany->AttachTo(crystalLogic);
        //G4cout << " Attaching biasing operator " << testMany->GetName()
        //<< " to logical volume " << crystalLogic->GetName() << G4endl;

        G4VSensitiveDetector* crystaldetector = new CrystalDetector("/crystaldetector");
        G4SDManager::GetSDMpointer()->AddNewDetector(crystaldetector);
        crystalLogic->SetSensitiveDetector(crystaldetector);
    }
        
    G4LogicalVolume* ssdLogic = G4LogicalVolumeStore::GetInstance()->GetVolume("ssd.logic");
    G4VSensitiveDetector* telescope = new SensitiveDetector("/telescope");
    G4SDManager::GetSDMpointer()->AddNewDetector(telescope);
    for(unsigned int i1=0;i1<4;i1++){
        ssdLogic->SetSensitiveDetector(telescope);
    }
    
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

