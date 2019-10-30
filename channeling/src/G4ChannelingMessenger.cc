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

#include "G4ChannelingMessenger.hh"
#include "G4Channeling.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4RunManager.hh"

#include "G4ios.hh"

G4ChannelingMessenger::
G4ChannelingMessenger(G4Channeling* mpga)
:fTarget(mpga){
    fChannelingDirectory = new G4UIdirectory("/process/channeling/");
    fChannelingDirectory->SetGuidance("Channeling Process control commands.");
    
    G4String command;
    
    command = "/process/channeling/set_min_energy";
    fMinimumEnergyCmd = new G4UIcmdWithADoubleAndUnit(command,this);
    fMinimumEnergyCmd->SetGuidance("Set minimum energy allowed.");
    fMinimumEnergyCmd->SetParameterName("min_energy",
                                        true);
    fMinimumEnergyCmd->SetDefaultValue(1.);
    fMinimumEnergyCmd->SetDefaultUnit("MeV");
    
    command = "/process/channeling/set_max_mom_ratio";
    fMaximumMomentumRatioCmd = new G4UIcmdWithADouble(command,this);
    fMaximumMomentumRatioCmd->SetGuidance("Set maximum momentum ratio allowed.");
    fMaximumMomentumRatioCmd->SetParameterName("max_mom_ratio",
                                               true);
    fMaximumMomentumRatioCmd->SetDefaultValue(0.01);
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ChannelingMessenger::
~G4ChannelingMessenger(){
    delete fMinimumEnergyCmd;
    delete fMaximumMomentumRatioCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ChannelingMessenger::SetNewValue(G4UIcommand *command,
                                        G4String newValue){
    if(command==fMinimumEnergyCmd){
        fTarget->SetMinimumEnergy(fMinimumEnergyCmd->GetNewDoubleValue(newValue));
    }
    if(command==fMaximumMomentumRatioCmd){
        fTarget->SetMaximumMomentumRatio(fMaximumMomentumRatioCmd->GetNewDoubleValue(newValue));
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4String G4ChannelingMessenger::GetCurrentValue(G4UIcommand * command){
    G4String cv;
    
    if( command==fMinimumEnergyCmd ){
        cv = fMinimumEnergyCmd->ConvertToString(fTarget->GetMinimumEnergy(),"MeV");
    }
    if( command==fMaximumMomentumRatioCmd ){
        cv = fMaximumMomentumRatioCmd->ConvertToString(fTarget->GetMaximumMomentumRatio());
    }
    return cv;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
