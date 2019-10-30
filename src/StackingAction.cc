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
// StackingAction.cc
// 
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "StackingAction.hh"

#include <iostream>
#include <iomanip>

#include "RunAction.hh"
#include "Run.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "StackingMessenger.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4Material.hh"
#include "G4RunManager.hh"

#include "HistoManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::StackingAction(EventAction* EA, DetectorConstruction* DE )
:eventaction(EA), detector(DE)
{
  killSecondary  = 0;
  stackMessenger = new StackingMessenger(this);
  rec=0;
  SecRec=0;
  absnbrec=0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::~StackingAction()
{
  delete stackMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double  StackingAction::DamageEnergy(G4double T,G4double A, G4double Z)
{
  //.................. T in  eV!!!!!!!!!!!!!
  G4double Z2= Z;
  G4double M2= A;
  G4double k_d;
  G4double epsilon_d;
  G4double g_epsilon_d;
  G4double E_nu;

  k_d=0.1334*std::pow(Z2,(2./3.))*std::pow(M2,(-1./2.));
  epsilon_d=0.01014*std::pow(Z2,(-7./3.))*(T/eV);
  g_epsilon_d= epsilon_d+0.40244*std::pow(epsilon_d,(3./4.))+3.4008*std::pow(epsilon_d,(1./6.));

  E_nu=1./(1.+ k_d*g_epsilon_d);

  return E_nu;//partition fraction!!!
}

//...................................................................

G4ClassificationOfNewTrack
StackingAction::ClassifyNewTrack(const G4Track* aTrack)
{
  G4int IDp= aTrack->GetParentID();
  //G4int TrID = aTrack->GetTrackID();

  G4double energy = aTrack->GetKineticEnergy();
  G4double charge = aTrack->GetDefinition()->GetPDGCharge();
  //G4String parName = aTrack->GetDefinition()->GetParticleName();
   //G4double xpoint  = (aTrack->GetPosition()).x();
   //G4double ypoint  = (aTrack->GetPosition()).y();
   //G4double zpoint  = (aTrack->GetPosition()).z(); 
   //G4double radius = sqrt((xpoint*xpoint)+(ypoint*ypoint));
   G4double Spoint  = (aTrack->GetPosition()).mag();

    // iskopijuota is if ciklo
    G4Material*  material= detector->GetMaterialM();
    G4double A2  = material->GetDensity()/(material->GetTotNbOfAtomsPerVolume()*amu);
    G4double Z2  = material->GetTotNbOfElectPerVolume()/material->GetTotNbOfAtomsPerVolume();

    //Lindhard partition                
    G4double LT = DamageEnergy(energy,A2,Z2);

    //daugyba
    G4double TL=energy*LT;



  ///**************** tesinys

  //keep primary particle
  if (IDp == 0) {

  return fUrgent;}
  Run* run = static_cast<Run*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun());


//tik absorberis
if (aTrack->GetTouchableHandle()->GetVolume() == detector->GetAbsorber()) {
  if (IDp > 0 && charge >0) {
    run->NumberRec(1);
    rec++;
	      
    //total sum of T*L(T)
    run->SumTL(energy*LT); 

    //total sum of T
    run->SumT(energy);

   G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
   analysisManager->FillH1(14, Spoint, TL);
    }

  else if (IDp > 1 && charge > 0 ) { 
   run->NumberSecRec(1);
   run->SumTt(energy);
   SecRec++;
  }}


  //stack or delete secondaries
  G4ClassificationOfNewTrack status = fUrgent;
  if (killSecondary) 
    {if (killSecondary == 1) {
      //add secondary energy before killing
      eventaction->AddEnergy(energy);
      eventaction->AddNonIonEnergy(energy);
      status = fKill;
        }
    }

  return status;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
