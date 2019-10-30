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
/// \file EventAction.hh
/// \brief Definition of the EventAction class
//
// $Id: EventAction.hh 76293 2013-11-08 13:11:23Z gcosmo $
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"


//chan
#include <iostream>
#include <fstream>
#include "G4PhysicsFreeVector.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class RunAction;

class EventAction : public G4UserEventAction
{
  public:
    EventAction();
   ~EventAction();

  public:
    virtual void BeginOfEventAction(const G4Event*);
    virtual void   EndOfEventAction(const G4Event*);
    
    void AddEdep (G4double Edep);
    void AddEflow(G4double Eflow);    

    void AddEnergy      (G4double edep)   {EnergyDeposit  += edep;};
    void AddNonIonEnergy(G4double enondep){NonIonEnergyDeposit  += enondep;};
    void AddTheta(G4double tet){theta+=tet;};	
    void AddTrakLenPrim(G4double length) {TrakLenPrim += length;};
    void CountStepsPrim()		 {nbStepsPrim++ ;};
    void AddTrakLenSec(G4double length)  {TrakLenSec  += length;};

    void AddTrueTrakLen(G4double trueLength) {TrueTrakLen += trueLength;};
    void AddProjTrakLen(G4double projLength) {ProjTrakLen += projLength;};

	
	           	  
  private:
    G4double fTotalEnergyDeposit;
    G4double fTotalEnergyFlow;   

    G4double EnergyDeposit;
    G4double NonIonEnergyDeposit;
    G4double theta;	

    G4double TrakLenPrim;
    G4double TrakLenSec;
    G4int nbStepsPrim;
    
    G4double TrueTrakLen;
    G4double ProjTrakLen;

    int hit;

    //chan bullshit
    G4int sdht_ID;
    G4int sdct_ID;


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
