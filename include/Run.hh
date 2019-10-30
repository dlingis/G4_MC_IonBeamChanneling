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
/// \file Run.hh
/// \brief Definition of the Run class
//
// $Id: Run.hh 71375 2013-06-14 07:39:33Z maire $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Run_h
#define Run_h 1

#include "G4Run.hh"
#include "G4VProcess.hh"
#include "globals.hh"
#include "G4THitsMap.hh"
#include <map>


//test
class G4ChannelingPhysics;
class DetectorConstruction;
class G4ParticleDefinition;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Run : public G4Run
{
  public:
    Run(DetectorConstruction*);
   ~Run();

  public:
    void SetPrimary(G4ParticleDefinition* particle, G4double energy);         
    void CountProcesses(const G4VProcess* process);
    void ParticleCount(G4String, G4double); 
    void AddEdep (G4double edep);
    void AddEflow (G4double eflow);                   
    void ParticleFlux(G4String, G4double);
    void AddEnergy (G4double);
    void AddNonIonEnergy (G4double);
    void SumTL (G4double);
    void SumT (G4double);
    void SumTt (G4double);
    void NumberRec (G4int);
    void NumberSecRec (G4int);
    void AddNumberOfSteps (G4int);
    void AddTheta (G4double);
    void AddTrakLenPrim (G4double);
    void AddTrakLenSec (G4double);
    //void AddHit() {tothit++;};
    void addh() {hit++;};	// sumuoja detektoriaus hitus
    void countEmerging() {partEmerging++;}; // sumuoja iseinusias is kristalo daleles
    void addkinen(G4double energy) { detKinEn += energy; detKinEn2 += energy*energy;}; // sumuoja kinetines energijas is 4 detektoriaus
    void addrbs() {rbs++;};
    void addsec() {second++;};
    void AddTotStep() {totstep++;};
    void AddNielStep() {nielstep++;};

    void AddProjectedRange(G4double z);



    // addition
    //virtual void RecordEvent(const G4Event*);
    //void DumpData(G4String&) const;
    virtual void Merge(const G4Run*);
    void EndOfRun();     


    void AddTrueRange (G4double l) { fTrueRange += l; fTrueRange2 += l*l;};
    void AddProjRange (G4double x) { fProjRange += x; fProjRange2 += x*x;};


   void AddNonIonisEnergy (G4double niel)   {totniel += niel; totniel2 += niel*niel;};

    //  absorber
   void absLEN (G4double length) {abslen += length; abslen2 = length*length;};
   void absSTP (G4double step)   {absstep += double(step); absstep2 += double(step)*double(step);};
   void absION (G4double edep)   {absfEdep += edep; absfEdep2 += edep*edep;};
   void absNON (G4double niel)   {absniel += niel; absniel2 += niel*niel;};
   void abssumTL(G4double energy){abssum_tl += energy; abssum_tl2 += energy*energy;};
   void abssumT(G4double energy) {abssum_t += energy; abssum_t2 += energy*energy;};
   void absnbrec(G4int rec)	 {absrec += rec;}


  
  private:
    struct ParticleData {
     ParticleData()
       : fCount(0), fEmean(0.), fEmin(0.), fEmax(0.) {}
     ParticleData(G4int count, G4double ekin, G4double emin, G4double emax)
       : fCount(count), fEmean(ekin), fEmin(emin), fEmax(emax) {}
     G4int     fCount;
     G4double  fEmean;
     G4double  fEmin;
     G4double  fEmax;
    };
     
  private:
    DetectorConstruction* fDetector;
    G4ParticleDefinition* fParticle;
    G4double              fEkin;

    //
    G4ChannelingPhysics* chphys;
   

    
    G4double fEnergyDeposit, fEnergyDeposit2;
    G4double fEnergyFlow,    fEnergyFlow2;            
    std::map<G4String,G4int>        fProcCounter;
    std::map<G4String,ParticleData> fParticleDataMap1;                    
    std::map<G4String,ParticleData> fParticleDataMap2;

    G4double        fTrueRange, fTrueRange2;             
    G4double        fProjRange, fProjRange2;

    G4double abslen,abslen2;
    G4double absfEdep,absfEdep2;
    G4double absniel,absniel2;
    G4double absstep,absstep2;
    G4double abssum_tl,abssum_tl2;
    G4double abssum_t,abssum_t2;
    G4int    absrec;

    G4double vInitialTime;



    G4double totniel, totniel2;



  G4double EnergyDeposit,  EnergyDeposit2;
  G4double NonIonEnergyDeposit,  NonIonEnergyDeposit2;
  G4double sum_TL, sum_TL2;
  G4double sum_T, sum_T2, sum_Tt, sum_Tt2;
  G4double Th; 
  G4int N_rec; 
  G4int N_Sec_Rec;
  G4double Nsteps, Nsteps2;
  G4double theta, theta2;
  G4double TrakLenPrim, TrakLenPrim2;
  G4double TrakLenSec, TrakLenSec2;
	G4int hit;
	G4int partEmerging;
	G4double detKinEn, detKinEn2;
	G4int rbs;
	G4int second;
	G4double nielstep;
	G4double totstep;

    G4double projectedR, projectedR2;

  //G4int TotalCount;

 void Print(const std::vector<G4String>& title,
               const std::map< G4int, std::vector<G4double> >&out,
               G4String&) const;
    
    std::map<G4int, G4THitsMap<G4double>* > fMap;
    G4String fOutputFileSpec;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

