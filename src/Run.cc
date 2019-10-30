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
/// \file Run.cc
/// \brief Implementation of the Run class
//
// $Id: Run.cc 71376 2013-06-14 07:44:50Z maire $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Run.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4BetheBlochModel.hh"
#include "G4BraggModel.hh"
#include "G4Proton.hh"

#include "G4hIonisation.hh"

//#include "G4MultiFunctionalDetector.hh"
//#include "G4SDManager.hh"
//#include "G4VPrimitiveScorer.hh"

//keisti neutronu srautui
//#include "RunMessenger.hh"

//#include <fstream>
#include <ctime>


#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run(DetectorConstruction* det)
: G4Run(),
  fDetector(det), fParticle(0), fEkin(0.),  fTrueRange(0.), fTrueRange2(0.),
  fProjRange(0.), fProjRange2(0.), fMap()
{
  //TotalCount = 0.;
  fEnergyDeposit = fEnergyDeposit2 = 0.;
  fEnergyFlow    = fEnergyFlow2    = 0.;  
  EnergyDeposit  = EnergyDeposit2  = 0.;
  NonIonEnergyDeposit = NonIonEnergyDeposit2 = 0.;
  sum_TL = sum_TL2 = 0.;
  sum_T = sum_T2 = 0.;
  sum_Tt = sum_Tt2 = 0.;
  Nsteps = Nsteps2 = 0.;
  theta = theta2 = 0.;
  TrakLenPrim = TrakLenPrim2 = 0.;
  TrakLenSec  = TrakLenSec2  = 0.;
  N_rec = 0.;
  N_Sec_Rec = 0.;
  Th=28.*eV;  
  fTrueRange=0.;
  fTrueRange2=0.;
  fProjRange=0.;
  fProjRange2=0.;
   hit = 0.;
	partEmerging = 0.;
	detKinEn = detKinEn2 = 0.;
	rbs = 0.;
	second = 0.;
	nielstep=0.;
	totstep=0.;
  //

   abslen=abslen2=0.;
   absfEdep=absfEdep2=0.;
   absniel=absniel2=0.;
   absstep=absstep2=0.;
   abssum_tl=abssum_tl2=0.;
   abssum_t=abssum_t2=0.;
   absrec=0.;

   totniel=totniel2=0.;

   projectedR=projectedR2=0.;
        

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run()
{ 
    // Important to clean up the map
    std::map<G4int, G4THitsMap<G4double>* >::iterator iter = fMap.begin();
    
    while (iter != fMap.end()) {
        delete iter->second;
        iter++;}
}
//========================================================
void Run::SetPrimary(G4ParticleDefinition* particle, G4double energy)
{ 
  fParticle = particle;
  fEkin = energy;
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::CountProcesses(const G4VProcess* process) 
{
  G4String procName = process->GetProcessName();
  std::map<G4String,G4int>::iterator it = fProcCounter.find(procName);
  if ( it == fProcCounter.end()) {
    fProcCounter[procName] = 1;
  }
  else {
    fProcCounter[procName]++; 
  }
}

void Run::AddProjectedRange (G4double z) 
{
  projectedR  += z;
  projectedR2 += z*z;
}


void Run::AddEnergy (G4double edep)
{ 
 EnergyDeposit += edep; 
 EnergyDeposit2 += edep*edep; 
}

void Run::AddNonIonEnergy (G4double enondep)
{ 
  NonIonEnergyDeposit += enondep; 
  NonIonEnergyDeposit2 += enondep*enondep; 
}

   // sum of secondary Kinetic energy*L(T):
void Run::SumTL(G4double energyL)
{
 sum_TL +=energyL;  
 sum_TL2 +=energyL*energyL;
}

  // sum of secondary Kinetic energy:
void Run::SumT(G4double energy)
{
  sum_T +=energy;  
  sum_T2 +=energy*energy; 
}
  // sum of tertiary particles energies
void Run::SumTt(G4double energy)
{
  sum_Tt +=energy;  
  sum_Tt2 +=energy*energy; 
}

  //number of recoils
void Run::NumberRec(G4int i)
{ 
  N_rec += i; 
}

  //number of tertiary particles
void Run::NumberSecRec(G4int i)
{
  N_Sec_Rec +=i;
}

  //number of steps
void Run::AddNumberOfSteps(G4int i )
{ 
  Nsteps+=double(i); 
  Nsteps2+=double(i)*double(i);
}
  //scattering angle	
void Run::AddTheta(G4double tet)
{ 
  theta+=tet; 
  theta2+=tet*tet;
}	
  // track length	
void Run::AddTrakLenPrim (G4double length)
{
  //TotalCount++;
  //
  TrakLenPrim += length; 
  TrakLenPrim2 += length*length;
}
void Run::AddTrakLenSec (G4double length)
{
  TrakLenSec += length; 
  TrakLenSec2 += length*length;
}
  
                  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::ParticleCount(G4String name, G4double Ekin)
{
  std::map<G4String, ParticleData>::iterator it = fParticleDataMap1.find(name);
  if ( it == fParticleDataMap1.end()) {
    fParticleDataMap1[name] = ParticleData(1, Ekin, Ekin, Ekin);
  }
  else {
    ParticleData& data = it->second;
    data.fCount++;
    data.fEmean += Ekin;
    //update min max
    G4double emin = data.fEmin;
    if (Ekin < emin) data.fEmin = Ekin;
    G4double emax = data.fEmax;
    if (Ekin > emax) data.fEmax = Ekin; 
  }   
}
                 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::AddEdep(G4double edep)
{ 
  fEnergyDeposit += edep;
  fEnergyDeposit2 += edep*edep;
}
                 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::AddEflow(G4double eflow)
{ 
  fEnergyFlow += eflow;
  fEnergyFlow2 += eflow*eflow;
}                  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::ParticleFlux(G4String name, G4double Ekin)
{
  std::map<G4String, ParticleData>::iterator it = fParticleDataMap2.find(name);
  if ( it == fParticleDataMap2.end()) {
    fParticleDataMap2[name] = ParticleData(1, Ekin, Ekin, Ekin);
  }
  else {
    ParticleData& data = it->second;
    data.fCount++;
    data.fEmean += Ekin;
    //update min max
    G4double emin = data.fEmin;
    if (Ekin < emin) data.fEmin = Ekin;
    G4double emax = data.fEmax;
    if (Ekin > emax) data.fEmax = Ekin; 
  }   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run* run)
{
  const Run* localRun = static_cast<const Run*>(run);
  
  //primary particle info
  //
  fParticle = localRun->fParticle;
  fEkin     = localRun->fEkin;
  
  // accumulate sums
  //
  fEnergyDeposit   += localRun->fEnergyDeposit;  
  fEnergyDeposit2  += localRun->fEnergyDeposit2;
  fEnergyFlow      += localRun->fEnergyFlow;
  fEnergyFlow2     += localRun->fEnergyFlow2;

  //additional parameters
  //TotalCount += localRun->TotalCount;
  //
  EnergyDeposit    += localRun->EnergyDeposit;  
  EnergyDeposit2   += localRun->EnergyDeposit2;  
  NonIonEnergyDeposit += localRun->NonIonEnergyDeposit;
  NonIonEnergyDeposit2 += localRun->NonIonEnergyDeposit2;
  sum_TL	   += localRun->sum_TL; 
  sum_TL2	   += localRun->sum_TL2;
  sum_T	  	   += localRun->sum_T; 
  sum_T2	   += localRun->sum_T2; 
  sum_Tt	   += localRun->sum_Tt; 
  sum_Tt2	   += localRun->sum_Tt2; 
  N_rec		   += localRun->N_rec; 
  N_Sec_Rec 	   += localRun->N_Sec_Rec;
  Nsteps	   += localRun->Nsteps;
  Nsteps2	   += localRun->Nsteps2;
  theta		   += localRun->theta;
  theta2	   += localRun->theta2;
  TrakLenPrim      += localRun->TrakLenPrim;
  TrakLenPrim2 	   += localRun->TrakLenPrim2;
  TrakLenSec       += localRun->TrakLenSec;
  TrakLenSec2 	   += localRun->TrakLenSec2;

  fTrueRange  += localRun->fTrueRange;
  fTrueRange2 += localRun->fTrueRange2;
  fProjRange  += localRun->fProjRange;
  fProjRange2 += localRun->fProjRange2;

  // ABSORBER
  abslen	+= localRun->abslen;
  abslen2 	+= localRun->abslen2;
  absstep	+= localRun->absstep;
  absstep2	+= localRun->absstep2;
  absfEdep      += localRun->absfEdep;
  absniel	+= localRun->absniel;
  absfEdep2     += localRun->absfEdep2;
  absniel2	+= localRun->absniel2;
  abssum_tl	+= localRun->abssum_tl;
  abssum_tl2	+= localRun->abssum_tl2;
  abssum_t	+= localRun->abssum_t;
  abssum_t2	+= localRun->abssum_t2;
  absrec 	+= localRun->absrec;


  totniel       += localRun->totniel;
  totniel2      += localRun->totniel2;

  hit		+= localRun->hit;
  partEmerging  += localRun->partEmerging;

  detKinEn 	+= localRun->detKinEn;
  detKinEn2	+= localRun->detKinEn2;

  rbs 		+= localRun->rbs;
  second 	+= localRun->second;
  nielstep	+= localRun->nielstep;
  totstep	+= localRun->totstep;


  projectedR	+=localRun->projectedR;
  projectedR2	+=localRun->projectedR2;

  //tothit 	+= localRun->hit;
  // END

  //map: processes count
  std::map<G4String,G4int>::const_iterator itp;
  for ( itp = localRun->fProcCounter.begin();
        itp != localRun->fProcCounter.end(); ++itp ) {

    G4String procName = itp->first;
    G4int localCount = itp->second;
    if ( fProcCounter.find(procName) == fProcCounter.end()) {
      fProcCounter[procName] = localCount;
    }
    else {
      fProcCounter[procName] += localCount;
    }  
  }
  
  //map: created particles count    
  std::map<G4String,ParticleData>::const_iterator itc;
  for (itc = localRun->fParticleDataMap1.begin(); 
       itc != localRun->fParticleDataMap1.end(); ++itc) {
    
    G4String name = itc->first;
    const ParticleData& localData = itc->second;   
    if ( fParticleDataMap1.find(name) == fParticleDataMap1.end()) {
      fParticleDataMap1[name]
       = ParticleData(localData.fCount, 
                      localData.fEmean, 
                      localData.fEmin, 
                      localData.fEmax);
    }
    else {
      ParticleData& data = fParticleDataMap1[name];   
      data.fCount += localData.fCount;
      data.fEmean += localData.fEmean;
      G4double emin = localData.fEmin;
      if (emin < data.fEmin) data.fEmin = emin;
      G4double emax = localData.fEmax;
      if (emax > data.fEmax) data.fEmax = emax; 
    }   
  }
  
  //map: particles flux count       
  std::map<G4String,ParticleData>::const_iterator itn;
  for (itn = localRun->fParticleDataMap2.begin(); 
       itn != localRun->fParticleDataMap2.end(); ++itn) {
    
    G4String name = itn->first;
    const ParticleData& localData = itn->second;   
    if ( fParticleDataMap2.find(name) == fParticleDataMap2.end()) {
      fParticleDataMap2[name]
       = ParticleData(localData.fCount, 
                      localData.fEmean, 
                      localData.fEmin, 
                      localData.fEmax);
    }
    else {
      ParticleData& data = fParticleDataMap2[name];   
      data.fCount += localData.fCount;
      data.fEmean += localData.fEmean;
      G4double emin = localData.fEmin;
      if (emin < data.fEmin) data.fEmin = emin;
      G4double emax = localData.fEmax;
      if (emax > data.fEmax) data.fEmax = emax; 
    }   
  }
    const std::map< G4int, G4THitsMap<G4double>* >& localMap = localRun->fMap;
    std::map< G4int, G4THitsMap<G4double>* >::const_iterator iter = localMap.begin();
    for ( ; iter != localMap.end() ; ++iter)
        (*(fMap[iter->first])) += (*(iter->second));


  G4Run::Merge(run); 
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::EndOfRun() 
{


	// papildoma isvestis i text faila

	//ofstream myfile("output.txt");


  G4int prec = 5, wid = prec + 2;  
  G4int dfprec = G4cout.precision(prec);
  //fDetector->GetMaterial()->GetName()
  //run condition
  //

  G4Material* material = fDetector->GetMaterialM();
  G4double density = material->GetDensity();
   
  G4String Particle = fParticle->GetParticleName();    
  G4cout << "\n The run is " << numberOfEvent << " "<< Particle << " of "
         << G4BestUnit(fEkin,"Energy") << " through "  << G4BestUnit(fDetector->GetLength(),"Length") << " length parallelpiped of " << material->GetName() << " (density: " << G4BestUnit(density,"Volumic Mass") << ")" << G4endl;

  if (numberOfEvent == 0) { G4cout.precision(dfprec);   return;}
             
  //frequency of processes
  //
  G4cout << "\n Process calls frequency :" << G4endl;
  G4int index = 0;
  std::map<G4String,G4int>::iterator it;    
  for (it = fProcCounter.begin(); it != fProcCounter.end(); it++) {
     G4String procName = it->first;
     G4int    count    = it->second;
     G4String space = " "; if (++index%3 == 0) space = "\n";
     G4cout << " " << std::setw(20) << procName << "="<< std::setw(7) << count
            << space;
  }
  G4cout << G4endl;
  
  //particles count
  //
  G4cout << "\n List of generated particles:" << G4endl;
     
 std::map<G4String,ParticleData>::iterator itc;               
 for (itc = fParticleDataMap1.begin(); itc != fParticleDataMap1.end(); itc++) { 
    G4String name = itc->first;
    ParticleData data = itc->second;
    G4int count = data.fCount;
    G4double eMean = data.fEmean/count;
    G4double eMin = data.fEmin;
    G4double eMax = data.fEmax;    
         
    G4cout << "  " << std::setw(13) << name << ": " << std::setw(7) << count
           << "  Emean = " << std::setw(wid) << G4BestUnit(eMean, "Energy")
           << "\t( "  << G4BestUnit(eMin, "Energy")
           << " --> " << G4BestUnit(eMax, "Energy") 
           << ")" << G4endl;           
 }
   
  // compute mean Energy deposited and rms
  //
  G4int TotNbofEvents = numberOfEvent;
  fEnergyDeposit /= TotNbofEvents; fEnergyDeposit2 /= TotNbofEvents;
  G4double rmsEdep = fEnergyDeposit2 - fEnergyDeposit*fEnergyDeposit;
  if (rmsEdep>0.) rmsEdep = std::sqrt(rmsEdep);
  else            rmsEdep = 0.;
  
  G4cout << "\n Mean energy deposit per event = "
         << G4BestUnit(fEnergyDeposit,"Energy") << ";  rms = "
         << G4BestUnit(rmsEdep,      "Energy") 
         << G4endl;
  
  // compute mean Energy flow and rms
  //
  fEnergyFlow /= TotNbofEvents; fEnergyFlow2 /= TotNbofEvents;
  G4double rmsEflow = fEnergyFlow2 - fEnergyFlow*fEnergyFlow;
  if (rmsEflow>0.) rmsEflow = std::sqrt(rmsEflow);
  else             rmsEflow = 0.;
  
  G4cout << " Mean energy flow per event    = "
         << G4BestUnit(fEnergyFlow,"Energy") << ";  rms = "
         << G4BestUnit(rmsEflow,   "Energy") 
         << G4endl;
                                
 //particles flux
 //
 G4cout << "\n List of particles emerging from the absorber :" << G4endl;
     
 std::map<G4String,ParticleData>::iterator itn;               
 for (itn = fParticleDataMap2.begin(); itn != fParticleDataMap2.end(); itn++) { 
    G4String name = itn->first;
    ParticleData data = itn->second;
    G4int count = data.fCount;
    G4double eMean = data.fEmean/count;
    G4double eMin = data.fEmin;
    G4double eMax = data.fEmax;
    G4double Eflow = data.fEmean/TotNbofEvents;        
         
    G4cout << "  " << std::setw(13) << name << ": " << std::setw(7) << count
           << "  Emean = " << std::setw(wid) << G4BestUnit(eMean, "Energy")
           << "\t( "  << G4BestUnit(eMin, "Energy")
           << " --> " << G4BestUnit(eMax, "Energy") 
           << ") \tEflow/event = " << G4BestUnit(Eflow, "Energy") << G4endl;
 }

  // compute mean and rms
  if (TotNbofEvents == 0) return;
  
  //total energy loss
  EnergyDeposit /= TotNbofEvents; EnergyDeposit2 /= TotNbofEvents;
  G4double FrmsEdep = EnergyDeposit2 - EnergyDeposit*EnergyDeposit;
  if(FrmsEdep >0.)FrmsEdep=std::sqrt(FrmsEdep/TotNbofEvents);
  else FrmsEdep =0;

  //nuclear energy loss
  NonIonEnergyDeposit /= TotNbofEvents; NonIonEnergyDeposit2 /= TotNbofEvents;
  G4double rmsEnondep = NonIonEnergyDeposit2 - NonIonEnergyDeposit*NonIonEnergyDeposit;
  if(rmsEnondep>0.) rmsEnondep= std::sqrt(rmsEnondep/TotNbofEvents);
  else rmsEnondep=0;

  //mean sum of T( kinetic energy of secondary)x L(T) partition energy 
  sum_TL/=TotNbofEvents;     sum_TL2/=TotNbofEvents;
  G4double rmssum_TL =sum_TL2- sum_TL*sum_TL;
  if(rmssum_TL>0.) rmssum_TL=std::sqrt(rmssum_TL/TotNbofEvents);
  else rmssum_TL =0;

  //mean kinetic energy of secondary particles (IDp==1) 
  G4double rmssum_T = 0.0;
  if(N_rec > 0) {
    sum_T/=N_rec;     sum_T2/=N_rec;
    rmssum_T =sum_T2- sum_T*sum_T;
    if(rmssum_T>0.) rmssum_T=std::sqrt(rmssum_T/N_rec);  }

  //mean kinetic energy of tertiary particles (IDp>1) 
  G4double rmssum_Tt = 0.0;
  if(N_Sec_Rec > 0) {
    sum_Tt/=N_Sec_Rec;     sum_Tt2/=N_Sec_Rec;
    rmssum_Tt =sum_Tt2- sum_Tt*sum_Tt;
    if(rmssum_Tt>0.) rmssum_Tt=std::sqrt(rmssum_Tt/N_Sec_Rec);  }

  //mean number of steps:
  Nsteps/=TotNbofEvents;  Nsteps2/=TotNbofEvents;
  G4double rmsSteps= Nsteps2 -Nsteps*Nsteps;
  if(rmsSteps>0) rmsSteps= std::sqrt(rmsSteps/TotNbofEvents);
  else rmsSteps=0;

  //scattering angle
  theta/=TotNbofEvents ; theta2/=TotNbofEvents;	
  G4double rmsTheta =theta2-theta*theta;
  if(rmsTheta>0.) rmsTheta =std::sqrt(rmsTheta/TotNbofEvents);
  else rmsTheta =0;
  
  //track length
  TrakLenPrim /= TotNbofEvents; TrakLenPrim2 /= TotNbofEvents;
  G4double rmsTLPrim = TrakLenPrim2 - TrakLenPrim*TrakLenPrim;
  if (rmsTLPrim>0.) rmsTLPrim = std::sqrt(rmsTLPrim/TotNbofEvents);
  else rmsTLPrim = 0.;
  // secondaries track length
  TrakLenSec /= N_rec; TrakLenSec2 /= N_rec;
  G4double rmsTLSec = TrakLenSec2 - TrakLenSec*TrakLenSec;
  if (rmsTLSec>0.) rmsTLSec = std::sqrt(rmsTLSec/N_rec);
  else rmsTLSec = 0.;

  // AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
  //************************************************
  // --------------------------------------------
  // absorber 
  //  edep
  absfEdep /= TotNbofEvents; absfEdep2 /= TotNbofEvents;
  G4double absrmsEdep = absfEdep2 - absfEdep*absfEdep;
  if(absrmsEdep >0.)absrmsEdep=std::sqrt(absrmsEdep/TotNbofEvents);
  else absrmsEdep =0;
  // track length 
  abslen /= TotNbofEvents; abslen2 /= TotNbofEvents;
  G4double absrmsTLPrim = abslen2 - abslen*abslen;
  if (absrmsTLPrim>0.) absrmsTLPrim = std::sqrt(absrmsTLPrim/TotNbofEvents);
  else absrmsTLPrim = 0.;
  //nuclear energy loss
  absniel /= TotNbofEvents; absniel2 /= TotNbofEvents;
  G4double absrmsEnondep = absniel2 - absniel*absniel;
  if(absrmsEnondep>0.) absrmsEnondep= std::sqrt(absrmsEnondep/TotNbofEvents);
  else absrmsEnondep=0;
  //mean sum of T( kinetic energy of secondary)x L(T) partition energy 
  abssum_tl/=TotNbofEvents;     abssum_tl2/=TotNbofEvents;
  G4double absrmssum_tl =abssum_tl2- abssum_tl*abssum_tl;
  if(absrmssum_tl>0.) absrmssum_tl=std::sqrt(absrmssum_tl/TotNbofEvents);
  else absrmssum_tl =0;
  //mean kinetic energy of secondary particles (IDp==1) 
  G4double absrmssum_t = 0.0;
  if(absrec > 0) {
    abssum_t/=absrec;     abssum_t2/=absrec;
    absrmssum_t =abssum_t2- abssum_t*abssum_t;
    if(absrmssum_t>0.) absrmssum_t=std::sqrt(absrmssum_t/absrec);  }
  //mean number of steps:
  absstep/=TotNbofEvents;  absstep2/=TotNbofEvents;
  G4double absrmsSteps= absstep2 -absstep*absstep;
  if(absrmsSteps>0) absrmsSteps= std::sqrt(absrmsSteps/TotNbofEvents);
  else absrmsSteps=0;
  // **************************************
  // end of absorber
  //*************************************************

  //nuclear energy loss
  totniel /= TotNbofEvents; totniel2 /= TotNbofEvents;
  G4double totrmsEnondep = totniel2 - totniel*totniel;
  if(totrmsEnondep>0.) totrmsEnondep= std::sqrt(totrmsEnondep/TotNbofEvents);
  else totrmsEnondep=0;


//totniel
  // PLEASE END MY MYSERY
  //..............................................................

  //G4double thickness  = fDetector->GetRadius();
  //Stopping Power and NIEL from simulation.
  // effective length
  G4double length=TrakLenPrim;

  // total energy loss  
  G4double meandEdx  = EnergyDeposit/length;
  // nuclear energy loss
  G4double meandEdx_nucl  = NonIonEnergyDeposit/length;
  // NIEL 
  G4double meandEdx_sumTL=sum_TL/length;

  //G4double RealMFP = TrakLenPrim/TotalCount;

	//[MeVcm2/g]
  G4double stopPower = meandEdx/density;  
  G4double stopPower_nucl = meandEdx_nucl/density;
  G4double stopPower_sumTL=meandEdx_sumTL/density;
 
 	//mean free path & corss section 
  G4double freepath= TrakLenPrim/Nsteps;
  G4double er1=rmsTLPrim/Nsteps;	
  G4double er2=freepath*rmsSteps/Nsteps;
  G4double rmsFreepath=std::sqrt(er1*er1+er2*er2);


  G4double NA  = material->GetTotNbOfAtomsPerVolume();
  G4double CrossSection =1./(NA*freepath); 
  G4double rmsCrossSection=rmsFreepath*CrossSection/freepath;


  // true and projected range from TestEm1
  fTrueRange /= TotNbofEvents; fTrueRange2 /= TotNbofEvents;
  G4double trueRms = fTrueRange2 - fTrueRange*fTrueRange;        
  if (trueRms>0.) trueRms = std::sqrt(trueRms); else trueRms = 0.;
        
  fProjRange /= TotNbofEvents; fProjRange2 /= TotNbofEvents;
  G4double projRms = fProjRange2 - fProjRange*fProjRange;        
  if (projRms>0.) projRms = std::sqrt(projRms); else projRms = 0.;

	// projected range from testem11
  //compute projected range of primary track
  //
  G4double rmsPR;
  projectedR /= TotNbofEvents; projectedR2 /= TotNbofEvents;
  rmsPR = projectedR2 - projectedR*projectedR;        
  if (rmsPR>0.) rmsPR = std::sqrt(rmsPR); else rmsPR = 0.;




  G4cout << "\n ============= Simulation statistics ==============\n";
  G4cout << "\n Primary Total track length in absorber:\n "
         << G4BestUnit(TrakLenPrim,"Length") << " +- "
         << G4BestUnit(rmsTLPrim,       "Length") << G4endl;

  G4cout << "\n Secondaries total track length in absorber:\n "
         << G4BestUnit(TrakLenSec,"Length") << " +- "
         << G4BestUnit(rmsTLSec,       "Length") << G4endl;

  G4cout << "\n ==================================================\n ";
  G4cout << "\n Primary particle statistics\n ";
  G4cout << "\n Mean Number of Steps:\n "<<Nsteps<< " +/- "<< rmsSteps<<G4endl;
								
  G4cout << "\n Mean Free Path :\n "<<G4BestUnit(freepath,"Length")<<  
				" +/- "<< G4BestUnit(rmsFreepath,"Length")<<G4endl;

  G4cout << "\n Mean Cross Section :\n "<<G4BestUnit(CrossSection,"Surface")<<
			" +/- "<<G4BestUnit(rmsCrossSection,"Surface")<<G4endl;

  G4cout << "\n Mean scattering angle :\n "<<G4BestUnit(theta,"Angle")<< " +/- "
			<< G4BestUnit(rmsTheta,"Angle")<<G4endl;
  
  G4cout << "\n Total energy deposit in absorber:\n "
         << G4BestUnit(EnergyDeposit,"Energy") << " +/- "
         << G4BestUnit(FrmsEdep,      "Energy") 
         << G4endl;
  G4cout << "-----> dE/dx total= " << meandEdx/(MeV/cm) << " MeV/cm"
         << "\t(" << stopPower/(MeV*cm2/g) << " MeV*cm2/g)"
         << G4endl;


  G4cout << "\n Nuclear energy deposit in absorber:\n "
         << G4BestUnit(NonIonEnergyDeposit,"Energy") << " +/- "
         << G4BestUnit(rmsEnondep,      "Energy")
         << G4endl;
  G4cout << "-----> dE/dx  nucl = " << meandEdx_nucl/(MeV/cm) << " MeV/cm"
         << "\t(" << stopPower_nucl/(MeV*cm2/g) << " MeV*cm2/g)"
         << G4endl;

  G4cout <<"\n NIEL in absorber (Th>"<<Th/eV <<" eV):\n "
         << G4BestUnit(sum_TL,"Energy") << " +/- "
         << G4BestUnit(rmssum_TL,      "Energy")
         << G4endl;
  G4cout << "-----> NIEL = " << meandEdx_sumTL/(MeV/cm) << " MeV/cm"
         << "\t(" << stopPower_sumTL/(MeV*cm2/g) << " MeV*cm2/g)"
         << G4endl;
/*
  // Srauto vertinimas 
  //G4double l=	1.;
  G4double a=   fDetector->GetLength();
  G4double area=a*a;//cm2
  G4double flu =TotNbofEvents/area; // #/cm2
  G4double flux =flu*cm2;
  G4double Tdam = sum_TL/eV;
  G4cout<<"\n============================================================="<<G4endl;
  G4cout<<"\n Total Number of primary knock-on atoms (PKA): PKA = "<< N_rec<<G4endl;
  G4cout<<" PKA per primary particle = "<< N_rec/TotNbofEvents<<G4endl;
  G4cout<<" Mean Kinetic Energy of PKA = "<<G4BestUnit(sum_T,"Energy")
                <<" +/- "  <<G4BestUnit(rmssum_T,"Energy")<<G4endl; 
  G4cout<<" Total flux model : " <<flux<< " #/cm2"<<G4endl;
  G4cout<<" T damage: " <<Tdam<< " eV "<<G4endl;
  G4cout<<" Total number of secondary knock-ons: " <<N_Sec_Rec<<G4endl;
  G4cout<<" Mean Kinetic Energy of SKA = "<<G4BestUnit(sum_Tt,"Energy")
                <<" +/- "  <<G4BestUnit(rmssum_Tt,"Energy")<<G4endl;*/
  G4cout << "\n ===========================================================\n";
  G4cout << "\n true Range = " << G4BestUnit(fTrueRange,"Length")
         << "   rms = "        << G4BestUnit(trueRms,  "Length");

  G4cout << "\n proj Range = " << G4BestUnit(fProjRange,"Length")
         << "   rms = "        << G4BestUnit(projRms,  "Length");
  G4cout << "\n ===========================================================\n";
  G4cout << " ================= GEOMETRIJOS parametrai==================";
  G4cout << "\n ===========================================================\n";
  G4cout << " Sluoksnio storis = " << G4BestUnit(fDetector->GetLength(),"Length")<<G4endl;
  G4cout << " Bendri detektoriaus matmenys : " << G4BestUnit(fDetector->GetSizes(), "Length") << G4endl;	
  G4cout << " Medziaga: " << fDetector->GetMaterial() << G4endl;
  G4cout << " **************************************************************\n";
  G4cout << " ***************** Detektoriai ********************************\n";
  G4cout << " Detektoriu matmenys: "<< G4BestUnit(fDetector->GetDetectorSizes(), "Length") << G4endl;
  G4cout << " Detektoriu medziaga: "<< fDetector->GetDetectorMaterial() << G4endl;
  G4cout << " Detektoriu pozicijos: "<< G4BestUnit(fDetector->GetDetectorDistance(0), "Length") << "  " << G4BestUnit(fDetector->GetDetectorDistance(1), "Length") << "  " << G4BestUnit(fDetector->GetDetectorDistance(2), "Length") << "  " << G4BestUnit(fDetector->GetDetectorDistance(3), "Length") << G4endl;


  //G4double vFinalTime = time(NULL); 
  //G4cout << " Calculation time in minutes: " << vFinalTime-vInitialTime/60. << G4endl; 
/*
	// dEdx skaiciavimas
  G4ParticleDefinition* prot = G4Proton::Proton();
  G4DataVector cuts;
  cuts.push_back(1*keV);


  G4double Z = material->GetZ();
  G4double A = material->GetA();
	G4double Emin, Emax, dE, Ecut;
	G4double Eminb, Emaxb, dEb, Ecutb;


  G4VEmModel* ioni = new G4BetheBlochModel();
  ioni->Initialise(prot, cuts);
  
  // compute CrossSection per atom and MeanFreePath
  //
  Emin = 1.01*MeV; Emax = 20.01*MeV; dE = 500*keV;
  Ecut = 100*keV;

  G4cout << "\n #### proton : CrossSection, MeanFreePath and StoppingPower"
         << " for " << material->GetName() 
         << ";\tEnergy cut = " << G4BestUnit (Ecut, "Energy") << G4endl;
         
  G4cout << "\n Energy \t ionization \t";
  G4cout <<           "\t ionization \t";
  G4cout <<           "\t ionization" << G4endl;
  
  for (G4double Energy = Emin; Energy <= Emax; Energy += dE) {
    G4cout << "\n " << G4BestUnit (Energy, "Energy")
     << "\t" 
     << G4BestUnit (ioni->ComputeCrossSectionPerAtom(prot,Energy,Z,A,Ecut),
                   "Surface")
     << "\t \t"         
     << G4BestUnit (ioni->ComputeMeanFreePath(prot,Energy,material,Ecut),
                   "Length")           
     << "\t \t"         
     << G4BestUnit (ioni->ComputeDEDXPerVolume(material,prot,Energy,Ecut),
                   "Energy/Length");
  }



  // low energy : Bragg Model
  ioni = new G4BraggModel(prot);
  ioni->Initialise(prot, cuts);
  
  // compute CrossSection per atom and MeanFreePath
  //
  Eminb = 1.*keV; Emaxb = 2.01*MeV; dEb = 200*keV;
  Ecutb = 10*keV;
  
  G4cout << "\n #### proton : low energy model (Bragg) "
         << ";\tEnergy cut = " << G4BestUnit (Ecutb, "Energy") << G4endl;
                           
  G4cout << "\n Energy \t ionization \t";
  G4cout <<           "\t ionization \t";
  G4cout <<           "\t ionization" << G4endl;
  
  for (G4double Energy = Eminb; Energy <= Emaxb; Energy += dEb) {
    G4cout << "\n " << G4BestUnit (Energy, "Energy")
     << "\t" 
     << G4BestUnit (ioni->ComputeCrossSectionPerAtom(prot,Energy,Z,A,Ecutb),
                   "Surface")
     << "\t \t"         
     << G4BestUnit (ioni->ComputeMeanFreePath(prot,Energy,material,Ecutb),
                   "Length")           
     << "\t \t"         
     << G4BestUnit (ioni->ComputeDEDXPerVolume(material,prot,Energy,Ecutb),
                   "Energy/Length");
  }
*/

  G4double chan = hit/4;
  G4double totchan = chan/TotNbofEvents;
  G4double totem = chan/partEmerging;
  // kinetines energijos
  // kinetic energy of particles reaching the 4th detector
  G4double rmssum_kinen = 0.0;
  if(chan > 0) {
    detKinEn/=hit;     detKinEn2/=hit;
    rmssum_kinen =detKinEn2- detKinEn*detKinEn;
    if(rmssum_kinen>0.) rmssum_kinen=std::sqrt(rmssum_kinen/chan);  }
  //


  G4cout << " \n**************************************************** " << G4endl;
  G4cout 
    << "\n Projected range = " << G4BestUnit(projectedR,"Length")
    << " +- " << G4BestUnit( rmsPR,"Length")<< G4endl;
    
  G4cout << " \n**************************************************** " << G4endl;
  G4cout << " # of particles that reaches last detector " << hit/4 << G4endl;
  G4cout << " part of total particles " << std::fixed << std::setprecision(2) << totchan*100 << " % " <<G4endl;
  G4cout << " part of emitted particles " << std::setprecision(2) << totem*100 << " % "<< G4endl;

  G4cout << " Mean KinEn of particles reaching 4th detector = "<<G4BestUnit(detKinEn,"Energy")
                <<" +/- "  <<G4BestUnit(rmssum_kinen,"Energy")<<G4endl; 
  G4cout << " ****************************************************** " << G4endl;
  G4cout << " # of backscattered primaries " << rbs << " see backscattered.txt for more info. " << G4endl;
  G4cout << " # of secondaries generated " << second << " see antrines.txt for more info. " << G4endl;
  G4cout << " ****************************************************** " << G4endl;
  G4cout << " ****************************************************** " << G4endl;
  G4cout << " ****************************************************** " << G4endl;
  G4cout << " DON'T FORGET TO FLUSH antrines AND backscattered FILES " << G4endl;
  G4cout << " ****************************************************** " << G4endl;
  G4cout << " ****************************************************** " << G4endl;
  G4cout << " ****************************************************** " << G4endl;

 
  //G4cout << " Kritinis kampas: " << G4BestUnit(GetCriticalAngle(), "Angle") << G4endl;
  //G4double stepMax = fPhysics->GetStepMaxProcess()->GetMaxStep();
  //G4double step = StepMax()->GetMaxStep(); 
  //G4double stepMax = chphys->GetStepMaxProcess()->GetMaxStep();
  //G4cout << " Step : " << stepMax << G4endl;

  //........................................................
 

  G4cout << " # of niel steps " << nielstep << G4endl;
  G4cout << " # of total steps " << totstep << G4endl;

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance(); 
	/// modifikuota is testEm7
  //if (analysisManager->IsActive()) {
   // analysisManager->OpenFile(); 
/*
    // histogram "1" is defined by the length of the target
    // zoomed histograms are defined by UI command  
    G4double DETlength  = fDetector->GetLength();
    //G4int nbmin = 1000;
    analysisManager->SetH1(17, 1000, 0.,DETlength,"um"); // mano pridetas
  	*/


  //normalize histograms
  //G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();        
  for (G4int ih=1; ih<16; ih++) {
	G4double fac;

    G4double binWidth = analysisManager->GetH1Width(ih);
    G4double unit     = analysisManager->GetH1Unit(ih);  
    /*G4double*/ fac = unit/binWidth;
    ///G4double fac = unit/(numberOfEvent*binWidth);
    if (ih == 2) fac = (1./(numberOfEvent*binWidth))*(mm/MeV);  
    //if (ih == 19) fac = (1./(numberOfEvent*binWidth))*(mm/MeV);  
    //if (ih == 16) fac = (1./(Nsteps*TotNbofEvents));
    if (ih == 15) {fac = (1./totstep);analysisManager->ScaleH1(ih,fac);};	
    if (ih == 16) {fac = (1./numberOfEvent); analysisManager->ScaleH1(ih,fac);};
    if (ih == 17) {fac = (1./(numberOfEvent*binWidth))*(mm/keV); analysisManager->ScaleH1(ih,fac);}; //buvo nielstep
    if (ih == 18) {fac = (1./(numberOfEvent*binWidth)); analysisManager->ScaleH1(ih,fac);};//buvo nielstep
    if (ih == 19) {fac = (1./(totstep*binWidth)); analysisManager->ScaleH1(ih,fac);};
    analysisManager->ScaleH1(ih,fac);
  }  

/*
		  " Sklaidos kampas ",						//15
		  " Projected range distribution ",				//16	->
		  " NIEL pagal atstuma, keV	",				//17 
		  " NIEL pagal energija",					//18
		  " EDEP pagal energija "					//19
           */
  //remove all contents in fProcCounter, fCount 
  fProcCounter.clear();
  fParticleDataMap2.clear();
                          
  //restore default format         
  G4cout.precision(dfprec);   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
