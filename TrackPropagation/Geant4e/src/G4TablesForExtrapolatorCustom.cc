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
// $Id: G4TablesForExtrapolatorCustom.cc 75145 2013-10-28 18:34:48Z vnivanch $
//
//---------------------------------------------------------------------------
//
// ClassName:    G4TablesForExtrapolatorCustom
//  
// Description:  This class provide calculation of energy loss, fluctuation, 
//               and msc angle
//
// Author:       09.12.04 V.Ivanchenko 
//
// Modification: 
// 08-04-05 Rename Propogator -> Extrapolator (V.Ivanchenko)
// 16-03-06 Add muon tables and fix bug in units (V.Ivanchenko)
// 21-03-06 Add verbosity defined in the constructor and Initialisation
//          start only when first public method is called (V.Ivanchenko)
// 03-05-06 Remove unused pointer G4Material* from number of methods (VI)
// 12-05-06 SEt linLossLimit=0.001 (VI)
//
//----------------------------------------------------------------------------
//
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "TrackPropagation/Geant4e/interface/G4TablesForExtrapolatorCustom.h"
#include "TrackPropagation/Geant4e/interface/Geant4ePropagator.h"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicsLogVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4EmParameters.hh"
#include "G4MollerBhabhaModel.hh"
#include "G4BetheBlochModel.hh"
#include "G4MuBetheBlochModel.hh"
#include "G4eBremsstrahlungRelModel.hh"
#include "G4MuPairProductionModel.hh"
#include "G4MuBremsstrahlungModel.hh"
#include "G4ProductionCuts.hh"
#include "G4LossTableBuilder.hh"
#include "G4WentzelVIModel.hh"
#include "G4ProductionCutsTable.hh"
#include "G4RunManagerKernel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4TablesForExtrapolatorCustom::G4TablesForExtrapolatorCustom(G4int verb, G4int bins, 
                                                 G4double e1, G4double e2, G4bool iononly)
  : verbose(verb), nbins(bins), emin(e1), emax(e2), ionOnly(iononly)
{
  Initialisation();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4TablesForExtrapolatorCustom:: ~G4TablesForExtrapolatorCustom()
{
  for(G4int i=0; i<nmat; i++) {delete couples[i];}

  dedxElectron->clearAndDestroy();
  dedxPositron->clearAndDestroy();
  dedxProton->clearAndDestroy();
  dedxMuon->clearAndDestroy();
  rangeElectron->clearAndDestroy();
  rangePositron->clearAndDestroy();
  rangeProton->clearAndDestroy();
  rangeMuon->clearAndDestroy();
  invRangeElectron->clearAndDestroy();
  invRangePositron->clearAndDestroy();
  invRangeProton->clearAndDestroy();
  invRangeMuon->clearAndDestroy();
  mscElectron->clearAndDestroy();

  delete dedxElectron;
  delete dedxPositron;
  delete dedxProton;
  delete dedxMuon;
  delete rangeElectron;
  delete rangePositron;
  delete rangeProton;
  delete rangeMuon;
  delete invRangeElectron;
  delete invRangePositron;
  delete invRangeProton;
  delete invRangeMuon;
  delete mscElectron;
  delete pcuts;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4PhysicsTable* 
G4TablesForExtrapolatorCustom::GetPhysicsTable(ExtTableType type) const
{
  const G4PhysicsTable* table = nullptr;
  switch (type) 
    {
    case fDedxElectron:
      table = dedxElectron;
      break;
    case fDedxPositron:
      table = dedxPositron;
      break;
    case fDedxProton:
      table = dedxProton;
      break;
    case fDedxMuon:
      table = dedxMuon;
      break;
    case fRangeElectron:
      table = rangeElectron;
      break;
    case fRangePositron:
      table = rangePositron;
      break;
    case fRangeProton:
      table = rangeProton;
      break;
    case fRangeMuon:
      table = rangeMuon;
      break;
    case fInvRangeElectron:
      table = invRangeElectron;
      break;
    case fInvRangePositron:
      table = invRangePositron;
      break;
    case fInvRangeProton:
      table = invRangeProton;
      break;
    case fInvRangeMuon:
      table = invRangeMuon;
      break;
    case fMscElectron:
      table = mscElectron;
    }
  return table;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4TablesForExtrapolatorCustom::Initialisation()
{
  if(verbose>1) {
    G4cout << "### G4TablesForExtrapolatorCustom::Initialisation" << G4endl;
  }
  currentParticle = nullptr;

  electron = G4Electron::Electron();
  positron = G4Positron::Positron();
  proton   = G4Proton::Proton();
  muonPlus = G4MuonPlus::MuonPlus();
  muonMinus= G4MuonMinus::MuonMinus();

  mass = charge2 = 0.0;

  nmat = G4Material::GetNumberOfMaterials();
  const G4MaterialTable* mtable = G4Material::GetMaterialTable();
  pcuts = new G4ProductionCuts();

  couples.resize(nmat,0);
  cuts.resize(nmat,DBL_MAX);

  for(G4int i=0; i<nmat; ++i) {
    couples[i] = new G4MaterialCutsCouple((*mtable)[i],pcuts);  
  }

  splineFlag = G4EmParameters::Instance()->Spline();

  dedxElectron     = PrepareTable();
  dedxPositron     = PrepareTable();
  dedxMuon         = PrepareTable();
  dedxProton       = PrepareTable();
  rangeElectron    = PrepareTable();
  rangePositron    = PrepareTable();
  rangeMuon        = PrepareTable();
  rangeProton      = PrepareTable();
  invRangeElectron = PrepareTable();
  invRangePositron = PrepareTable();
  invRangeMuon     = PrepareTable();
  invRangeProton   = PrepareTable();
  mscElectron      = PrepareTable();

  G4LossTableBuilder builder; 

  if(verbose>1) {
    G4cout << "### G4TablesForExtrapolatorCustom Builds electron tables" 
	   << G4endl;
  }
  ComputeElectronDEDX(electron, dedxElectron);
  builder.BuildRangeTable(dedxElectron,rangeElectron);  
  builder.BuildInverseRangeTable(rangeElectron, invRangeElectron);  

  if(verbose>1) {
    G4cout << "### G4TablesForExtrapolatorCustom Builds positron tables" 
	   << G4endl;
  }
  ComputeElectronDEDX(positron, dedxPositron);
  builder.BuildRangeTable(dedxPositron, rangePositron);  
  builder.BuildInverseRangeTable(rangePositron, invRangePositron);  

  if(verbose>1) {
    G4cout << "### G4TablesForExtrapolatorCustom Builds muon tables" << G4endl;
  }
  ComputeMuonDEDX(muonPlus, dedxMuon);
  builder.BuildRangeTable(dedxMuon, rangeMuon);  
  builder.BuildInverseRangeTable(rangeMuon, invRangeMuon);  
  /*
  G4cout << "DEDX MUON" << G4endl
  G4cout << *dedxMuon << G4endl;
  G4cout << "RANGE MUON" << G4endl
  G4cout << *rangeMuon << G4endl;
  G4cout << "INVRANGE MUON" << G4endl
  G4cout << *invRangeMuon << G4endl;
  */
  if(verbose>1) {
    G4cout << "### G4TablesForExtrapolatorCustom Builds proton tables" 
	   << G4endl;
  }
  ComputeProtonDEDX(proton, dedxProton);
  builder.BuildRangeTable(dedxProton, rangeProton);  
  builder.BuildInverseRangeTable(rangeProton, invRangeProton);  

  ComputeTrasportXS(electron, mscElectron);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsTable* G4TablesForExtrapolatorCustom::PrepareTable()
{
  G4PhysicsTable* table = new G4PhysicsTable();

  for(G4int i=0; i<nmat; ++i) {  

    G4PhysicsVector* v = new G4PhysicsLogVector(emin, emax, nbins);
    v->SetSpline(splineFlag);
    table->push_back(v);
  }
  return table;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4TablesForExtrapolatorCustom::ComputeElectronDEDX(
                                  const G4ParticleDefinition* part,
				  G4PhysicsTable* table) 
{
  G4MollerBhabhaModel* ioni = new G4MollerBhabhaModel();
  G4eBremsstrahlungRelModel* brem = new G4eBremsstrahlungRelModel();
  ioni->Initialise(part, cuts);
  brem->Initialise(part, cuts);

  mass    = electron_mass_c2;
  charge2 = 1.0;
  currentParticle = part;

  const G4MaterialTable* mtable = G4Material::GetMaterialTable();
  if(0<verbose) {
    G4cout << "G4TablesForExtrapolatorCustom::ComputeElectronDEDX for " 
           << part->GetParticleName() 
           << G4endl;
  }
  for(G4int i=0; i<nmat; ++i) {  

    const G4Material* mat = (*mtable)[i];
    if(1<verbose) {
      G4cout << "i= " << i << "  mat= " << mat->GetName() << G4endl;
    }
    const G4MaterialCutsCouple* couple = couples[i];
    G4PhysicsVector* aVector = (*table)[i];

    for(G4int j=0; j<=nbins; ++j) {
        
       G4double e = aVector->Energy(j);
       G4double dedx = ioni->ComputeDEDX(couple,part,e,e) 
	 + brem->ComputeDEDX(couple,part,e,e);
       if(1<verbose) {
         G4cout << "j= " << j
                << "  e(MeV)= " << e/MeV  
                << " dedx(Mev/cm)= " << dedx*cm/MeV
                << " dedx(Mev.cm2/g)= " 
		<< dedx/((MeV*mat->GetDensity())/(g/cm2)) 
		<< G4endl;
       }
       aVector->PutValue(j,dedx);
    }
    if(splineFlag) { aVector->FillSecondDerivatives(); }
  }
  delete ioni;
  delete brem;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4TablesForExtrapolatorCustom::ComputeMuonDEDX(const G4ParticleDefinition* part, 
					 G4PhysicsTable* table)
{
  G4BetheBlochModel* ioni = new G4BetheBlochModel();
  G4MuBetheBlochModel* ionialt = new G4MuBetheBlochModel();
//   G4MuBetheBlochModel* ioni = new G4MuBetheBlochModel();
  G4MuPairProductionModel* pair = new G4MuPairProductionModel();
  G4MuBremsstrahlungModel* brem = new G4MuBremsstrahlungModel();
  ioni->Initialise(part, cuts);
  ionialt->Initialise(part, cuts);
  pair->Initialise(part, cuts);
  brem->Initialise(part, cuts);

  mass    = part->GetPDGMass();
  charge2 = 1.0;
  currentParticle = part;

  const G4MaterialTable* mtable = G4Material::GetMaterialTable();


//   G4ProductionCutsTable::GetProductionCutsTable()->UpdateCoupleTable(G4RunManagerKernel::GetRunManagerKernel()->GetCurrentWorld());

//   const G4ProductionCutsTable* theCoupleTable=
//     G4ProductionCutsTable::GetProductionCutsTable();

//   std::cout << "couple table size = " << theCoupleTable->GetTableSize() << std::endl;

//   const G4DataVector *theCuts =
//     static_cast<const G4DataVector*>(theCoupleTable->GetEnergyCutsVector(1));

  if(0<verbose) {
    G4cout << "G4TablesForExtrapolatorCustom::ComputeMuonDEDX for " 
	   << part->GetParticleName() 
           << G4endl;
  }
 
  const double mass = 0.1056583745;
  const double massev = mass*1e9;
  G4double eMass     = 0.51099906 / GeV;
  G4double massRatio = eMass / mass;
 
  for(G4int i=0; i<nmat; ++i) {  

    const G4Material* mat = (*mtable)[i];
    if(1<verbose) {
      G4cout << "i= " << i << "  mat= " << mat->GetName() << G4endl;
    }
    const G4MaterialCutsCouple* couple = couples[i];
    G4PhysicsVector* aVector = (*table)[i];

//     const size_t icouple = couple->GetIndex();
//     const G4double cut  = (*theCuts)[icouple];
//     const G4double cut  = 0.;

//     int icouple = -1;
//     G4double cut = 0.;
//     for (unsigned int j=0; j<theCoupleTable->GetTableSize(); ++j) {
//       if (theCoupleTable->GetMaterialCutsCouple(j)->GetMaterial() == mat) {
//         icouple = j;
// //         cut = (*theCuts)[j];
//         cut = theCoupleTable->GetMaterialCutsCouple(j)->GetProductionCuts()->GetProductionCut("e-");
//         break;
//       }
//     }
//
//     std::cout << "i= " << i << "  mat= " << mat->GetName() << " icouple = " << icouple << " cut = " << cut << std::endl;
//     std::cout << "i= " << i << "  mat= " << mat->GetName() << std::endl;
    
// //     std::cout << "material = " << mat->GetName() << std::endl;
// //     std::cout << "computing dedxmin" << std::endl;
//     double dedxionimin = std::numeric_limits<double>::max();
// //     double emin = 0.;
// //     double dedxioniminalt = std::numeric_limits<double>::max();
// //     double eminalt = 0.;
//     for(G4int j=0; j<=nbins; j++) {
//       G4double e = aVector->Energy(j);
//       const double dedxioni = ioni->ComputeDEDX(couple,part,e,e);
// //       const double dedxionialt = ionialt->ComputeDEDX(couple,part,e,e);
// //       std::cout << "e = " << e << " dedxioni = " << dedxioni << " dedxionialt = " << dedxionialt << std::endl;
//       dedxionimin = std::min(dedxionimin, dedxioni);
// //       if (dedxioni < dedxionimin) {
// //         dedxionimin = dedxioni;
// //         emin = e;
// //       }
// //       if (dedxionialt < dedxioniminalt) {
// //         dedxioniminalt = dedxionialt;
// //         eminalt = e;
// //       }
//     }
    
    
//     std::cout << mat->GetName() << " dedxmin = " << dedxionimin << " emin = " << emin << " dedxioniminalt = " << dedxioniminalt << " eminalt = " << eminalt << std::endl;

//     const double cut = 0.121694;
    
    G4double effZ, effA;
    Geant4ePropagator::CalculateEffectiveZandA(mat, effZ, effA);
    G4double I = 16.*pow(effZ,0.9);
    const double f2 = effZ <= 2. ? 0. : 2./effZ;
    const double f1 = 1. - f2;
    const double e2 = 10.*effZ*effZ;
    const double e1 = pow(I/pow(e2,f2),1./f1);
    const double r = 0.4;

    for(G4int j=0; j<=nbins; j++) {
        
       G4double e = aVector->Energy(j);
       
       
       G4double pgev = e / GeV;
       G4double Etot = sqrt(pgev*pgev + mass*mass);
       G4double beta = pgev/Etot;
       G4double gamma = Etot / mass;
       G4double eta       = beta * gamma;
       G4double etasq     = eta * eta;
       G4double F1        = 2 * eMass * etasq;
       G4double F2        = 1. + 2. * massRatio * gamma + massRatio * massRatio;
       G4double Emax      = 1.E+6 * F1 / F2;  // now in keV
        
        
       const double emaxev = Emax*1e3; // keV -> eV
       

       const double sigma1partial = f1*(log(2.*massev*beta*beta*gamma*gamma/e1) - beta*beta)/e1/(log(2.*massev*beta*beta*gamma*gamma/I) - beta*beta)*(1.-r);

       const double sigma2partial = f1*(log(2.*massev*beta*beta*gamma*gamma/e2) - beta*beta)/e2/(log(2.*massev*beta*beta*gamma*gamma/I) - beta*beta)*(1.-r);

       const double sigma3partial = emaxev/I/(emaxev+I)/log((emaxev+I)/I)*r;

       const double e3med = I/(1. - 0.5*emaxev/(emaxev + I));
       const double e3mean = I*(emaxev +I)*log((emaxev+I)/I)/emaxev;
       const double e3mode = I;

       const double emed = sigma1partial*e1 + sigma2partial*e2 + sigma3partial*e3med;
       const double emean = sigma1partial*e1 + sigma2partial*e2 + sigma3partial*e3mean;
       const double emode = sigma1partial*e1 + sigma2partial*e2 + sigma3partial*e3mode;

       const double dedxratio = emed/emean;
       const double moderatio = emode/emean;
       const double alpha = 0.996;
//        const double alpha = 0.999;
       
       const double ealpha = I/(1. - alpha*emaxev/(emaxev + I));
//        const double ealphamev = ealpha*1e-6;

//        G4double dedx = std::min(ioni->ComputeDEDX(couple,part,e,e), dedxionimin) +
// 	               pair->ComputeDEDX(couple,part,e,e) +
// 	               brem->ComputeDEDX(couple,part,e,e);
//        G4double dedx = ioni->ComputeDEDX(couple,part,e,cut) +
// 	               pair->ComputeDEDX(couple,part,e,e) +
// 	               brem->ComputeDEDX(couple,part,e,e);
//        const double cut = i == 18 ? 0.121694 : e; // energy cut only for silicon for now
//        const double cut = e;
//        const double cut = 100.;
//        const double cut = ealphamev;
       const double dedxioni = e > 1000 ? ionialt->ComputeDEDX(couple,part,e,e) : ioni->ComputeDEDX(couple,part,e,e);
//        const double dedxioni = ioni->ComputeDEDX(couple,part,e,cut);
//        const double dedxioni = ioni->ComputeDEDX(couple,part,e,e);
//        const double dedxioni = dedxratio*ioni->ComputeDEDX(couple,part,e,e);
//        const double dedxioni = ioni->ComputeDEDX(couple,part,e,ealpha);
       G4double dedx = ionOnly ? dedxioni : dedxioni + pair->ComputeDEDX(couple,part,e,e) +
	               brem->ComputeDEDX(couple,part,e,e);
//        G4double dedx = dedxioni;
//        const double dedxionicut = ioni->ComputeDEDX(couple,part,e,cut);
//        const double dedxioni = ioni->ComputeDEDX(couple,part,e,e);
//        G4double dedx = ioni->ComputeDEDX(couple,part,e,e) +
// 	               pair->ComputeDEDX(couple,part,e,e) +
// 	               brem->ComputeDEDX(couple,part,e,e);
//         std::cout << "material = " << mat->GetName() << " e = " << e << " cut = " << cut << " dedxioni = " << dedxioni << " dedx = " << dedx << std::endl;
        std::cout << "material = " << mat->GetName() << " e = " << e << " dedxratio = " << dedxratio << " moderatio = " << moderatio << " dedxioni = " << dedxioni << " dedx = " << dedx << std::endl;


//        std::cout << "e = " << e << " dedxionicut = " << dedxionicut << " dedxioni = " << dedxioni << " dedx = " << dedx << std::endl;

//        dedx *= 1e-6;
//        dedx *= 1e-12;
       aVector->PutValue(j,dedx);
       if(1<verbose) {
         G4cout << "j= " << j
                << "  e(MeV)= " << e/MeV  
                << " dedx(Mev/cm)= " << dedx*cm/MeV
                << " dedx(Mev/(g/cm2)= " 
                << dedx/((MeV*mat->GetDensity())/(g/cm2))
		<< G4endl;
       }
    }
    if(splineFlag) { aVector->FillSecondDerivatives(); }
  }
  delete ioni;
  delete ionialt;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4TablesForExtrapolatorCustom::ComputeProtonDEDX(const G4ParticleDefinition* part, 
                        		   G4PhysicsTable* table)
{
  G4BetheBlochModel* ioni = new G4BetheBlochModel();
  ioni->Initialise(part, cuts);

  mass    = part->GetPDGMass();
  charge2 = 1.0;
  currentParticle = part;

  const G4MaterialTable* mtable = G4Material::GetMaterialTable();

  if(0<verbose) {
    G4cout << "G4TablesForExtrapolatorCustom::ComputeProtonDEDX for " 
	   << part->GetParticleName() 
           << G4endl;
  }
 
  for(G4int i=0; i<nmat; ++i) {  

    const G4Material* mat = (*mtable)[i];
    if(1<verbose)
      G4cout << "i= " << i << "  mat= " << mat->GetName() << G4endl;
    const G4MaterialCutsCouple* couple = couples[i];
    G4PhysicsVector* aVector = (*table)[i];
    for(G4int j=0; j<=nbins; j++) {
        
       G4double e = aVector->Energy(j);
       G4double dedx = ioni->ComputeDEDX(couple,part,e,e);
       aVector->PutValue(j,dedx);
       if(1<verbose) {
         G4cout << "j= " << j
                << "  e(MeV)= " << e/MeV  
                << " dedx(Mev/cm)= " << dedx*cm/MeV
                << " dedx(Mev.cm2/g)= " << dedx/((mat->GetDensity())/(g/cm2)) 
		<< G4endl;
       }
    }
    if(splineFlag) { aVector->FillSecondDerivatives(); }
  }
  delete ioni;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4TablesForExtrapolatorCustom::ComputeTrasportXS(const G4ParticleDefinition* part, 
					       G4PhysicsTable* table)
{
  G4WentzelVIModel* msc = new G4WentzelVIModel();
  msc->SetPolarAngleLimit(CLHEP::pi);
  msc->Initialise(part, cuts);

  mass    = part->GetPDGMass();
  charge2 = 1.0;
  currentParticle = part;

  const G4MaterialTable* mtable = G4Material::GetMaterialTable();

  if(0<verbose) {
    G4cout << "G4TablesForExtrapolatorCustom::ComputeProtonDEDX for " 
	   << part->GetParticleName() 
           << G4endl;
  }
 
  for(G4int i=0; i<nmat; ++i) {  

    const G4Material* mat = (*mtable)[i];
    msc->SetCurrentCouple(couples[i]);
    if(1<verbose)
      G4cout << "i= " << i << "  mat= " << mat->GetName() << G4endl;
    G4PhysicsVector* aVector = (*table)[i];
    for(G4int j=0; j<=nbins; j++) {
        
       G4double e = aVector->Energy(j);
       G4double xs = msc->CrossSectionPerVolume(mat,part,e);
       aVector->PutValue(j,xs);
       if(1<verbose) {
         G4cout << "j= " << j << "  e(MeV)= " << e/MeV  
                << " xs(1/mm)= " << xs*mm << G4endl;
       }
    }
    if(splineFlag) { aVector->FillSecondDerivatives(); }
  }
  delete msc;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

