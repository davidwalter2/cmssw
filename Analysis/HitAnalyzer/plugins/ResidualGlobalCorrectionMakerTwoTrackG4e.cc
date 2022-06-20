#include "ResidualGlobalCorrectionMakerBase.h"
#include "MagneticFieldOffset.h"

// required for Transient Tracks
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
// required for vtx fitting
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/MultiTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/MultiTrackPointingKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/CombinedKinematicConstraint.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "Math/Vector4Dfwd.h"

#include "TrackPropagation/Geant4e/interface/Geant4ePropagator.h"

#include "FWCore/Common/interface/TriggerNames.h"


class ResidualGlobalCorrectionMakerTwoTrackG4e : public ResidualGlobalCorrectionMakerBase
{
public:
  explicit ResidualGlobalCorrectionMakerTwoTrackG4e(const edm::ParameterSet &);
  ~ResidualGlobalCorrectionMakerTwoTrackG4e() {}

//   static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:

  virtual void beginStream(edm::StreamID) override;
  virtual void produce(edm::Event &, const edm::EventSetup &) override;
  
  bool doVtxConstraint_;
  bool doMassConstraint_;
  double massConstraint_;
  double massConstraintWidth_;
  
  float Jpsi_d;
  float Jpsi_x;
  float Jpsi_y;
  float Jpsi_z;
  float Jpsi_pt;
  float Jpsi_eta;
  float Jpsi_phi;
  float Jpsi_mass;
  
  float Jpsi_sigmamass;
  
  float Muplus_pt;
  float Muplus_eta;
  float Muplus_phi;
  
  float Muminus_pt;
  float Muminus_eta;
  float Muminus_phi;
  
  float Jpsikin_x;
  float Jpsikin_y;
  float Jpsikin_z;
  float Jpsikin_pt;
  float Jpsikin_eta;
  float Jpsikin_phi;
  float Jpsikin_mass;
  
  float Mupluskin_pt;
  float Mupluskin_eta;
  float Mupluskin_phi;
  
  float Muminuskin_pt;
  float Muminuskin_eta;
  float Muminuskin_phi;
  
  float Jpsicons_d;
  float Jpsicons_x;
  float Jpsicons_y;
  float Jpsicons_z;
  float Jpsicons_pt;
  float Jpsicons_eta;
  float Jpsicons_phi;
  float Jpsicons_mass;
  
  float Mupluscons_pt;
  float Mupluscons_eta;
  float Mupluscons_phi;
  
  float Muminuscons_pt;
  float Muminuscons_eta;
  float Muminuscons_phi;
  
  float Jpsikincons_x;
  float Jpsikincons_y;
  float Jpsikincons_z;
  float Jpsikincons_pt;
  float Jpsikincons_eta;
  float Jpsikincons_phi;
  float Jpsikincons_mass;
  
  float Mupluskincons_pt;
  float Mupluskincons_eta;
  float Mupluskincons_phi;
  
  float Muminuskincons_pt;
  float Muminuskincons_eta;
  float Muminuskincons_phi;
  
  float Jpsigen_x;
  float Jpsigen_y;
  float Jpsigen_z;
  float Jpsigen_pt;
  float Jpsigen_eta;
  float Jpsigen_phi;
  float Jpsigen_mass;
  
  float Muplusgen_pt;
  float Muplusgen_eta;
  float Muplusgen_phi;
  
  float Muminusgen_pt;
  float Muminusgen_eta;
  float Muminusgen_phi;
  
  std::array<float, 3> Muplus_refParms;
  std::array<float, 3> Muminus_refParms;
  
  std::vector<float> Muplus_jacRef;
  std::vector<float> Muminus_jacRef;
  
  unsigned int Muplus_nhits;
  unsigned int Muplus_nvalid;
  unsigned int Muplus_nvalidpixel;
  unsigned int Muplus_nmatchedvalid;
  unsigned int Muplus_nambiguousmatchedvalid;
  
  unsigned int Muminus_nhits;
  unsigned int Muminus_nvalid;
  unsigned int Muminus_nvalidpixel;
  unsigned int Muminus_nmatchedvalid;
  unsigned int Muminus_nambiguousmatchedvalid;
  
  bool Muplus_isMuon;
  bool Muplus_muonLoose;
  bool Muplus_muonMedium;
  bool Muplus_muonTight;
  bool Muplus_muonIsPF;
  bool Muplus_muonIsTracker;
  bool Muplus_muonIsGlobal;
  bool Muplus_muonIsStandalone;
  bool Muplus_muonInnerTrackBest;

  bool Muminus_isMuon;
  bool Muminus_muonLoose;
  bool Muminus_muonMedium;
  bool Muminus_muonTight;
  bool Muminus_muonIsPF;
  bool Muminus_muonIsTracker;
  bool Muminus_muonIsGlobal;
  bool Muminus_muonIsStandalone;
  bool Muminus_muonInnerTrackBest;

  float edmval_cons0;
  int niter_cons0;

  float dmassconvval = 0.;
  float dinvmasssqconvval = 0.;
  
  float dmassconvval_cons0 = 0.;
  float dinvmasssqconvval_cons0 = 0.;
  
//   std::vector<float> hessv;
  

  
};


ResidualGlobalCorrectionMakerTwoTrackG4e::ResidualGlobalCorrectionMakerTwoTrackG4e(const edm::ParameterSet &iConfig) : ResidualGlobalCorrectionMakerBase(iConfig) 
{
  doVtxConstraint_ = iConfig.getParameter<bool>("doVtxConstraint");
  doMassConstraint_ = iConfig.getParameter<bool>("doMassConstraint");
  massConstraint_ = iConfig.getParameter<double>("massConstraint");
  massConstraintWidth_ = iConfig.getParameter<double>("massConstraintWidth");
}

void ResidualGlobalCorrectionMakerTwoTrackG4e::beginStream(edm::StreamID streamid)
{
  ResidualGlobalCorrectionMakerBase::beginStream(streamid);
  
  if (fillTrackTree_) {
    tree->Branch("Jpsi_d", &Jpsi_d);
    tree->Branch("Jpsi_x", &Jpsi_x);
    tree->Branch("Jpsi_y", &Jpsi_y);
    tree->Branch("Jpsi_z", &Jpsi_z);
    tree->Branch("Jpsi_pt", &Jpsi_pt);
    tree->Branch("Jpsi_eta", &Jpsi_eta);
    tree->Branch("Jpsi_phi", &Jpsi_phi);
    tree->Branch("Jpsi_mass", &Jpsi_mass);
    
    tree->Branch("Jpsi_sigmamass", &Jpsi_sigmamass);
    
    tree->Branch("Muplus_pt", &Muplus_pt);
    tree->Branch("Muplus_eta", &Muplus_eta);
    tree->Branch("Muplus_phi", &Muplus_phi);
    
    tree->Branch("Muminus_pt", &Muminus_pt);
    tree->Branch("Muminus_eta", &Muminus_eta);
    tree->Branch("Muminus_phi", &Muminus_phi);
    
    tree->Branch("Jpsikin_x", &Jpsikin_x);
    tree->Branch("Jpsikin_y", &Jpsikin_y);
    tree->Branch("Jpsikin_z", &Jpsikin_z);
    tree->Branch("Jpsikin_pt", &Jpsikin_pt);
    tree->Branch("Jpsikin_eta", &Jpsikin_eta);
    tree->Branch("Jpsikin_phi", &Jpsikin_phi);
    tree->Branch("Jpsikin_mass", &Jpsikin_mass);
    
    tree->Branch("Mupluskin_pt", &Mupluskin_pt);
    tree->Branch("Mupluskin_eta", &Mupluskin_eta);
    tree->Branch("Mupluskin_phi", &Mupluskin_phi);
    
    tree->Branch("Muminuskin_pt", &Muminuskin_pt);
    tree->Branch("Muminuskin_eta", &Muminuskin_eta);
    tree->Branch("Muminuskin_phi", &Muminuskin_phi);
    
    tree->Branch("Jpsicons_d", &Jpsicons_d);
    tree->Branch("Jpsicons_x", &Jpsicons_x);
    tree->Branch("Jpsicons_y", &Jpsicons_y);
    tree->Branch("Jpsicons_z", &Jpsicons_z);
    tree->Branch("Jpsicons_pt", &Jpsicons_pt);
    tree->Branch("Jpsicons_eta", &Jpsicons_eta);
    tree->Branch("Jpsicons_phi", &Jpsicons_phi);
    tree->Branch("Jpsicons_mass", &Jpsicons_mass);
    
    tree->Branch("Mupluscons_pt", &Mupluscons_pt);
    tree->Branch("Mupluscons_eta", &Mupluscons_eta);
    tree->Branch("Mupluscons_phi", &Mupluscons_phi);
    
    tree->Branch("Muminuscons_pt", &Muminuscons_pt);
    tree->Branch("Muminuscons_eta", &Muminuscons_eta);
    tree->Branch("Muminuscons_phi", &Muminuscons_phi);
    
    tree->Branch("Jpsikincons_x", &Jpsikincons_x);
    tree->Branch("Jpsikincons_y", &Jpsikincons_y);
    tree->Branch("Jpsikincons_z", &Jpsikincons_z);
    tree->Branch("Jpsikincons_pt", &Jpsikincons_pt);
    tree->Branch("Jpsikincons_eta", &Jpsikincons_eta);
    tree->Branch("Jpsikincons_phi", &Jpsikincons_phi);
    tree->Branch("Jpsikincons_mass", &Jpsikincons_mass);
    
    tree->Branch("Mupluskincons_pt", &Mupluskincons_pt);
    tree->Branch("Mupluskincons_eta", &Mupluskincons_eta);
    tree->Branch("Mupluskincons_phi", &Mupluskincons_phi);
    
    tree->Branch("Muminuskincons_pt", &Muminuskincons_pt);
    tree->Branch("Muminuskincons_eta", &Muminuskincons_eta);
    tree->Branch("Muminuskincons_phi", &Muminuskincons_phi);
    
    tree->Branch("Jpsigen_x", &Jpsigen_x);
    tree->Branch("Jpsigen_y", &Jpsigen_y);
    tree->Branch("Jpsigen_z", &Jpsigen_z);
    tree->Branch("Jpsigen_pt", &Jpsigen_pt);
    tree->Branch("Jpsigen_eta", &Jpsigen_eta);
    tree->Branch("Jpsigen_phi", &Jpsigen_phi);
    tree->Branch("Jpsigen_mass", &Jpsigen_mass);
    
    tree->Branch("Muplusgen_pt", &Muplusgen_pt);
    tree->Branch("Muplusgen_eta", &Muplusgen_eta);
    tree->Branch("Muplusgen_phi", &Muplusgen_phi);
    
    tree->Branch("Muminusgen_pt", &Muminusgen_pt);
    tree->Branch("Muminusgen_eta", &Muminusgen_eta);
    tree->Branch("Muminusgen_phi", &Muminusgen_phi);
    
    tree->Branch("Muplus_refParms", Muplus_refParms.data(), "Muplus_refParms[3]/F");
    tree->Branch("Muminus_refParms", Muminus_refParms.data(), "Muminus_refParms[3]/F");
    
    if (fillJac_) {
      tree->Branch("Muplus_jacRef", &Muplus_jacRef);
      tree->Branch("Muminus_jacRef", &Muminus_jacRef);
    }
    
    tree->Branch("Muplus_nhits", &Muplus_nhits);
    tree->Branch("Muplus_nvalid", &Muplus_nvalid);
    tree->Branch("Muplus_nvalidpixel", &Muplus_nvalidpixel);
    tree->Branch("Muplus_nmatchedvalid", &Muplus_nmatchedvalid);
    tree->Branch("Muplus_nambiguousmatchedvalid", &Muplus_nambiguousmatchedvalid);
    
    tree->Branch("Muminus_nhits", &Muminus_nhits);
    tree->Branch("Muminus_nvalid", &Muminus_nvalid);
    tree->Branch("Muminus_nvalidpixel", &Muminus_nvalidpixel);
    tree->Branch("Muminus_nmatchedvalid", &Muminus_nmatchedvalid);
    tree->Branch("Muminus_nambiguousmatchedvalid", &Muminus_nambiguousmatchedvalid);

    tree->Branch("Muplus_isMuon", &Muplus_isMuon);
    tree->Branch("Muplus_muonLoose", &Muplus_muonLoose);
    tree->Branch("Muplus_muonMedium", &Muplus_muonMedium);
    tree->Branch("Muplus_muonTight", &Muplus_muonTight);
    tree->Branch("Muplus_muonIsPF", &Muplus_muonIsPF);
    tree->Branch("Muplus_muonIsTracker", &Muplus_muonIsTracker);
    tree->Branch("Muplus_muonIsGlobal", &Muplus_muonIsGlobal);
    tree->Branch("Muplus_muonIsStandalone", &Muplus_muonIsStandalone);
    tree->Branch("Muplus_muonInnerTrackBest", &Muplus_muonInnerTrackBest);

    tree->Branch("Muminus_isMuon", &Muminus_isMuon);
    tree->Branch("Muminus_muonLoose", &Muminus_muonLoose);
    tree->Branch("Muminus_muonMedium", &Muminus_muonMedium);
    tree->Branch("Muminus_muonTight", &Muminus_muonTight);
    tree->Branch("Muminus_muonIsPF", &Muminus_muonIsPF);
    tree->Branch("Muminus_muonIsTracker", &Muminus_muonIsTracker);
    tree->Branch("Muminus_muonIsGlobal", &Muminus_muonIsGlobal);
    tree->Branch("Muminus_muonIsStandalone", &Muminus_muonIsStandalone);
    tree->Branch("Muminus_muonInnerTrackBest", &Muminus_muonInnerTrackBest);

    tree->Branch("edmval_cons0", &edmval_cons0);
    tree->Branch("niter_cons0", &niter_cons0);

    tree->Branch("dmassconvval", &dmassconvval);
    tree->Branch("dinvmasssqconvval", &dinvmasssqconvval);
    tree->Branch("dmassconvval_cons0", &dmassconvval_cons0);
    tree->Branch("dinvmasssqconvval_cons0", &dinvmasssqconvval_cons0);

//     tree->Branch("hessv", &hessv);
    
  }
}


// ------------ method called for each event  ------------
void ResidualGlobalCorrectionMakerTwoTrackG4e::produce(edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  
  const bool dogen = fitFromGenParms_;
 
//   const bool dolocalupdate = true;
  const bool dolocalupdate = false;

  
  using namespace edm;

  Handle<reco::TrackCollection> trackOrigH;
  iEvent.getByToken(inputTrackOrig_, trackOrigH);


  // loop over gen particles

  edm::ESHandle<GlobalTrackingGeometry> globalGeometry;
  iSetup.get<GlobalTrackingGeometryRecord>().get(globalGeometry);
    
  edm::ESHandle<TrackerTopology> trackerTopology;
  iSetup.get<TrackerTopologyRcd>().get(trackerTopology);
  
  
  edm::ESHandle<TransientTrackingRecHitBuilder> ttrh;
  iSetup.get<TransientRecHitRecord>().get("WithAngleAndTemplate",ttrh);
  
  ESHandle<Propagator> thePropagator;
  iSetup.get<TrackingComponentsRecord>().get("Geant4ePropagator", thePropagator);

  const MagneticField* field = thePropagator->magneticField();
  const Geant4ePropagator *g4prop = dynamic_cast<const Geant4ePropagator*>(thePropagator.product());
  
//   Handle<std::vector<reco::GenParticle>> genPartCollection;
  Handle<edm::View<reco::Candidate>> genPartCollection;
  Handle<GenEventInfoProduct> genEventInfo;
  Handle<std::vector<int>> genPartBarcodes;
  Handle<std::vector<PileupSummaryInfo>> pileupSummary;
  if (doGen_) {
    iEvent.getByToken(GenParticlesToken_, genPartCollection);
    iEvent.getByToken(genEventInfoToken_, genEventInfo);
    iEvent.getByToken(pileupSummaryToken_, pileupSummary);
  }
  
  std::vector<Handle<std::vector<PSimHit>>> simHits(inputSimHits_.size());
  edm::Handle<std::vector<SimTrack>> simTracks;
  if (doSim_) {
    iEvent.getByToken(genParticlesBarcodeToken_, genPartBarcodes);
    for (unsigned int isimhit = 0; isimhit<inputSimHits_.size(); ++isimhit) {
      iEvent.getByToken(inputSimHits_[isimhit], simHits[isimhit]);
    }
    iEvent.getByToken(inputSimTracks_, simTracks);
  }
  
  Handle<edm::View<reco::Muon> > muons;
  if (doMuons_) {
    iEvent.getByToken(inputMuons_, muons);
  }
  
  Handle<edm::TriggerResults> triggerResults;
  if (doTrigger_)  {
    iEvent.getByToken(inputTriggerResults_, triggerResults);
  }

  KFUpdator updator;
  TkClonerImpl const& cloner = static_cast<TkTransientTrackingRecHitBuilder const *>(ttrh.product())->cloner();
  
  edm::ESHandle<TransientTrackBuilder> TTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", TTBuilder);
  KinematicParticleFactoryFromTransientTrack pFactory;

  Handle<reco::BeamSpot> bsH;
  iEvent.getByToken(inputBs_, bsH);
  
  constexpr double mmu = 0.1056583745;
  constexpr double mmuerr = 0.0000000024;

  VectorXd gradfull;
  MatrixXd hessfull;
  
  VectorXd dxfull;
  MatrixXd dxdparms;
  VectorXd grad;
  MatrixXd hess;
  LDLT<MatrixXd> Cinvd;
//   MatrixXd covstate;
  Matrix<double, 6, 6> covrefmom;
//   FullPivLU<MatrixXd> Cinvd;
//   ColPivHouseholderQR<MatrixXd> Cinvd;
  
  std::array<MatrixXd, 2> jacarr;
  
  run = iEvent.run();
  lumi = iEvent.luminosityBlock();
  event = iEvent.id().event();

  genweight = 1.;
  if (doGen_) {
    genweight = genEventInfo->weight();
    
    Pileup_nPU = pileupSummary->front().getPU_NumInteractions();
    Pileup_nTrueInt = pileupSummary->front().getTrueNumInteractions();
  }
  
  // trigger bits
  if (doTrigger_) {
    auto const &triggerNames = iEvent.triggerNames(*triggerResults);
    
    if (triggerNames.parameterSetID() != triggerNamesId_) {
      //trigger menu changed, update list of trigger path idxs
      
      triggerIdxs_.clear();
      
      for (auto const &trigger : triggers_) {
        const std::string basename = trigger + "_v";
//         std::cout << "basename = " << basename << std::endl;
        std::size_t idx = triggerNames.size();
        for (std::size_t itrig = 0; itrig < triggerNames.size(); ++itrig) {
//           auto findres = triggerNames.triggerName(itrig).find(basename);
//           std::cout << triggerNames.triggerName(itrig) << " findres = " << findres << std::endl;
          if (triggerNames.triggerName(itrig).find(basename) == 0) {
            idx = itrig;
            break;
          }
        }
        triggerIdxs_.push_back(idx);
      }
      
//       for (std::size_t itrig = 0; itrig < triggerIdxs_.size(); ++itrig) {
//         std::cout << "itrig = " << itrig << ", idx = " << triggerIdxs_[itrig] << std::endl;
//       }
      
      triggerNamesId_ = triggerNames.parameterSetID();
    }
    
    
    
    
    
    // set trigger decision bits
    for (std::size_t itrig = 0; itrig < triggerIdxs_.size(); ++itrig) {
      const std::size_t idx = triggerIdxs_[itrig];
//       if (idx < triggerResults->size()) {
//         std::cout << "itrig = " << itrig << " idx = " << idx << " accept = " << triggerResults->accept(idx) << std::endl;
//       }
      triggerDecisions_[itrig] = idx < triggerResults->size() ? triggerResults->accept(idx) : false;
    }
    
    
//     for (unsigned int i = 0; i < triggerNames.size(); ++i) {
//       std::cout << i << " " << triggerNames.triggerName(i) << " accept = " << triggerResults->accept(i) << std::endl;
//     }
//     auto const idx = iEvent.triggerNames(*triggerResults).triggerIndex("HLT_Mu7p5_Track2_Jpsi");
//     auto const idx = triggerNames.triggerIndex("HLT_Mu7p5_Track3p5_Jpsi_v4");
//     auto const idx2 = triggerNames.triggerIndex("HLT_eawgawe");
//     std::cout << "trigger names size = " << triggerNames.size() << std::endl;
//     std::cout << "trigger index = " << idx << std::endl;
//     std::cout << "trigger index2 = " << idx2 << std::endl;
  
  }
  
  // loop over combinatorics of track pairs
  for (auto itrack = trackOrigH->begin(); itrack != trackOrigH->end(); ++itrack) {
    if (itrack->isLooper()) {
      continue;
    }

    const reco::TransientTrack itt = TTBuilder->build(*itrack);


    const reco::Muon *matchedmuon0 = nullptr;
    if (doMuons_) {
      for (auto const &muon : *muons) {
        if (muon.bestTrack()->algo() == itrack->algo()) {
          if ( (muon.bestTrack()->momentum() - itrack->momentum()).mag2() < 1e-3 ) {
            matchedmuon0 = &muon;
          }
        }
        else if (muon.innerTrack().isNonnull() && muon.innerTrack()->algo() == itrack->algo()) {
          if ( (muon.innerTrack()->momentum() - itrack->momentum()).mag2() < 1e-3 ) {
            matchedmuon0 = &muon;
          }
        }
      }
    }

    for (auto jtrack = itrack + 1; jtrack != trackOrigH->end(); ++jtrack) {
      if (jtrack->isLooper()) {
        continue;
      }

      const reco::TransientTrack jtt = TTBuilder->build(*jtrack);

      const reco::Muon *matchedmuon1 = nullptr;
      if (doMuons_) {
        for (auto const &muon : *muons) {
          if (muon.bestTrack()->algo() == jtrack->algo()) {
            if ( (muon.bestTrack()->momentum() - jtrack->momentum()).mag2() < 1e-3 ) {
              matchedmuon1 = &muon;
            }
          }
          else if (muon.innerTrack().isNonnull() && muon.innerTrack()->algo() == jtrack->algo()) {
            if ( (muon.innerTrack()->momentum() - jtrack->momentum()).mag2() < 1e-3 ) {
              matchedmuon1 = &muon;
            }
          }
        }
      }

      const reco::Candidate *mu0gen = nullptr;
      const reco::Candidate *mu1gen = nullptr;
      
      double massconstraintval = massConstraint_;
      if (doGen_ && !doSim_) {
        for (auto const &genpart : *genPartCollection) {
          if (genpart.status() != 1) {
            continue;
          }
          if (std::abs(genpart.pdgId()) != 13) {
            continue;
          }
          
//           float dR0 = deltaR(genpart.phi(), itrack->phi(), genpart.eta(), itrack->eta());
          float dR0 = deltaR(genpart, *itrack);
          if (dR0 < 0.1 && genpart.charge() == itrack->charge()) {
            mu0gen = &genpart;
          }
          
//           float dR1 = deltaR(genpart.phi(), jtrack->phi(), genpart.eta(), jtrack->eta());
          float dR1 = deltaR(genpart, *jtrack);
          if (dR1 < 0.1 && genpart.charge() == jtrack->charge()) {
            mu1gen = &genpart;
          }
        }
        
//         if (mu0gen != nullptr && mu1gen != nullptr) {
//           auto const jpsigen = ROOT::Math::PtEtaPhiMVector(mu0gen->pt(), mu0gen->eta(), mu0gen->phi(), mmu) +
//                             ROOT::Math::PtEtaPhiMVector(mu1gen->pt(), mu1gen->eta(), mu1gen->phi(), mmu);
// 
//           massconstraintval = jpsigen.mass();
//         }
//         else {
//           continue;
//         }
        
      }
      
//       std::cout << "massconstraintval = " << massconstraintval << std::endl;
    
      std::array<TransientTrackingRecHit::RecHitContainer, 2> hitsarr;
      
      // prepare hits
      for (unsigned int id = 0; id < 2; ++id) {
        const reco::Track &track = id == 0 ? *itrack : *jtrack;
        auto &hits = hitsarr[id];
        hits.reserve(track.recHitsSize());
      
      
        for (auto it = track.recHitsBegin(); it != track.recHitsEnd(); ++it) {
          if ((*it)->geographicalId().det() != DetId::Tracker) {
            continue;
          }

          const GeomDet* detectorG = globalGeometry->idToDet((*it)->geographicalId());
          const GluedGeomDet* detglued = dynamic_cast<const GluedGeomDet*>(detectorG);
          
          // split matched invalid hits
          if (detglued != nullptr && !(*it)->isValid()) {
//             bool order = detglued->stereoDet()->surface().position().mag() > detglued->monoDet()->surface().position().mag();
            
            const auto stereopos = detglued->stereoDet()->surface().position();
            const auto monopos = detglued->monoDet()->surface().position();
            
            const Eigen::Vector3d stereoposv(stereopos.x(), stereopos.y(), stereopos.z());
            const Eigen::Vector3d monoposv(monopos.x(), monopos.y(), monopos.z());
            const Eigen::Vector3d trackmomv(track.momentum().x(), track.momentum().y(), track.momentum().z());
            
            bool order = (stereoposv - monoposv).dot(trackmomv) > 0.;
            
            const GeomDetUnit* detinner = order ? detglued->monoDet() : detglued->stereoDet();
            const GeomDetUnit* detouter = order ? detglued->stereoDet() : detglued->monoDet();
            
            hits.push_back(TrackingRecHit::RecHitPointer(new InvalidTrackingRecHit(*detinner, (*it)->type())));
            hits.push_back(TrackingRecHit::RecHitPointer(new InvalidTrackingRecHit(*detouter, (*it)->type())));
          }
          else {
            // apply hit quality criteria
            const bool ispixel = GeomDetEnumerators::isTrackerPixel(detectorG->subDetector());
            bool hitquality = false;
            if (applyHitQuality_ && (*it)->isValid()) {
              const TrackerSingleRecHit* tkhit = dynamic_cast<const TrackerSingleRecHit*>(*it);
              assert(tkhit != nullptr);
              
              if (ispixel) {
                const SiPixelRecHit *pixhit = dynamic_cast<const SiPixelRecHit*>(tkhit);
                const SiPixelCluster& cluster = *tkhit->cluster_pixel();
                assert(pixhit != nullptr);
                
//                 std::cout << "getSplitClusterErrorX = " << cluster.getSplitClusterErrorX() << std::endl;
                
//                 const double jpsieta = dimu_kinfit->currentState().freeTrajectoryState().momentum().eta();
//                 const double jpsipt = dimu_kinfit->currentState().freeTrajectoryState().momentum().perp();
//                 if (std::abs(jpsieta)>2. && jpsipt>20. && it == track.recHitsBegin()) {
//                   std::cout << "id = " << id << " detid = " << (*it)->geographicalId().rawId() << " minPixelRow = " << cluster.minPixelRow() << " maxPixelRow = " << cluster.maxPixelRow() << " minPixelCol = " << cluster.minPixelCol() << " maxPixelCol = " << cluster.maxPixelCol() << std::endl;
//                 }
                
                hitquality = !pixhit->isOnEdge() && cluster.sizeX() > 1;
//                 hitquality = !pixhit->isOnEdge() && cluster.sizeX() > 1 && pixhit->qBin() < 2;
//                 hitquality = !pixhit->isOnEdge() && cluster.sizeX() > 1 && cluster.sizeY() > 1;
              }
              else {
                assert(tkhit->cluster_strip().isNonnull());
                const SiStripCluster& cluster = *tkhit->cluster_strip();
                const StripTopology* striptopology = dynamic_cast<const StripTopology*>(&(detectorG->topology()));
                assert(striptopology);
                
                const uint16_t firstStrip = cluster.firstStrip();
                const uint16_t lastStrip = cluster.firstStrip() + cluster.amplitudes().size() - 1;
                const bool isOnEdge = firstStrip == 0 || lastStrip == (striptopology->nstrips() - 1);
                
    //             if (isOnEdge) {
    //               std::cout << "strip hit isOnEdge" << std::endl;
    //             }
                
//                 hitquality = !isOnEdge;
                hitquality = true;
              }
              
            }
            else {
              hitquality = true;
            }

            
            if (hitquality) {
              hits.push_back((*it)->cloneForFit(*detectorG));              
            }
            else {
              hits.push_back(TrackingRecHit::RecHitPointer(new InvalidTrackingRecHit(*detectorG, TrackingRecHit::inactive)));
            }
          }          
        }
      }
      
      unsigned int nhits = 0;
      unsigned int nvalid = 0;
      unsigned int nvalidpixel = 0;
      unsigned int nvalidalign2d = 0;
      
      
      std::array<unsigned int, 2> nhitsarr = {{ 0, 0 }};
      std::array<unsigned int, 2> nvalidarr = {{ 0, 0 }};
      std::array<unsigned int, 2> nvalidpixelarr = {{ 0, 0 }};
      std::array<unsigned int, 2> nmatchedvalidarr = {{ 0, 0 }};
      std::array<unsigned int, 2> nambiguousmatchedvalidarr = {{ 0, 0 }};
      
      // second loop to count hits
      for (unsigned int id = 0; id < 2; ++id) {
        std::unordered_map<unsigned int, unsigned int> trackidmap;
        auto const &hits = hitsarr[id];
//         layerStatesarr[id].reserve(hits.size());
        for (auto const &hit : hits) {
          ++nhits;
          ++nhitsarr[id];
          if (hit->isValid()) {
            ++nvalid;
            ++nvalidarr[id];
            
            const bool ispixel = GeomDetEnumerators::isTrackerPixel(hit->det()->subDetector());
            if (ispixel) {
              ++nvalidpixel;
              ++nvalidpixelarr[id];
            }
            
            
            const bool align2d = detidparms.count(std::make_pair(1, hit->geographicalId()));
            if (align2d) {
              ++nvalidalign2d;
            }
            
//             std::cout << hit->localPosition() << std::endl;
            
//             std::cout << "one hit:" << std:: endl;
            
            // count matching hits from sim tracks
            std::unordered_set<unsigned int> trackidset;
            if (doSim_) {
              for (auto const& simhith : simHits) {
                for (const PSimHit& simHit : *simhith) {
//                   if (simHit.detUnitId() == hit->geographicalId()) {
// //                     std::cout << "trackId = " << simHit.trackId() << " particleType = " << simHit.particleType() << "localpos = " << simHit.localPosition() << std::endl;
//                   }
                  
//                   if (simHit.detUnitId() == hit->geographicalId()) {
                  if (simHit.detUnitId() == hit->geographicalId() && std::abs(simHit.particleType()) == 13) {
                    //only count each simtrack once on a given detid
                    if (trackidset.count(simHit.trackId())) {
                      continue;
                    }
                    else {
                      trackidset.insert(simHit.trackId());
                    }
                    if (trackidmap.count(simHit.trackId())) {
                      ++trackidmap[simHit.trackId()];
                    }
                    else {
                      trackidmap[simHit.trackId()] = 1;
                    }
                  }
                }
              }
              if (trackidset.size() > 1) {
                ++nambiguousmatchedvalidarr[id];
              }
              
            }
            
          } 
        }
        
        // check for 50%+1 match
        for (auto const &pair : trackidmap) {
          const unsigned int trackid = pair.first;
          const unsigned int nmatch = pair.second;
          if (nmatch > nvalidarr[id]/2) {
            nmatchedvalidarr[id] = nmatch;
            //50% + 1 match found, find the sim track
            for (auto const& simTrack : *simTracks) {
              if (simTrack.trackId() == trackid) {
                // now find corresponding gen particle
                for (auto g = genPartCollection->begin(); g != genPartCollection->end(); ++g) {
                  const int genBarcode = (*genPartBarcodes)[g - genPartCollection->begin()];
                  if (genBarcode == simTrack.genpartIndex()) {
                    if (id == 0) {
                      mu0gen = &(*g);
                    }
                    else if (id == 1) {
                      mu1gen = &(*g);
                    }
                    break;
                  }
                }
                break;
              }
            }
            break;
          }
        }
        
      }

      if (nhitsarr[0] == 0 || nhitsarr[1] == 0) {
        continue;
      }
      
//       if (mu0gen == nullptr || mu1gen == nullptr || mu0gen->eta()<2.2 || mu1gen->eta()<2.2) {
//         continue;
//       }
      
      
      AlgebraicSymMatrix55 null55;
      const CurvilinearTrajectoryError nullerr(null55);

      
//       const unsigned int nparsAlignment = 2*nvalid + nvalidalign2d;
//       const unsigned int nparsAlignment = 6*nvalid;
      const unsigned int nparsAlignment = 5*nvalid + nvalidalign2d;
      const unsigned int nparsBfield = nhits;
      const unsigned int nparsEloss = nhits;
//       const unsigned int nparsEloss = nhits + 2;
      const unsigned int npars = nparsAlignment + nparsBfield + nparsEloss;
      
      const unsigned int nstateparms = 10 + 5*nhits;
      const unsigned int nparmsfull = nstateparms + npars;
      
    
      
      bool valid = true;
      
      
      const unsigned int nicons = doMassConstraint_ ? 2 : 1;
//       const unsigned int nicons = doMassConstraint_ ? 3 : 1;
      
      for (unsigned int icons = 0; icons < nicons; ++icons) {
        

        // common vertex fit
        std::vector<RefCountedKinematicParticle> parts;
        
        float masserr = mmuerr;
        float chisq = 0.;
        float ndf = 0.;
        parts.push_back(pFactory.particle(itt, mmu, chisq, ndf, masserr));
        parts.push_back(pFactory.particle(jtt, mmu, chisq, ndf, masserr));
        
        RefCountedKinematicTree kinTree;
        if (icons > 0) {
//           double kinconstraintval = 1./std::sqrt(massconstraintval);
          double kinconstraintval = massconstraintval;
          TwoTrackMassKinematicConstraint constraint(kinconstraintval);
          KinematicConstrainedVertexFitter vtxFitter;
          kinTree = vtxFitter.fit(parts, &constraint);
        }
        else {
          KinematicParticleVertexFitter vtxFitter;
          kinTree = vtxFitter.fit(parts);
        }
        
        if (kinTree->isEmpty() || !kinTree->isConsistent()) {
//           continue;
          valid = false;
          break;
        }
        
        kinTree->movePointerToTheTop();
        RefCountedKinematicParticle dimu_kinfit = kinTree->currentParticle();
        const double m0 = dimu_kinfit->currentState().mass();
        
        if (false) {
          // debug output
  //         kinTree->movePointerToTheTop();
          
  //         RefCountedKinematicParticle dimu_kinfit = kinTree->currentParticle();
          RefCountedKinematicVertex dimu_vertex = kinTree->currentDecayVertex();
          
          std::cout << dimu_kinfit->currentState().mass() << std::endl;
          std::cout << dimu_vertex->position() << std::endl;
        }
        
        const std::vector<RefCountedKinematicParticle> outparts = kinTree->finalStateParticles();
//         std::array<Matrix<double, 7, 1>, 2> refftsarr = {{ outparts[0]->currentState().freeTrajectoryState(),
//                                                           outparts[1]->currentState().freeTrajectoryState() }};
              
        std::array<Matrix<double, 7, 1>, 2> refftsarr;
        
        if (fitFromGenParms_ && mu0gen != nullptr && mu1gen != nullptr) {
          refftsarr[0][0] = mu0gen->vertex().x();
          refftsarr[0][1] = mu0gen->vertex().y();
          refftsarr[0][2] = mu0gen->vertex().z();
          refftsarr[0][3] = mu0gen->momentum().x();
          refftsarr[0][4] = mu0gen->momentum().y();
          refftsarr[0][5] = mu0gen->momentum().z();
          refftsarr[0][6] = mu0gen->charge();
          
          refftsarr[1][0] = mu1gen->vertex().x();
          refftsarr[1][1] = mu1gen->vertex().y();
          refftsarr[1][2] = mu1gen->vertex().z();
          refftsarr[1][3] = mu1gen->momentum().x();
          refftsarr[1][4] = mu1gen->momentum().y();
          refftsarr[1][5] = mu1gen->momentum().z();
          refftsarr[1][6] = mu1gen->charge();
        }
        else {
          for (unsigned int id = 0; id < 2; ++id) {
            const GlobalPoint pos = outparts[id]->currentState().freeTrajectoryState().position();
            const GlobalVector mom = outparts[id]->currentState().freeTrajectoryState().momentum();
            
            refftsarr[id][0] = pos.x();
            refftsarr[id][1] = pos.y();
            refftsarr[id][2] = pos.z();
            refftsarr[id][3] = mom.x();
            refftsarr[id][4] = mom.y();
            refftsarr[id][5] = mom.z();
            refftsarr[id][6] = outparts[id]->currentState().freeTrajectoryState().charge();
          }
          
        }
        
        std::array<std::vector<Matrix<double, 7, 1>>, 2> layerStatesarr;
        for (unsigned int id = 0; id < 2; ++id) {
          auto const &hits = hitsarr[id];
          layerStatesarr[id].reserve(hits.size());
        }
        

        double chisqvalold = std::numeric_limits<double>::max();

        std::array<unsigned int, 2> trackstateidxarr;
        std::array<int, 2> muchargearr;
        
  //       constexpr unsigned int niters = 1;
//         constexpr unsigned int niters = 3;
//         constexpr unsigned int niters = 5;
//         constexpr unsigned int niters = 10;
        
//         constexpr unsigned int niters = 1;
        const unsigned int niters = (dogen && !dolocalupdate) ? 1 : 10;


//         const unsigned int niters = icons == 0 ? 10 : 1;
        

        for (unsigned int iiter=0; iiter<niters; ++iiter) {
          
          gradfull = VectorXd::Zero(nparmsfull);
          hessfull = MatrixXd::Zero(nparmsfull, nparmsfull);

          globalidxv.clear();
          globalidxv.resize(npars, 0);
          
//           nParms = npars;
//           if (fillTrackTree_) {
//             tree->SetBranchAddress("globalidxv", globalidxv.data());
//           }
          
          std::array<Matrix<double, 5, 7>, 2> FdFmrefarr;
//           std::array<unsigned int, 2> trackstateidxarr;
          std::array<unsigned int, 2> trackparmidxarr;
          
          unsigned int trackstateidx = 10;
          unsigned int parmidx = 0;
          unsigned int alignmentparmidx = 0;
          
          double chisq0val = 0.;
          
          if (iiter > 0) {
            //update current state from reference point state
            const Matrix<double, 10, 1> statepca = twoTrackCart2pca(refftsarr[0], refftsarr[1]);
            const Matrix<double, 10, 1> statepcaupd = statepca + dxfull.head<10>();

            refftsarr = twoTrackPca2cart(statepcaupd);
          }

          
  //         const bool firsthitshared = hitsarr[0][0]->sharesInput(&(*hitsarr[1][0]), TrackingRecHit::some);
          
  //         std::cout << "firsthitshared = " << firsthitshared << std::endl;
          
          double dbetavalref = 0.;
          std::array<double, 2> dbetavalrefarr;
          for (unsigned int id = 0; id < 2; ++id) {
              auto &hits = hitsarr[id];
              auto const& hit = hits[0];
              
              const uint32_t gluedid = trackerTopology->glued(hit->det()->geographicalId());
              const bool isglued = gluedid != 0;
              const DetId parmdetid = isglued ? DetId(gluedid) : hit->geographicalId();

              const unsigned int bfieldglobalidx = detidparms.at(std::make_pair(6, parmdetid));                
              
              const double dbetaval = corparms_[bfieldglobalidx];
              dbetavalref += 0.5*dbetaval;
              dbetavalrefarr[id] = dbetaval;
          }

          const Matrix<double, 10, 10> twotrackpca2curvref = twoTrackPca2curvJacobianD(refftsarr[0], refftsarr[1], field, dbetavalrefarr[0], dbetavalrefarr[1]);

          
          for (unsigned int id = 0; id < 2; ++id) {
    //         FreeTrajectoryState refFts = outparts[id]->currentState().freeTrajectoryState();
//             FreeTrajectoryState &refFts = refftsarr[id];
            
            Matrix<double, 7, 1> &refFts = refftsarr[id];
            auto &hits = hitsarr[id];
            
            std::vector<Matrix<double, 7, 1>> &layerStates = layerStatesarr[id];
                      
            trackstateidxarr[id] = trackstateidx;
            trackparmidxarr[id] = parmidx;
            
            const unsigned int tracknhits = hits.size();
            
//             if (iiter > 0 || icons > 0) {
//             if (iiter > 0) {
//               //update current state from reference point state (errors not needed beyond first iteration)
//
//               const double qbp = refFts[6]/refFts.segment<3>(3).norm();
//               const double lam = std::atan(refFts[5]/std::sqrt(refFts[3]*refFts[3] + refFts[4]*refFts[4]));
//               const double phi = std::atan2(refFts[4], refFts[3]);
//
//               const double qbpupd = qbp + dxfull(trackstateidx);
//               const double lamupd = lam + dxfull(trackstateidx + 1);
//               const double phiupd = phi + dxfull(trackstateidx + 2);
//
//               const double charge = std::copysign(1., qbpupd);
//               const double pupd = std::abs(1./qbpupd);
//
//               const double pxupd = pupd*std::cos(lamupd)*std::cos(phiupd);
//               const double pyupd = pupd*std::cos(lamupd)*std::sin(phiupd);
//               const double pzupd = pupd*std::sin(lamupd);
//
//               refFts.head<3>() += dxfull.head<3>();
//               refFts[3] = pxupd;
//               refFts[4] = pyupd;
//               refFts[5] = pzupd;
//               refFts[6] = charge;
//             }
            
//             // initialize with zero uncertainty
//             refFts = FreeTrajectoryState(refFts.parameters(), nullerr);
//             
//             const ROOT::Math::PxPyPzMVector momtmp(refFts[3], refFts[4], refFts[5], mmu);
//             const Matrix<double, 5, 1> Felossadhoc = elossAdHocJacobianD(refFts, mmu);
//             const unsigned int etaphiidx = hetaphi->FindFixBin(momtmp.eta(), momtmp.phi());
//             
//   //           std::cout << "refFts:" << std::endl;
//   //           std::cout << refFts.position() << std::endl;
//   //           std::cout << refFts.momentum() << std::endl;
//   //           std::cout << refFts.charge() << std::endl;
//             
//   //           auto const &surface0 = *hits[0]->surface();
//             auto const &surface0 = *surfacemap_.at(hits[0]->geographicalId());
//   //           auto propresult = fPropagator->geometricalPropagator().propagateWithPath(refFts, surface0);
//       //       auto propresult = fPropagator->geometricalPropagator().propagateWithPath(refFts, *beampipe);
// //             auto const propresultref = g4prop->propagateGenericWithJacobian(refFts, surface0);
//             auto const propresultref = g4prop->propagateGenericWithJacobianAlt(refFts, surface0);
//             if (!std::get<0>(propresultref).isValid()) {
//               std::cout << "Abort: Propagation of reference state Failed!" << std::endl;
//               valid = false;
//               break;
//             }
//             TrajectoryStateOnSurface updtsos = std::get<0>(propresultref);
//             
//             const Matrix<double, 5, 6> hybrid2curvref = hybrid2curvJacobian(refFts);

//             
//   //           JacobianCartesianToCurvilinear cart2curvref(refFts.parameters());
//   //           auto const &jacCart2CurvRef = Map<const Matrix<double, 5, 6, RowMajor>>(cart2curvref.jacobian().Array());
//             
//             Matrix<double, 5, 7> FdFm = Map<const Matrix<double, 5, 7, RowMajor>>(std::get<1>(propresultref).Array());
//             
//             FdFmrefarr[id] = FdFm;
//             
//             double dEdxlast = std::get<3>(propresultref);

            

//             const Matrix<double, 5, 6> hybrid2curvref = hybrid2curvJacobianD(refFts, field, dbetavalref);

//             std::cout << "id = " << id << " twotrackpca2curvref:\n" << twotrackpca2curvref << "\nhybrid2curvref:\n" << hybrid2curvref << std::endl;

            Matrix<double, 7, 1> updtsos = refFts;
            
            
            if (bsConstraint_) {
              // apply beamspot constraint
              // TODO add residual corrections for beamspot parameters?
              // TODO impose constraints on individual tracks when not applying common vertex constraint?
              
              constexpr unsigned int nlocalvtx = 3;
              
              constexpr unsigned int nlocal = nlocalvtx;
              
              constexpr unsigned int localvtxidx = 0;
              
              constexpr unsigned int fullvtxidx = 7;
              
              using BSScalar = AANT<double, nlocal>;
              
              const double sigb1 = bsH->BeamWidthX();
              const double sigb2 = bsH->BeamWidthY();
              const double sigb3 = bsH->sigmaZ();
              const double dxdz = bsH->dxdz();
              const double dydz = bsH->dydz();
              const double x0 = bsH->x0();
              const double y0 = bsH->y0();
              const double z0 = bsH->z0();
              
              
              // covariance matrix of luminous region in global coordinates
              // taken from https://github.com/cms-sw/cmssw/blob/abc1f17b230effd629c9565fb5d95e527abcb294/RecoVertex/BeamSpotProducer/src/FcnBeamSpotFitPV.cc#L63-L90

              // FIXME xy correlation is not stored and assumed to be zero
              const double corrb12 = 0.;
              
              const double varb1 = sigb1*sigb1;
              const double varb2 = sigb2*sigb2;
              const double varb3 = sigb3*sigb3;
              
              Matrix<double, 3, 3> covBS = Matrix<double, 3, 3>::Zero();
              // parametrisation: rotation (dx/dz, dy/dz); covxy
              covBS(0,0) = varb1;
              covBS(1,0) = covBS(0,1) = corrb12*sigb1*sigb2;
              covBS(1,1) = varb2;
              covBS(2,0) = covBS(0,2) = dxdz*(varb3-varb1)-dydz*covBS(1,0);
              covBS(2,1) = covBS(1,2) = dydz*(varb3-varb2)-dxdz*covBS(1,0);
              covBS(2,2) = varb3;

      //         std::cout << "covBS:" << std::endl;
      //         std::cout << covBS << std::endl;
              
              Matrix<BSScalar, 3, 1> dvtx = Matrix<BSScalar, 3, 1>::Zero();
              for (unsigned int j=0; j<dvtx.size(); ++j) {
                init_twice_active_var(dvtx[j], nlocal, localvtxidx + j);
              }
              
              Matrix<BSScalar, 3, 1> dbs0;
              dbs0[0] = BSScalar(refFts[0] - x0);
              dbs0[1] = BSScalar(refFts[1] - y0);
              dbs0[2] = BSScalar(refFts[2] - z0);
              
      //         std::cout << "dposition / d(qop, lambda, phi) (should be 0?):" << std::endl;
      //         std::cout << Map<const Matrix<double, 6, 5, RowMajor>>(jac.Array()).topLeftCorner<3,3>() << std::endl;
              
              const Matrix<BSScalar, 3, 3> covBSinv = covBS.inverse().cast<BSScalar>();
              
              const Matrix<BSScalar, 3, 1> dbs = dbs0 + dvtx;
              const BSScalar chisq = dbs.transpose()*covBSinv*dbs;
              
              chisq0val += chisq.value().value();
              
              auto const& gradlocal = chisq.value().derivatives();
              //fill local hessian
              Matrix<double, nlocal, nlocal> hesslocal;
              for (unsigned int j=0; j<nlocal; ++j) {
                hesslocal.row(j) = chisq.derivatives()[j].derivatives();
              }
              
              //fill global gradient
              gradfull.segment<nlocalvtx>(fullvtxidx) += gradlocal.head<nlocalvtx>();
              //fill global hessian (upper triangular blocks only)
              hessfull.block<nlocalvtx,nlocalvtx>(fullvtxidx, fullvtxidx) += hesslocal.topLeftCorner<nlocalvtx,nlocalvtx>();
              
            }
            
            

            for (unsigned int ihit = 0; ihit < hits.size(); ++ihit) {
      //         std::cout << "ihit " << ihit << std::endl;
              auto const& hit = hits[ihit];
              
              const uint32_t gluedid = trackerTopology->glued(hit->det()->geographicalId());
              const bool isglued = gluedid != 0;
              const DetId parmdetid = isglued ? DetId(gluedid) : hit->geographicalId();
              const GeomDet* parmDet = isglued ? globalGeometry->idToDet(parmdetid) : hit->det();
              const double xifraction = isglued ? hit->det()->surface().mediumProperties().xi()/parmDet->surface().mediumProperties().xi() : 1.;

              const unsigned int bfieldglobalidx = detidparms.at(std::make_pair(6, parmdetid));                
              const unsigned int elossglobalidx = detidparms.at(std::make_pair(7, parmdetid));
              
              const double dbetaval = corparms_[bfieldglobalidx];
              const double dxival = corparms_[elossglobalidx];
              
//               if (ihit > 0) {
//                 if (std::abs(updtsos.globalMomentum().eta()) > 4.0) {
//                   std::cout << "WARNING:  Invalid state!!!" << std::endl;
//                   valid = false;
//                   break;
//                 }
//                 
//                 
//       //        auto const &surfaceip1 = *hits[ihit+1]->surface();
//       //           auto const &surface = *hit->surface();
//       //           const Plane &surface = *hit->surface();
//                 auto const &surface = *surfacemap_.at(hit->geographicalId());
//       //           auto propresult = thePropagator->propagateWithPath(updtsos, surface);
// //                 auto const propresult = g4prop->propagateGenericWithJacobian(*updtsos.freeState(), surface);
//                 auto const propresult = g4prop->propagateGenericWithJacobianAlt(*updtsos.freeState(), surface);
//       //           propresult = fPropagator->geometricalPropagator().propagateWithPath(updtsos, *hits[ihit+1]->surface());
//                 if (!std::get<0>(propresult).isValid()) {
//                   std::cout << "Abort: Propagation Failed!" << std::endl;
//                   valid = false;
//                   break;
//                 }
//                 
//                 FdFm = Map<const Matrix<double, 5, 7, RowMajor>>(std::get<1>(propresult).Array());
//       //           FdFm = localTransportJacobian(updtsos, propresult, false);
//                 updtsos = std::get<0>(propresult);
//                 
//                 dEdxlast = std::get<3>(propresult);
//               }
              
              
              const GloballyPositioned<double> &surface = surfacemapD_.at(hit->geographicalId());
              
              const Point3DBase<double, GlobalTag> crosspostmp(updtsos[0], updtsos[1], updtsos[2]);
              const Vector3DBase<double, GlobalTag> crossmomtmp(updtsos[3], updtsos[4], updtsos[5]);
              
              if (surface.toLocal(crosspostmp).z() * surface.toLocal(crossmomtmp).z() > 0) {
                std::cout << "Abort: wrong propagation direction!\n";
                valid = false;
                break;
              }
              
              auto propresult = g4prop->propagateGenericWithJacobianAltD(updtsos, surface, dbetaval, dxival);
              if (!std::get<0>(propresult)) {
                std::cout << "Abort: Propagation Failed!" << std::endl;
                valid = false;
                break;
              }
              
//               auto propresultalt = g4prop->propagateGenericWithJacobianAltD(updtsos, surface, 10.);
//               const Matrix<double, 5, 5> Qcurvalt = std::get<2>(propresultalt);
              
              updtsos = std::get<1>(propresult);
              const Matrix<double, 5, 5> Qcurv = std::get<2>(propresult);
              const Matrix<double, 5, 7> FdFm = std::get<3>(propresult);
              const double dEdxlast = std::get<4>(propresult);

//               std::cout << "Qcurv" << std::endl;
//               std::cout << Qcurv << std::endl;
//               std::cout << "Qcurvalt" << std::endl;
//               std::cout << Qcurvalt << std::endl;

              
              // compute convolution correction in local coordinates (BEFORE material effects are applied)
      //         const Matrix<double, 2, 1> dxlocalconv = localPositionConvolution(updtsos);
              
              // curvilinear to local jacobian
//               JacobianCurvilinearToLocal curv2localm(updtsos.surface(), updtsos.localParameters(), *updtsos.magneticField());
//               const AlgebraicMatrix55& curv2localjacm = curv2localm.jacobian();
//               const Matrix<double, 5, 5> Hm = Map<const Matrix<double, 5, 5, RowMajor>>(curv2localjacm.Array()); 
//               const Matrix<double, 5, 5> Hm = curv2localJacobianAlt(updtsos);
//               const Matrix<double, 5, 5> Hm = curv2localJacobianAlteloss(updtsos, dEdxlast, mmu);
              
              const Matrix<double, 5, 5> Hm = curv2localJacobianAltelossD(updtsos, field, surface, dEdxlast, mmu, dbetaval);
              
              //get the process noise matrix
//               AlgebraicMatrix55 const Qmat = updtsos.localError().matrix();
//               const Map<const Matrix<double, 5, 5, RowMajor>>Q(Qmat.Array());
              
//               const Matrix<double, 5, 5> Q = Hm*Qcurv*Hm.transpose();
              const Matrix<double, 5, 5> Q = dolocalupdate ? Hm*Qcurv*Hm.transpose() : Qcurv;

              Matrix<double, 5, 1> localparmsprop;
              const Point3DBase<double, GlobalTag> posprop(updtsos[0], updtsos[1], updtsos[2]);
              const Vector3DBase<double, GlobalTag> momprop(updtsos[3], updtsos[4], updtsos[5]);
              
              const Point3DBase<double, LocalTag> localpos = surface.toLocal(posprop);
              const Vector3DBase<double, LocalTag> localmom = surface.toLocal(momprop);
              
              localparmsprop[0] = updtsos[6]/updtsos.segment<3>(3).norm();
              localparmsprop[1] = localmom.x()/localmom.z();
              localparmsprop[2] = localmom.y()/localmom.z();
              localparmsprop[3] = localpos.x();
              localparmsprop[4] = localpos.y();        
              
              Matrix<double, 5, 1> localparms = localparmsprop;

              Matrix<double, 5, 1> idx0 = Matrix<double, 5, 1>::Zero();
              if (dolocalupdate) {
                if (iiter==0) {
        //         if (true) {
                  layerStates.push_back(updtsos);
        //           layerStatesStart.push_back(updtsos);
                  

                }
                else {
                  //current state from previous state on this layer
                  //save current parameters  
                  
                  Matrix<double, 7, 1>& oldtsos = layerStates[ihit];
                  const Matrix<double, 5, 5> Hold = curv2localJacobianAltelossD(oldtsos, field, surface, dEdxlast, mmu, dbetaval);
  //                 const Matrix<double, 5, 1> dxlocal = Hold*dxfull.segment<5>(5*(ihit+1));
                  const Matrix<double, 5, 1> dxlocal = Hold*dxfull.segment<5>(trackstateidx + 3 + 5*ihit);
                  
                  const Point3DBase<double, GlobalTag> pos(oldtsos[0], oldtsos[1], oldtsos[2]);
                  const Point3DBase<double, LocalTag> localpos = surface.toLocal(pos);
                  
                  const Point3DBase<double, LocalTag> localposupd(localpos.x() + dxlocal[3], localpos.y() + dxlocal[4], localpos.z());
                  const Point3DBase<double, GlobalTag> posupd = surface.toGlobal(localposupd);
                
                  
                  const Vector3DBase<double, GlobalTag> mom(oldtsos[3], oldtsos[4], oldtsos[5]);
                  const Vector3DBase<double, LocalTag> localmom = surface.toLocal(mom);
                  
                  const double dxdz = localmom.x()/localmom.z();
                  const double dydz = localmom.y()/localmom.z();
                  
                  
                  
                  const double dxdzupd = dxdz + dxlocal[1];
                  const double dydzupd = dydz + dxlocal[2];
                  
                  const double qop = oldtsos[6]/oldtsos.segment<3>(3).norm();
                  const double qopupd = qop + dxlocal[0];
                  
                  const double pupd = std::abs(1./qopupd);
                  const double charge = std::copysign(1., qopupd);
                  
                  const double signpz = std::copysign(1., localmom.z());
                  const double localmomfact = signpz/std::sqrt(1. + dxdzupd*dxdzupd + dydzupd*dydzupd);
                  const Vector3DBase<double, LocalTag> localmomupd(pupd*dxdzupd*localmomfact, pupd*dydzupd*localmomfact, pupd*localmomfact);
                  const Vector3DBase<double, GlobalTag> momupd = surface.toGlobal(localmomupd);
                            
                  oldtsos[0] = posupd.x();
                  oldtsos[1] = posupd.y();
                  oldtsos[2] = posupd.z();
                  oldtsos[3] = momupd.x();
                  oldtsos[4] = momupd.y();
                  oldtsos[5] = momupd.z();
                  oldtsos[6] = charge;          
                  
                  updtsos = oldtsos;
                            
                  localparms[0] = qopupd;
                  localparms[1] = dxdzupd;
                  localparms[2] = dydzupd;
                  localparms[3] = localposupd.x();
                  localparms[4] = localposupd.y();
                  
                  idx0 = localparms - localparmsprop;
                  
                }
              }
              
              
//               // update state from previous iteration
//               //momentum kink residual
//               AlgebraicVector5 idx0(0., 0., 0., 0., 0.);
// //               if (iiter==0 && icons==0) {
//               if (iiter==0) {
//                 updtsos.update(updtsos.localParameters(),
//                                 LocalTrajectoryError(0.,0.,0.,0.,0.),
//                                 updtsos.surface(),
//                                 updtsos.magneticField(),
//                                 updtsos.surfaceSide());
//                 layerStates.push_back(updtsos);
//               }
//               else {          
//                 //current state from previous state on this layer
//                 //save current parameters          
//                 TrajectoryStateOnSurface& oldtsos = layerStates[ihit];
//                 
// //                 JacobianCurvilinearToLocal curv2localold(oldtsos.surface(), oldtsos.localParameters(), *oldtsos.magneticField());
// //                 const AlgebraicMatrix55& curv2localjacold = curv2localold.jacobian();
// //                 const Matrix<double, 5, 5> Hold = Map<const Matrix<double, 5, 5, RowMajor>>(curv2localjacold.Array());
// //                 const Matrix<double, 5, 5> Hold = curv2localJacobianAlt(oldtsos);
//                 const Matrix<double, 5, 5> Hold = curv2localJacobianAlteloss(oldtsos, dEdxlast, mmu);
// 
//                 
//                 const AlgebraicVector5 local = oldtsos.localParameters().vector();
//                 
//                 auto const& dxlocal = Hold*dxfull.segment<5>(trackstateidx + 3 + 5*ihit);
//                 const Matrix<double, 5, 1> localupd = Map<const Matrix<double, 5, 1>>(local.Array()) + dxlocal;
//                 AlgebraicVector5 localvecupd(localupd[0],localupd[1],localupd[2],localupd[3],localupd[4]);
//                 
//                 idx0 = localvecupd - updtsos.localParameters().vector();
//                 
//                 const LocalTrajectoryParameters localparms(localvecupd, oldtsos.localParameters().pzSign());
//                 
//       //           std::cout << "before update: oldtsos:" << std::endl;
//       //           std::cout << oldtsos.localParameters().vector() << std::endl;
//   //               oldtsos.update(localparms, oldtsos.surface(), field, oldtsos.surfaceSide());
//                 oldtsos.update(localparms, LocalTrajectoryError(0.,0.,0.,0.,0.), oldtsos.surface(), field, oldtsos.surfaceSide());
//       //           std::cout << "after update: oldtsos:" << std::endl;
//       //           std::cout << oldtsos.localParameters().vector() << std::endl;
//                 updtsos = oldtsos;
// 
//               }
              
              //apply measurement update if applicable
      //         std::cout << "constructing preciseHit" << std::endl;
//               auto const& preciseHit = hit->isValid() ? cloner.makeShared(hit, updtsos) : hit;

              const GlobalPoint pos(updtsos[0], updtsos[1], updtsos[2]);
              const GlobalVector mom(updtsos[3], updtsos[4], updtsos[5]);
              const GlobalTrajectoryParameters glob(pos, mom, updtsos[6], field);
              const TrajectoryStateOnSurface tsostmp(glob, *hit->surface());

              auto const& preciseHit = hit->isValid() ? cloner.makeShared(hit, tsostmp) : hit;

              if (hit->isValid() && !preciseHit->isValid()) {
                std::cout << "Abort: Failed updating hit" << std::endl;
                valid = false;
                break;
              }
              
      //         const uint32_t gluedid = trackerTopology->glued(preciseHit->det()->geographicalId());
      //         const bool isglued = gluedid != 0;
      //         const DetId parmdetid = isglued ? DetId(gluedid) : preciseHit->geographicalId();
      //         const bool align2d = detidparms.count(std::make_pair(1, parmdetid));
      //         const GeomDet* parmDet = isglued ? globalGeometry->idToDet(parmdetid) : preciseHit->det();
              
              const bool align2d = detidparms.count(std::make_pair(1, preciseHit->geographicalId()));
              
              // curvilinear to local jacobian
//               JacobianCurvilinearToLocal curv2localp(updtsos.surface(), updtsos.localParameters(), *updtsos.magneticField());
//               const AlgebraicMatrix55& curv2localjacp = curv2localp.jacobian();
//               const Matrix<double, 5, 5> Hp = Map<const Matrix<double, 5, 5, RowMajor>>(curv2localjacp.Array());
//               const Matrix<double, 5, 5> Hp = curv2localJacobianAlt(updtsos);
//               const Matrix<double, 5, 5> Hp = curv2localJacobianAlteloss(updtsos, dEdxlast, mmu);
              const Matrix<double, 5, 5> Hp = curv2localJacobianAltelossD(updtsos, field, surface, dEdxlast, mmu, dbetaval);
              
  //             const Matrix<double, 5, 5> Hpalt = curv2localJacobianAlt(updtsos);
  //             
  //             std::cout << "Hp" << std::endl;
  //             std::cout << Hp << std::endl;
  //             std::cout << "Hpalt" << std::endl;
  //             std::cout << Hpalt << std::endl;
              
              
              if (true) {                
    //             std::cout << "EdE first hit:" << std::endl;
    //             std::cout << EdE << std::endl;
    //             
    //             std::cout << "xival = " << xival << std::endl;
                
  //               AlgebraicVector5 idx0(0., 0., 0., 0., 0.);
//                 const Vector5d dx0 = Map<const Vector5d>(idx0.Array());
                const Vector5d dx0 = idx0;

                if (ihit == 0) {
                  constexpr unsigned int nvtxstate = 10;
                  constexpr unsigned int nlocalstate = 5;
                  constexpr unsigned int nlocalbfield = 1;
                  constexpr unsigned int nlocaleloss = 1;
                  constexpr unsigned int nlocalparms = nlocalbfield + nlocaleloss;
                  
                  constexpr unsigned int nlocal = nvtxstate + nlocalstate + nlocalparms;
                  
                  constexpr unsigned int localvtxidx = 0;
                  constexpr unsigned int localstateidx = localvtxidx + nvtxstate;
                  constexpr unsigned int localparmidx = localstateidx + nlocalstate;
                  
                  constexpr unsigned int fullvtxidx = 0;
                  const unsigned int fullstateidx = trackstateidx;
                  const unsigned int fullparmidx = nstateparms + parmidx;
                  
                  using MSScalar = AANT<double, nlocal>;
                                  
                  const unsigned int vtxjacidx = 5*id;

                  Matrix<MSScalar, 5, 10> Fvtx = (FdFm.leftCols<5>()*twotrackpca2curvref.middleRows<5>(vtxjacidx)).cast<MSScalar>();

                  Matrix<MSScalar, 5, 1> Fb = FdFm.col(5).cast<MSScalar>();
                  Matrix<MSScalar, 5, 1> Fxi = FdFm.col(6).cast<MSScalar>();
                  
                  Matrix<MSScalar, 5, 5> Hmstate = Hm.cast<MSScalar>();
                  Matrix<MSScalar, 5, 5> Hpstate = Hp.cast<MSScalar>();
                  
                  Matrix<MSScalar, 5, 5> Qinv = Q.inverse().cast<MSScalar>();
                  
                  // initialize active scalars for common vertex parameters
                  Matrix<MSScalar, 10, 1> dvtx = Matrix<MSScalar, 10, 1>::Zero();
                  for (unsigned int j=0; j<dvtx.size(); ++j) {
                    init_twice_active_var(dvtx[j], nlocal, localvtxidx + j);
                  }
                  
                  Matrix<MSScalar, 5, 1> du = Matrix<MSScalar, 5, 1>::Zero();
                  for (unsigned int j=0; j<du.size(); ++j) {
                    init_twice_active_var(du[j], nlocal, localstateidx + j);
                  }
                  
                  // initialize active scalars for correction parameters                  
//                   MSScalar dbeta(corparms_[bfieldglobalidx]);
                  MSScalar dbeta(0.);
                  init_twice_active_var(dbeta, nlocal, localparmidx);
                  
//                   MSScalar dxi(corparms_[elossglobalidx]);
                  MSScalar dxi(0.);
                  init_twice_active_var(dxi, nlocal, localparmidx + 1);
                  
//                   const Matrix<MSScalar, 5, 1> dprop = dx0.cast<MSScalar>() + Hpstate*du - Hmstate*Fvtx*dvtx - Hmstate*Fmom*dmom - Hmstate*Fb*dbeta - Hmstate*Fxi*dxi;
                  Matrix<MSScalar, 5, 1> dprop;
                  if (dolocalupdate) {
                    dprop = dx0.cast<MSScalar>() + Hpstate*du - Hmstate*Fvtx*dvtx - Hmstate*Fb*dbeta - Hmstate*Fxi*dxi;
                  }
                  else {
                    dprop = du - Fvtx*dvtx - Fb*dbeta - Fxi*dxi;
                  }
                  const MSScalar chisq = dprop.transpose()*Qinv*dprop;
                  
                  chisq0val += chisq.value().value();
                  
                  auto const& gradlocal = chisq.value().derivatives();
                  //fill local hessian
                  Matrix<double, nlocal, nlocal> hesslocal;
                  for (unsigned int j=0; j<nlocal; ++j) {
                    hesslocal.row(j) = chisq.derivatives()[j].derivatives();
                  }
                  
                  constexpr std::array<unsigned int, 3> localsizes = {{ nvtxstate, nlocalstate, nlocalparms }};
                  constexpr std::array<unsigned int, 3> localidxs = {{ localvtxidx, localstateidx, localparmidx }};
                  const std::array<unsigned int, 3> fullidxs = {{ fullvtxidx, fullstateidx, fullparmidx }};
                  
                  for (unsigned int iidx = 0; iidx < localidxs.size(); ++iidx) {
                    gradfull.segment(fullidxs[iidx], localsizes[iidx]) += gradlocal.segment(localidxs[iidx], localsizes[iidx]);
                    for (unsigned int jidx = 0; jidx < localidxs.size(); ++jidx) {
                      hessfull.block(fullidxs[iidx], fullidxs[jidx], localsizes[iidx], localsizes[jidx]) += hesslocal.block(localidxs[iidx], localidxs[jidx], localsizes[iidx], localsizes[jidx]);
                    }
                  }
                                  
                }
                else {
                  //TODO statejac stuff
                  
                  constexpr unsigned int nlocalstate = 10;
                  constexpr unsigned int nlocalbfield = 1;
                  constexpr unsigned int nlocaleloss = 1;
                  constexpr unsigned int nlocalparms = nlocalbfield + nlocaleloss;
                  
                  constexpr unsigned int nlocal = nlocalstate + nlocalparms;
                  
                  constexpr unsigned int localstateidx = 0;
                  constexpr unsigned int localparmidx = localstateidx + nlocalstate;
                  
                  const unsigned int fullstateidx = trackstateidx + 5*(ihit - 1);
                  const unsigned int fullparmidx = nstateparms + parmidx;
                  
                  using MSScalar = AANT<double, nlocal>;

                  Matrix<MSScalar, 5, 5> Fstate = FdFm.leftCols<5>().cast<MSScalar>();
                  Matrix<MSScalar, 5, 1> Fb = FdFm.col(5).cast<MSScalar>();
                  Matrix<MSScalar, 5, 1> Fxi = FdFm.col(6).cast<MSScalar>();
                  
                  Matrix<MSScalar, 5, 5> Hmstate = Hm.cast<MSScalar>();
                  Matrix<MSScalar, 5, 5> Hpstate = Hp.cast<MSScalar>();
                  
                  Matrix<MSScalar, 5, 5> Qinv = Q.inverse().cast<MSScalar>();
                                          
                  // initialize active scalars for state parameters
                  Matrix<MSScalar, 5, 1> dum = Matrix<MSScalar, 5, 1>::Zero();
                  //suppress gradients of reference point parameters when fitting with gen constraint
                  for (unsigned int j=0; j<dum.size(); ++j) {
                    init_twice_active_var(dum[j], nlocal, localstateidx + j);
                  }

                  Matrix<MSScalar, 5, 1> du = Matrix<MSScalar, 5, 1>::Zero();
                  for (unsigned int j=0; j<du.size(); ++j) {
                    init_twice_active_var(du[j], nlocal, localstateidx + 5 + j);
                  }
                  
//                   MSScalar dbeta(corparms_[bfieldglobalidx]);
                  MSScalar dbeta(0.);
                  init_twice_active_var(dbeta, nlocal, localparmidx);
                  
//                   MSScalar dxi(corparms_[elossglobalidx]);
                  MSScalar dxi(0.);
                  init_twice_active_var(dxi, nlocal, localparmidx + 1);
                              
//                   const Matrix<MSScalar, 5, 1> dprop = dx0.cast<MSScalar>() + Hpstate*du - Hmstate*Fstate*dum - Hmstate*Fb*dbeta - Hmstate*Fxi*dxi;
                  Matrix<MSScalar, 5, 1> dprop;
                  if (dolocalupdate) {
                    dprop = dx0.cast<MSScalar>() + Hpstate*du - Hmstate*Fstate*dum - Hmstate*Fb*dbeta - Hmstate*Fxi*dxi;
                  }
                  else {
                    dprop = du - Fstate*dum - Fb*dbeta - Fxi*dxi;
                  }
                  const MSScalar chisq = dprop.transpose()*Qinv*dprop;
                  
                  chisq0val += chisq.value().value();
                    
                  auto const& gradlocal = chisq.value().derivatives();
                  //fill local hessian
                  Matrix<double, nlocal, nlocal> hesslocal;
                  for (unsigned int j=0; j<nlocal; ++j) {
                    hesslocal.row(j) = chisq.derivatives()[j].derivatives();
                  }
                  
                  constexpr std::array<unsigned int, 2> localsizes = {{ nlocalstate, nlocalparms }};
                  constexpr std::array<unsigned int, 2> localidxs = {{ localstateidx, localparmidx }};
                  const std::array<unsigned int, 2> fullidxs = {{ fullstateidx, fullparmidx }};
                  
                  for (unsigned int iidx = 0; iidx < localidxs.size(); ++iidx) {
                    gradfull.segment(fullidxs[iidx], localsizes[iidx]) += gradlocal.segment(localidxs[iidx], localsizes[iidx]);
                    for (unsigned int jidx = 0; jidx < localidxs.size(); ++jidx) {
                      hessfull.block(fullidxs[iidx], fullidxs[jidx], localsizes[iidx], localsizes[jidx]) += hesslocal.block(localidxs[iidx], localidxs[jidx], localsizes[iidx], localsizes[jidx]);
                    }
                  }
                  
                }
                
//                 const unsigned int bfieldglobalidx = detidparms.at(std::make_pair(6, parmdetid));
                globalidxv[parmidx] = bfieldglobalidx;
                parmidx++;
                
//                 const unsigned int elossglobalidx = detidparms.at(std::make_pair(7, parmdetid));
                globalidxv[parmidx] = elossglobalidx;
                parmidx++;

              }
              
              if (preciseHit->isValid()) {

                auto fillAlignGrads = [&](auto Nalign) {
                  constexpr unsigned int nlocalstate = 2;
                  constexpr unsigned int localstateidx = 0;
                  constexpr unsigned int localalignmentidx = nlocalstate;
                  constexpr unsigned int localparmidx = localalignmentidx;

                  // abusing implicit template argument to pass
                  // a template value via std::integral_constant
                  constexpr unsigned int nlocalalignment = Nalign();
                  constexpr unsigned int nlocalparms = nlocalalignment;
                  constexpr unsigned int nlocal = nlocalstate + nlocalparms;
                  
    //               std::cout << "idx = " << id << ", ihit = " << ihit << ", alignmentparmidx = " << alignmentparmidx << ", nlocalalignment = " << nlocalalignment << std::endl;

                  using AlignScalar = AANT<double, nlocal>;
                  
                  const unsigned int fullstateidx = trackstateidx + 5*ihit + 3;
                  const unsigned int fullparmidx = nstateparms + nparsBfield + nparsEloss + alignmentparmidx;

                  const bool ispixel = GeomDetEnumerators::isTrackerPixel(preciseHit->det()->subDetector());

                  //TODO add hit validation stuff
                  //TODO add simhit stuff

//                   const bool hit1d = preciseHit->dimension() == 1 || ispixel;
                  const bool hit1d = preciseHit->dimension() == 1;
                  
                  Matrix<AlignScalar, 2, 2> Hu = Hp.bottomRightCorner<2,2>().cast<AlignScalar>();

                  Matrix<AlignScalar, 2, 1> dy0;
                  Matrix<AlignScalar, 2, 2> Vinv;
                  // rotation from module to strip coordinates
      //             Matrix<AlignScalar, 2, 2> R;
                  Matrix2d R;

                  const double lxcor = localparms[3];
                  const double lycor = localparms[4];

      //             const double lxcor = localparms[3] - localconv[0];
      //             const double lycor = localparms[4] - localconv[1];

                  const Topology &topology = preciseHit->det()->topology();

                  // undo deformation correction
                  const LocalPoint lpnull(0., 0.);
                  const MeasurementPoint mpnull = topology.measurementPosition(lpnull);
                  const Topology::LocalTrackPred pred(tsostmp.localParameters().vector());

                  auto const defcorr = topology.localPosition(mpnull, pred) - topology.localPosition(mpnull);

                  const double hitx = preciseHit->localPosition().x() - defcorr.x();
                  const double hity = preciseHit->localPosition().y() - defcorr.y();


//                   if (preciseHit->dimension() == 1) {
                  if (hit1d) {
//                     dy0[0] = AlignScalar(preciseHit->localPosition().x() - localparms[3]);
                    dy0[0] = AlignScalar(hitx - lxcor);
                    dy0[1] = AlignScalar(0.);
                    
      //               bool simvalid = false;
      //               for (auto const& simhith : simHits) {
      //                 for (const PSimHit& simHit : *simhith) {
      //                   if (simHit.detUnitId() == preciseHit->geographicalId()) {                      
      //                     
      //                     dy0[0] = AlignScalar(simHit.localPosition().x() - updtsos.localPosition().x());
      //                     
      //                     simvalid = true;
      //                     break;
      //                   }
      //                 }
      //                 if (simvalid) {
      //                   break;
      //                 }
      //               }
                    
                    Vinv = Matrix<AlignScalar, 2, 2>::Zero();
                    Vinv(0,0) = 1./preciseHit->localPositionError().xx();
                    
      //               R = Matrix<AlignScalar, 2, 2>::Identity();
                    R = Matrix2d::Identity();
                  }
                  else {
                    // 2d hit
                    Matrix2d iV;
                    iV << preciseHit->localPositionError().xx(), preciseHit->localPositionError().xy(),
                          preciseHit->localPositionError().xy(), preciseHit->localPositionError().yy();
                    if (ispixel) {
//                     if (true) {
                      //take 2d hit as-is for pixels
//                       dy0[0] = AlignScalar(preciseHit->localPosition().x() - localparms[3]);
//                       dy0[1] = AlignScalar(preciseHit->localPosition().y() - localparms[4]);

                      dy0[0] = AlignScalar(hitx - lxcor);
                      dy0[1] = AlignScalar(hity - lycor);
                    
                      Vinv = iV.inverse().cast<AlignScalar>();
                      //FIXME various temporary hacks;
                      
      //                 dy0[1] = AlignScalar(0.);
      //                 Vinv = Matrix<AlignScalar, 2, 2>::Zero();
      //                 Vinv(0,0) = 1./preciseHit->localPositionError().xx();
                      
      //                 if (GeomDetEnumerators::isEndcap(preciseHit->det()->subDetector())) {
      //                 if (GeomDetEnumerators::isBarrel(preciseHit->det()->subDetector())) {
      //                   PXBDetId detidtest(preciseHit->det()->geographicalId());
      //                   int layertest = detidtest.layer();
      //                   
      //                   if (layertest > 1) {
      //                     Vinv = Matrix<AlignScalar, 2, 2>::Zero();
      //                   }
      //                   
      // //                   Vinv = Matrix<AlignScalar, 2, 2>::Zero();
      // //                   dy0[0] = AlignScalar(0.);
      // //                   dy0[1] = AlignScalar(0.);
      //                 }
                      
      //                 bool simvalid = false;
      //                 for (auto const& simhith : simHits) {
      //                   for (const PSimHit& simHit : *simhith) {
      //                     if (simHit.detUnitId() == preciseHit->geographicalId()) {                      
      //                       
      //                       if (GeomDetEnumerators::isBarrel(preciseHit->det()->subDetector())) {
      //                         dy0[0] = AlignScalar(simHit.localPosition().x() - updtsos.localPosition().x());
      //                         dy0[1] = AlignScalar(simHit.localPosition().y() - updtsos.localPosition().y());
      //                       }
      //                       
      // //                       dy0[1] = AlignScalar(0.);
      //                       
      //                       simvalid = true;
      //                       break;
      //                     }
      //                   }
      //                   if (simvalid) {
      //                     break;
      //                   }
      //                 }
                      
                      
      //                 R = Matrix<AlignScalar, 2, 2>::Identity();
                      R = Matrix2d::Identity();
                    }
                    else {

                      constexpr bool dopolar = true;
//                       constexpr bool dopolar = false;

                      if (dopolar) {
                        // transform to polar coordinates to end the madness
                        //TODO handle the module deformations consistently here (currently equivalent to dropping/undoing deformation correction)

                        const ProxyStripTopology *proxytopology = dynamic_cast<const ProxyStripTopology*>(&(preciseHit->det()->topology()));

//                         // undo deformation correction
//                         const Topology::LocalTrackPred pred(tsostmp.localParameters().vector());
//
//                         auto const defcorr = proxytopology->localPosition(0., pred) - proxytopology->localPosition(0.);
//
//                         const double hitx = preciseHit->localPosition().x() - defcorr.x();
//                         const double hity = preciseHit->localPosition().y() - defcorr.y();

                        const TkRadialStripTopology *radialtopology = dynamic_cast<const TkRadialStripTopology*>(&proxytopology->specificTopology());

                        const double rdir = radialtopology->yAxisOrientation();
                        const double radius = radialtopology->originToIntersection();

                        const double phihit = rdir*std::atan2(hitx, rdir*hity + radius);

                        // invert original calculation of covariance matrix to extract variance on polar angle
                        const double detHeight = radialtopology->detHeight();
                        const double radsigma = detHeight*detHeight/12.;

                        const double t1 = std::tan(phihit);
                        const double t2 = t1*t1;

                        const double tt = preciseHit->localPositionError().xx() - t2*radsigma;

                        const double phierr2 = tt / std::pow(radialtopology->centreToIntersection(), 2);

                        // TODO apply (inverse) corrections for module deformations here? (take into account for jacobian?)
                        const double phistate = rdir*std::atan2(lxcor, rdir*lycor + radius);

                        Vinv = Matrix<AlignScalar, 2, 2>::Zero();
                        Vinv(0, 0) = 1./phierr2;

                        R = Matrix2d::Zero();

                        const double yp = rdir*lycor + radius;
                        const double invden = 1./(lxcor*lxcor + yp*yp);

                        // dphi / dx
                        R(0, 0) = rdir*yp*invden;
                        // dphi / dy
                        R(0, 1) = -lxcor*invden;


                        dy0[0] = phihit - phistate;
                        dy0[1] = 0.;

                      }
                      else {
                        // standard treatment

                        // diagonalize and take only smallest eigenvalue for 2d hits in strip wedge modules,
                        // since the constraint parallel to the strip is spurious
                        SelfAdjointEigenSolver<Matrix2d> eigensolver(iV);

                        R = eigensolver.eigenvectors().transpose();
                        if (R(0,0) < 0.) {
                          R.row(0) *= -1.;
                        }
                        if (R(1,1) <0.) {
                          R.row(1) *= -1.;
                        }

                        Matrix<double, 2, 1> dy0local;

                        dy0local[0] = preciseHit->localPosition().x() - lxcor;
                        dy0local[1] = preciseHit->localPosition().y() - lycor;


                        const Matrix<double, 2, 1> dy0eig = R*dy0local;


                        //TODO deal properly with rotations (rotate back to module local coords?)
                        dy0[0] = AlignScalar(dy0eig[0]);
                        dy0[1] = AlignScalar(0.);

                        Vinv = Matrix<AlignScalar, 2, 2>::Zero();
                        Vinv(0,0) = AlignScalar(1./eigensolver.eigenvalues()[0]);
                      }

                    }
                  }
                  
    //               rxfull.row(ivalidhit) = R.row(0).cast<float>();
    //               ryfull.row(ivalidhit) = R.row(1).cast<float>();
                  
    //               validdxeigjac.block<2,2>(2*ivalidhit, 3*(ihit+1)) = R*Hp.bottomRightCorner<2,2>();
                  
                  const Matrix<AlignScalar, 2, 2> Ralign = R.cast<AlignScalar>();
                  
                  Matrix<AlignScalar, 2, 1> dx = Matrix<AlignScalar, 2, 1>::Zero();
                  for (unsigned int j=0; j<dx.size(); ++j) {
                    init_twice_active_var(dx[j], nlocal, localstateidx + j);
                  }

                  Matrix<AlignScalar, 6, 1> dalpha = Matrix<AlignScalar, 6, 1>::Zero();
                  // order in which to use parameters, especially relevant in case nlocalalignment < 6
//                   constexpr std::array<unsigned int, 6> alphaidxs = {{5, 0, 1, 2, 3, 4}};
                  constexpr std::array<unsigned int, 6> alphaidxs = {{0, 2, 3, 4, 5, 1}};


//                   for (unsigned int idim=0; idim<nlocalalignment; ++idim) {
//                     const unsigned int xglobalidx = detidparms.at(std::make_pair(alphaidxs[idim], preciseHit->geographicalId()));
//                     dalpha[alphaidxs[idim]] = AlignScalar(corparms_[xglobalidx]);
//                   }
                  
                  for (unsigned int idim=0; idim<nlocalalignment; ++idim) {
      //               init_twice_active_var(dalpha[idim], nlocal, localalignmentidx+idim);
                    init_twice_active_var(dalpha[alphaidxs[idim]], nlocal, localalignmentidx+idim);
                  }
                  
                  // alignment jacobian
                  Matrix<AlignScalar, 2, 6> A = Matrix<AlignScalar, 2, 6>::Zero();

                  const double localqopval = localparms[0];
                  const double localdxdzval = localparms[1];
                  const double localdydzval = localparms[2];
                  const double localxval = localparms[3];
                  const double localyval = localparms[4];
                              
                  // dx/dx
                  A(0,0) = AlignScalar(1.);
                  // dy/dy
                  A(1,1) = AlignScalar(1.);
                  // dx/dz
                  A(0,2) = localdxdzval;
                  // dy/dz
                  A(1,2) = localdydzval;
                  // dx/dtheta_x
                  A(0,3) = -localyval*localdxdzval;
                  // dy/dtheta_x
                  A(1,3) = -localyval*localdydzval;
                  // dx/dtheta_y
                  A(0,4) = -localxval*localdxdzval;
                  // dy/dtheta_y
                  A(1,4) = -localxval*localdydzval;
                  // dx/dtheta_z
                  A(0,5) = -localyval;
                  // dy/dtheta_z
                  A(1,5) = localxval;
                  
                  
      //             std::cout << "strip local z shift gradient: " << (Ralign*A.col(2))[0].value().value() << std::endl;
                  
                  // rotation from alignment basis to module local coordinates
      //             Matrix<AlignScalar, 2, 2> A;
      //             if (isglued) {
      //               const GlobalVector modx = preciseHit->det()->surface().toGlobal(LocalVector(1.,0.,0.));
      //               const GlobalVector mody = preciseHit->det()->surface().toGlobal(LocalVector(0.,1.,0.));
      //               
      //               const GlobalVector gluedx = parmDet->surface().toGlobal(LocalVector(1.,0.,0.));
      //               const GlobalVector gluedy = parmDet->surface().toGlobal(LocalVector(0.,1.,0.));
      //               
      //               A(0,0) = AlignScalar(modx.dot(gluedx));
      //               A(0,1) = AlignScalar(modx.dot(gluedy));
      //               A(1,0) = AlignScalar(mody.dot(gluedx));
      //               A(1,1) = AlignScalar(mody.dot(gluedy));
      //             }
      //             else {
      //               A = Matrix<AlignScalar, 2, 2>::Identity();
      //             }
      // 
      //             Matrix<AlignScalar, 2, 1> dh = dy0 - R*Hu*dx - R*A*dalpha;
                  
                            
                  double thetaincidence = std::asin(1./std::sqrt(std::pow(localdxdzval,2) + std::pow(localdydzval,2) + 1.));

                  
      //             bool morehitquality = applyHitQuality_ ? thetaincidence > 0.25 : true;
                  bool morehitquality = true;
                  
                  if (morehitquality) {
                    nValidHitsFinal++;
                    if (ispixel) {
                      nValidPixelHitsFinal++;
                    }
                  }
                  else {
                    Vinv = Matrix<AlignScalar, 2, 2>::Zero();
                  }

                  Matrix<AlignScalar, 2, 1> dh = dy0 - Ralign*Hu*dx - Ralign*A*dalpha;
                  AlignScalar chisq = dh.transpose()*Vinv*dh;
                  
//                   std::cout << "icons = " << icons << " iiter = " << iiter << " ihit = " << ihit << " chisq = " << chisq.value().value() << std::endl;
                  
                  chisq0val += chisq.value().value();

                  auto const& gradlocal = chisq.value().derivatives();
                  //fill local hessian
                  Matrix<double, nlocal, nlocal> hesslocal;
                  for (unsigned int j=0; j<nlocal; ++j) {
                    hesslocal.row(j) = chisq.derivatives()[j].derivatives();
                  }
                  
                  constexpr std::array<unsigned int, 2> localsizes = {{ nlocalstate, nlocalparms }};
                  constexpr std::array<unsigned int, 2> localidxs = {{ localstateidx, localparmidx }};
                  const std::array<unsigned int, 2> fullidxs = {{ fullstateidx, fullparmidx }};
                  
                  for (unsigned int iidx = 0; iidx < localidxs.size(); ++iidx) {
                    gradfull.segment(fullidxs[iidx], localsizes[iidx]) += gradlocal.segment(localidxs[iidx], localsizes[iidx]);
                    for (unsigned int jidx = 0; jidx < localidxs.size(); ++jidx) {
                      hessfull.block(fullidxs[iidx], fullidxs[jidx], localsizes[iidx], localsizes[jidx]) += hesslocal.block(localidxs[iidx], localidxs[jidx], localsizes[iidx], localsizes[jidx]);
                    }
                  }
                  
      //             Matrix<double, nlocal, 1> gradloctest0;
      //             Matrix<double, 1, 1> gradloctest1;
      //             Matrix<double, 2, 1> gradloctest2;
                  
      //             std::cout << "nlocalalignment: " << nlocalalignment << " nlocal: " << nlocal << std::endl;
      //             std::cout << "gradlocal type: " << typeid(gradlocal).name() << std::endl;
      //             std::cout << "gradloctest0 type: " << typeid(gradloctest0).name() << std::endl;
      //             std::cout << "gradloctest1 type: " << typeid(gradloctest1).name() << std::endl;
      //             std::cout << "gradloctest2 type: " << typeid(gradloctest2).name() << std::endl;
      //             
      //             std::cout << "nhits: " << nhits << " nvalid: " << nvalid << " nvalidalign2d: " << nvalidalign2d << " ihit: " << ihit << std::endl;
      //             std::cout << "gradfull.size(): " << gradfull.size() << " nlocalstate: " << nlocalstate << " fullstateidx: " << fullstateidx << " nlocalparms: " << nlocalparms << " fullparmidx: " << fullparmidx << std::endl;
                  
                  // FIXME the templated block functions don't work here for some reason
                  //fill global gradient
    //               gradfull.segment<nlocalstate>(fullstateidx) += gradlocal.head(nlocalstate);
    //               gradfull.segment<nlocalparms>(fullparmidx) += gradlocal.segment(localparmidx, nlocalparms);
    // 
    //               //fill global hessian (upper triangular blocks only)
    //               hessfull.block<nlocalstate,nlocalstate>(fullstateidx, fullstateidx) += hesslocal.topLeftCorner(nlocalstate,nlocalstate);
    //               hessfull.block<nlocalstate,nlocalparms>(fullstateidx, fullparmidx) += hesslocal.topRightCorner(nlocalstate, nlocalparms);
    //               hessfull.block<nlocalparms, nlocalparms>(fullparmidx, fullparmidx) += hesslocal.bottomRightCorner(nlocalparms, nlocalparms);
                  
                  for (unsigned int idim=0; idim<nlocalalignment; ++idim) {
                    const unsigned int xglobalidx = detidparms.at(std::make_pair(alphaidxs[idim], preciseHit->geographicalId()));
                    globalidxv[nparsBfield + nparsEloss + alignmentparmidx] = xglobalidx;
                    alignmentparmidx++;
                    if (alphaidxs[idim]==0) {
                      hitidxv.push_back(xglobalidx);
                    }
                  }
                };
                
                
//                 if (align2d) {
//                   fillAlignGrads(std::integral_constant<unsigned int, 3>());
//                 }
//                 else {
//                   fillAlignGrads(std::integral_constant<unsigned int, 2>());
//                 }
                if (align2d) {
                  fillAlignGrads(std::integral_constant<unsigned int, 6>());
                }
                else {
                  fillAlignGrads(std::integral_constant<unsigned int, 5>());
                }
//                 fillAlignGrads(std::integral_constant<unsigned int, 6>());
                
              }
              
            }

            if (!valid) {
              break;
            }
            
            trackstateidx += 5*tracknhits;
          }
          
          if (!valid) {
            break;
          }
          
    //       MatrixXd massjac;

          const Matrix<double, 7, 1> &refFts0 = refftsarr[0];
          const Matrix<double, 7, 1> &refFts1 = refftsarr[1];    
          
          const Matrix<double, 6, 6> mhess = massHessianAltD(refFts0, refFts1, mmu);
          const Matrix<double, 6, 6> mhessinvsq = massinvsqHessianAltD(refFts0, refFts1, mmu);
          
          const double dmassconv = iiter > 0 ? 0.5*(mhess*covrefmom).trace() : 0.;
          const double dmassconvinvsq = iiter > 0 ? 0.5*(mhessinvsq*covrefmom).trace() : 0.;
          
          
          dmassconvval = dmassconv;
          dinvmasssqconvval = dmassconvinvsq;
          
          if (icons == 0) {
            dmassconvval_cons0 = dmassconv;
            dinvmasssqconvval_cons0 = dmassconvinvsq;
          }
          
          // add mass constraint to gbl fit
//           if (icons == 1 && iiter > 3) {
          if (icons > 0) {
//           if (false) {
            //TODO simplify this to treat the 6 parameters contiguously (now that they are contiguous in the original vector)

            constexpr unsigned int nlocalstate0 = 3;
            constexpr unsigned int nlocalstate1 = 3;
            
            constexpr unsigned int nlocal = nlocalstate0 + nlocalstate1;
            
            constexpr unsigned int localstateidx0 = 0;
            constexpr unsigned int localstateidx1 = localstateidx0 + nlocalstate0;
            
            const unsigned int fullstateidx0 = 0;
            const unsigned int fullstateidx1 = 3;
            
            using MScalar = AANT<double, nlocal>;
            
            //TODO optimize to avoid recomputation of FTS
    //         const FreeTrajectoryState refFts0 = outparts[0]->currentState().freeTrajectoryState();
    //         const FreeTrajectoryState refFts1 = outparts[1]->currentState().freeTrajectoryState();
            
//             const FreeTrajectoryState &refFts0 = refftsarr[0];
//             const FreeTrajectoryState &refFts1 = refftsarr[1];
            
            const ROOT::Math::PxPyPzMVector mom0(refFts0[3],
                                                    refFts0[4],
                                                    refFts0[5],
                                                    mmu);
            
            const ROOT::Math::PxPyPzMVector mom1(refFts1[3],
                                                    refFts1[4],
                                                    refFts1[5],
                                                    mmu);
            
            const double massval = (mom0 + mom1).mass();
            
            const Matrix<double, 1, 6> mjacalt = massJacobianAltD(refFts0, refFts1, mmu);
//             const Matrix<double, 1, 6> mjacaltinvsq = massinvsqJacobianAltD(refFts0, refFts1, mmu);
            
  //           const double mrval = 1./massval/massval;
            
  //           const double mrconstraintval = 1./massconstraintval/massconstraintval;
            
  //           const double 
            
            const Matrix<double, 5, 7> &FdFmref0 = FdFmrefarr[0];
            const Matrix<double, 5, 7> &FdFmref1 = FdFmrefarr[1];
            
//             JacobianCurvilinearToCartesian curv2cartref0(refFts0.parameters());
//             auto const &jacCurv2Cartref0 = Map<const Matrix<double, 6, 5, RowMajor>>(curv2cartref0.jacobian().Array());
//             
//             JacobianCurvilinearToCartesian curv2cartref1(refFts1.parameters());
//             auto const &jacCurv2Cartref1 = Map<const Matrix<double, 6, 5, RowMajor>>(curv2cartref1.jacobian().Array());
            
  //           const Matrix<double, 1, 6> m2jac = massJacobian(refFts0, refFts1, mmu);
            
  //           std::cout << "massval = " << massval << std::endl;
            

            
//             double dmassconv = 0.;
//             
//             if (iiter > 0) {
// 
//               
//               //mhess = Q Lambda Q^-1
//               // < dx^T Q sqrt(Lambda) sqrt(Lambda) Q^-1 dx >
//               // y = sqrt(Lambda) Q^-1 dx
//               // sigmaY = sqrt(Lambda) Q^-1 sigmax Q sqrt(Lambda)
//               // <> = Tr(sqrt(Lambda) Q^-1 sigmax Q sqrt(Lambda))
//               // <> = Tr(Lambda Q^-1 sigmax Q)
//               
//               
//               SelfAdjointEigenSolver<Matrix<double, 6, 6>> eshess(mhess);
//               
// //               std::cout << "eshess:\n" << eshess.eigenvalues() << std::endl;
//               
//               dmassconv = 0.5*(eshess.eigenvalues().asDiagonal()*eshess.eigenvectors().transpose()*covrefmom*eshess.eigenvectors()).trace();
//               
//               const double dmassconvalt = 0.5*(mhess*covrefmom).trace();
// 
// //               std::cout << "iiter = " << iiter << " dmassconv = " << dmassconv << " relative = " << dmassconv/massconstraintval << std::endl;
//               std::cout << "iiter = " << iiter << " dmassconv = " << dmassconv << " relative = " << dmassconv/massconstraintval << " dmassconvalt = " << dmassconvalt << std::endl;
// 
// 
//             }
//             
//             dmassconvval = dmassconv;
            
//             const Matrix<double, 1, 6> mjacalt = (1./massval)*massJacobianAlt(refFts0, refFts1, mmu);
            
            
  //           const Matrix<double, 1, 6> mjacalt = mrJacobian(refFts0, refFts1, mmu);
            
  //           const Matrix<double, 1, 3> mjac0 = m2jac.leftCols<3>()*jacCurv2Cartref0.bottomLeftCorner<3, 3>();
  //           const Matrix<double, 1, 3> mjac1 = m2jac.rightCols<3>()*jacCurv2Cartref1.bottomLeftCorner<3, 3>();
            
  //           const Matrix<double, 1, 3> mjacalt0 = mjacalt.leftCols<3>();
  //           const Matrix<double, 1, 3> mjacalt1 = mjacalt.rightCols<3>();
            
  //           std::cout << "mjac0" << std::endl;
  //           std::cout << mjac0 << std::endl;
  //           std::cout << "mjacalt0" << std::endl;
  //           std::cout << mjacalt0 << std::endl;
  // 
  //           std::cout << "mjac1" << std::endl;
  //           std::cout << mjac1 << std::endl;
  //           std::cout << "mjacalt1" << std::endl;
  //           std::cout << mjacalt1 << std::endl;
                              
            // initialize active scalars for common vertex parameters
//             Matrix<MScalar, 3, 1> dvtx = Matrix<MScalar, 3, 1>::Zero();
//             for (unsigned int j=0; j<dvtx.size(); ++j) {
//               init_twice_active_var(dvtx[j], nlocal, localvtxidx + j);
//             }

            // initialize active scalars for state parameters
            // (first track is index together with vertex parameters)
            
            Matrix<MScalar, 3, 1> dmomcurv0 = Matrix<MScalar, 3, 1>::Zero();
            for (unsigned int j=0; j<dmomcurv0.size(); ++j) {
              init_twice_active_var(dmomcurv0[j], nlocal, localstateidx0 + j);
            }
            
            Matrix<MScalar, 3, 1> dmomcurv1 = Matrix<MScalar, 3, 1>::Zero();
            for (unsigned int j=0; j<dmomcurv1.size(); ++j) {
              init_twice_active_var(dmomcurv1[j], nlocal, localstateidx1 + j);
            }
            
  //           const Matrix<MScalar, 3, 1> dmom0 = jacCurv2Cartref0.bottomLeftCorner<3, 3>().cast<MScalar>()*dmomcurv0;
  //           const Matrix<MScalar, 3, 1> dmom1 = jacCurv2Cartref1.bottomLeftCorner<3, 3>().cast<MScalar>()*dmomcurv1;
            
            // resonance width
    //         const MScalar invSigmaMsq(0.25/massConstraint_/massConstraint_/massConstraintWidth_/massConstraintWidth_);
    //         const MScalar dmsq0 = MScalar(m0*m0 - massConstraint_*massConstraint_);
            
//             const MScalar invSigmaMsq(1./massConstraintWidth_/massConstraintWidth_);
//             const MScalar dmsq0 = MScalar(massval - massconstraintval - dmassconv);

            const MScalar invSigmaMsq(1./massConstraintWidth_/massConstraintWidth_);
//             const MScalar dmsq0 = MScalar(1./massval/massval - massconstraintval - dmassconvinvsq);
//             const MScalar dmsq0 = MScalar(1./massval/massval - massconstraintval);
            
            
            const MScalar dmsq0 = MScalar(massval - massconstraintval - dmassconv);
//             const MScalar dmsq0 = MScalar(massval - massconstraintval);
            
//             const double sigmalnm = massConstraintWidth_/massconstraintval;
//             const MScalar invSigmaMsq(1./sigmalnm/sigmalnm);
//             const MScalar dmsq0 = MScalar(std::log(massval) - std::log(massconstraintval));
            
  //           const MScalar invSigmaMsq(0.25*std::pow(massval, 6)/massConstraintWidth_/massConstraintWidth_);
  //           const MScalar dmsq0 = MScalar(mrval - mrconstraintval);
  //           const MScalar dmsq0 = MScalar(massval - massConstraint_);
  //           const MScalar dmsq0 = MScalar(0.);
    //         const MScalar dmsq0 = MScalar(m0 - massConstraint_);
            
            
    //         std::cout << "invSigmaMsq = " << invSigmaMsq.value().value() << std::endl;
            
  //           const MScalar dmsqtrack0 = (m2jac.leftCols<3>().cast<MScalar>()*dmom0)[0];
  //           const MScalar dmsqtrack1 = (m2jac.rightCols<3>().cast<MScalar>()*dmom1)[0];
            
            const MScalar dmsqtrack0 = (mjacalt.leftCols<3>().cast<MScalar>()*dmomcurv0)[0];
            const MScalar dmsqtrack1 = (mjacalt.rightCols<3>().cast<MScalar>()*dmomcurv1)[0];

//             const MScalar dmsqtrack0 = (mjacaltinvsq.leftCols<3>().cast<MScalar>()*dmomcurv0)[0];
//             const MScalar dmsqtrack1 = (mjacaltinvsq.rightCols<3>().cast<MScalar>()*dmomcurv1)[0];
            
            const MScalar dmsq = dmsq0 + dmsqtrack0 + dmsqtrack1;
            
    //         const MScalar chisq = dmsq*invSigmaMsq*dmsq;
            const MScalar chisq = invSigmaMsq*dmsq*dmsq;
            
            chisq0val += chisq.value().value();
            
            auto const& gradlocal = chisq.value().derivatives();
            //fill local hessian
            Matrix<double, nlocal, nlocal> hesslocal;
            for (unsigned int j=0; j<nlocal; ++j) {
              hesslocal.row(j) = chisq.derivatives()[j].derivatives();
            }
            
            constexpr std::array<unsigned int, 2> localsizes = {{ nlocalstate0, nlocalstate0 }};
            constexpr std::array<unsigned int, 2> localidxs = {{ localstateidx0, localstateidx1 }};
            const std::array<unsigned int, 2> fullidxs = {{ fullstateidx0, fullstateidx1 }};
            
            for (unsigned int iidx = 0; iidx < localidxs.size(); ++iidx) {
              gradfull.segment(fullidxs[iidx], localsizes[iidx]) += gradlocal.segment(localidxs[iidx], localsizes[iidx]);
              for (unsigned int jidx = 0; jidx < localidxs.size(); ++jidx) {
                hessfull.block(fullidxs[iidx], fullidxs[jidx], localsizes[iidx], localsizes[jidx]) += hesslocal.block(localidxs[iidx], localidxs[jidx], localsizes[iidx], localsizes[jidx]);
              }
            }
            
          }

         

          
          
    //       std::cout << nhits << std::endl;
    //       std::cout << nvalid << std::endl;
    //       std::cout << nvalidalign2d << std::endl;
    //       std::cout << nparsAlignment << std::endl;
    //       std::cout << alignmentparmidx << std::endl;
    // 
    //       std::cout << nparsBfield << std::endl;
    //       std::cout << nparsEloss << std::endl;
    //       std::cout << parmidx << std::endl;
          
          assert(trackstateidx == nstateparms);
          assert(parmidx == (nparsBfield + nparsEloss));
          assert(alignmentparmidx == nparsAlignment);
          
    //       if (nhits != nvalid) {
    //         continue;
    //       }

          auto freezeparm = [&](unsigned int idx) {
            gradfull[idx] = 0.;
            hessfull.row(idx) *= 0.;
            hessfull.col(idx) *= 0.;
            hessfull(idx,idx) = 1e6;
          };
          
          if (fitFromGenParms_) {
            for (unsigned int i=0; i<10; ++i) {
              freezeparm(i);
            }
          }

          if (doVtxConstraint_) {
            // freeze track pca distance to 0 to impose common vertex constraint
            freezeparm(6);
          }
          
//           if (fitFromGenParms_) {
//             freezeparm(2);
//             for (unsigned int id = 0; id < 2; ++id) {
//               freezeparm(trackstateidxarr[id] + 1);
//             }
//           }
          
  //         if (fitFromGenParms_) {
  //           for (unsigned int id = 0; id < 2; ++id) {
  //             for (unsigned int i=1; i<3; ++i) {
  //               freezeparm(trackstateidxarr[id] + i);
  //             }
  //           }
  //         }
          
          //now do the expensive calculations and fill outputs
          
          //symmetrize the matrix (previous block operations do not guarantee that the needed blocks are filled)
          //TODO handle this more efficiently?
    //       hessfull.triangularView<StrictlyLower>() = hessfull.triangularView<StrictlyUpper>().transpose();
          
    //       for (unsigned int i=0; i<3; ++i) {
    //         gradfull[i] = 0.;
    //         hessfull.row(i) *= 0.;
    //         hessfull.col(i) *= 0.;
    //         hessfull(i,i) = 1e6;
    //       }
          
    //       for (auto trackstateidx : trackstateidxarr) {
    //         for (unsigned int i = trackstateidx; i < (trackstateidx + 1); ++i) {
    //           gradfull[i] = 0.;
    //           hessfull.row(i) *= 0.;
    //           hessfull.col(i) *= 0.;
    //           hessfull(i,i) = 1e6;
    //         }
    //       }
          
    //       {
    //         unsigned int i = trackstateidxarr[1];
    //         gradfull[i] = 0.;
    //         hessfull.row(i) *= 0.;
    //         hessfull.col(i) *= 0.;
    //         hessfull(i,i) = 1e6; 
    //       }
    //       
          
    //       std::cout << "gradfull:" << std::endl;
    //       std::cout << gradfull << std::endl;
    //       
    //       std::cout << "gradfull.head(nstateparms):" << std::endl;
    //       std::cout << gradfull.head(nstateparms) << std::endl;
    // 
    //       std::cout << "gradfull.tail(npars):" << std::endl;
    //       std::cout << gradfull.tail(npars) << std::endl;
    //       
    //       std::cout << "hessfull.diagonal():" << std::endl;
    //       std::cout << hessfull.diagonal() << std::endl;
          
          auto const& dchisqdx = gradfull.head(nstateparms);
          auto const& dchisqdparms = gradfull.tail(npars);
          
          auto const& d2chisqdx2 = hessfull.topLeftCorner(nstateparms, nstateparms);
          auto const& d2chisqdxdparms = hessfull.topRightCorner(nstateparms, npars);
          auto const& d2chisqdparms2 = hessfull.bottomRightCorner(npars, npars);
          

          
          Cinvd.compute(d2chisqdx2);
          
          dxfull = -Cinvd.solve(dchisqdx);
//           dxdparms = -Cinvd.solve(d2chisqdxdparms).transpose();
          
          
//           std::cout << "dxfull vtx: " << dxfull.head<3>() << std::endl;
          
  //         dxdparms = -Cinvd.solve(d2chisqdxdparms).transpose();
          
      //     if (debugprintout_) {
      //       std::cout << "dxrefdparms" << std::endl;
      //       std::cout << dxdparms.leftCols<5>() << std::endl;
      //     }
          
  //         grad = dchisqdparms + dxdparms*dchisqdx;
          //TODO check the simplification
      //     hess = d2chisqdparms2 + 2.*dxdparms*d2chisqdxdparms + dxdparms*d2chisqdx2*dxdparms.transpose();
  //         hess = d2chisqdparms2 + dxdparms*d2chisqdxdparms;
          
          const Matrix<double, 1, 1> deltachisq = dchisqdx.transpose()*dxfull + 0.5*dxfull.transpose()*d2chisqdx2*dxfull;
          
//           std::cout << "iiter = " << iiter << ", deltachisq = " << deltachisq[0] << std::endl;
//   //         
//           SelfAdjointEigenSolver<MatrixXd> es(d2chisqdx2, EigenvaluesOnly);
//           const double condition = es.eigenvalues()[nstateparms-1]/es.eigenvalues()[0];
//           std::cout << "eigenvalues:" << std::endl;
//           std::cout << es.eigenvalues().transpose() << std::endl;
//           std::cout << "condition: " << condition << std::endl;
          
          chisqval = chisq0val + deltachisq[0];
          
          deltachisqval = chisq0val + deltachisq[0] - chisqvalold;
          
          chisqvalold = chisq0val + deltachisq[0];
          
//           ndof = 5*nhits + nvalid + nvalidalign2d - nstateparms;
          ndof = 5*nhits + nvalid + nvalidpixel - nstateparms;
          
          if (bsConstraint_) {
            ndof += 3;
          }

          if (doVtxConstraint_) {
            ++ndof;
          }
          
          if (icons == 1) {
            ++ndof;
          }
          
//           std::cout << "icons = " << icons << " iiter =" << iiter << " dx = " << refftsarr[0].position() - refftsarr[1].position() << std::endl;
          
    //       std::cout << "dchisqdparms.head<6>()" << std::endl;
    //       std::cout << dchisqdparms.head<6>() << std::endl;
    //       
    //       std::cout << "grad.head<6>()" << std::endl;
    //       std::cout << grad.head<6>() << std::endl;
    //       
    //       std::cout << "d2chisqdparms2.topLeftCorner<6, 6>():" << std::endl;
    //       std::cout << d2chisqdparms2.topLeftCorner<6, 6>() << std::endl;
    //       std::cout << "hess.topLeftCorner<6, 6>():" << std::endl;
    //       std::cout << hess.topLeftCorner<6, 6>() << std::endl;
    //       
    //       std::cout << "dchisqdparms.segment<6>(nparsBfield+nparsEloss)" << std::endl;
    //       std::cout << dchisqdparms.segment<6>(nparsBfield+nparsEloss) << std::endl;
    //       
    //       std::cout << "grad.segment<6>(nparsBfield+nparsEloss)" << std::endl;
    //       std::cout << grad.segment<6>(nparsBfield+nparsEloss) << std::endl;
    //       
    //       std::cout << "d2chisqdparms2.block<6, 6>(nparsBfield+nparsEloss, nparsBfield+nparsEloss):" << std::endl;
    //       std::cout << d2chisqdparms2.block<6, 6>(nparsBfield+nparsEloss, nparsBfield+nparsEloss) << std::endl;
    //       std::cout << "hess.block<6, 6>(nparsBfield+nparsEloss, nparsBfield+nparsEloss):" << std::endl;
    //       std::cout << hess.block<6, 6>(nparsBfield+nparsEloss, nparsBfield+nparsEloss) << std::endl;
    // //       
    //       
    //       std::cout << "d2chisqdparms2.block<6, 6>(trackparmidxarr[1], trackparmidxarr[1]):" << std::endl;
    //       std::cout << d2chisqdparms2.block<6, 6>(trackparmidxarr[1], trackparmidxarr[1]) << std::endl;
    //       std::cout << "hess.block<6, 6>(trackparmidxarr[1], trackparmidxarr[1]):" << std::endl;
    //       std::cout << hess.block<6, 6>(trackparmidxarr[1], trackparmidxarr[1]) << std::endl;
    //       
    //       std::cout << "d2chisqdparms2.bottomRightCorner<6, 6>():" << std::endl;
    //       std::cout << d2chisqdparms2.bottomRightCorner<6, 6>() << std::endl;
    //       std::cout << "hess.bottomRightCorner<6, 6>():" << std::endl;
    //       std::cout << hess.bottomRightCorner<6, 6>() << std::endl;

    //       const double 
    // //       const double corxi0plusminus = hess(1, trackparmidxarr[1] + 1)/std::sqrt(hess(1,1)*hess(trackparmidxarr[1] + 1, trackparmidxarr[1] + 1));
    // //       const double corxi1plusminus = hess(3, trackparmidxarr[1] + 3)/std::sqrt(hess(3,3)*hess(trackparmidxarr[1] + 3, trackparmidxarr[1] + 3));
    //       
    //       const double cor01plus = hess(1, 3)/std::sqrt(hess(1, 1)*hess(3, 3));
    // //       const double cor01minus = hess(trackparmidxarr[1] + 1, trackparmidxarr[1] + 3)/std::sqrt(hess(trackparmidxarr[1] + 1, trackparmidxarr[1] + 1)*hess(trackparmidxarr[1] + 3, trackparmidxarr[1] + 3));
    // 
    //       const double cor12plus = hess(3, 5)/std::sqrt(hess(3, 3)*hess(5, 5));
    // //       const double cor12minus = hess(trackparmidxarr[1] + 3, trackparmidxarr[1] + 5)/std::sqrt(hess(trackparmidxarr[1] + 3, trackparmidxarr[1] + 3)*hess(trackparmidxarr[1] + 5, trackparmidxarr[1] + 5));
    //       
    // //       std::cout << "corxi0plusminus = " << corxi0plusminus << std::endl;
    // //       std::cout << "corxi1plusminus = " << corxi1plusminus << std::endl;
    //       std::cout << "cor01plus = " << cor01plus << std::endl;
    // //       std::cout << "cor01minus = " << cor01minus << std::endl;
    //       std::cout << "cor12plus = " << cor12plus << std::endl;
    // //       std::cout << "cor12minus = " << cor12minus << std::endl;
          
    //       std::cout << "hess(1, 1)" << std::endl;
    //       std::cout << hess(1, 1) << std::endl;
    //       std::cout << "hess(trackparmidxarr[1] + 1, trackparmidxarr[1] + 1)" << std::endl;
    //       std::cout << hess(trackparmidxarr[1] + 1, trackparmidxarr[1] + 1) << std::endl;
    //       std::cout << "hess(1, trackparmidxarr[1] + 1)" << std::endl;
    //       std::cout << hess(1, trackparmidxarr[1] + 1) << std::endl;
          
          // compute final kinematics
          
          kinTree->movePointerToTheTop();
          RefCountedKinematicVertex dimu_vertex = kinTree->currentDecayVertex();
          
          if (icons == 0) {
            Jpsikin_x = dimu_vertex->position().x();
            Jpsikin_y = dimu_vertex->position().y();
            Jpsikin_z = dimu_vertex->position().z();
          }
          else {
            Jpsikincons_x = dimu_vertex->position().x();
            Jpsikincons_y = dimu_vertex->position().y();
            Jpsikincons_z = dimu_vertex->position().z(); 
          }
          
          // apply the GBL fit results to the vertex position
          const Matrix<double, 10, 1> statepca = twoTrackCart2pca(refftsarr[0], refftsarr[1]);
          const Matrix<double, 10, 1> statepcaupd = statepca + dxfull.head<10>();

//           std::cout << "statepcaupd d = " << statepcaupd[6] << std::endl;

          const bool firstplus = statepcaupd[0] > 0.;

          if (icons == 0) {
            // define sign of d wrt charge of tracks
            Jpsi_d = firstplus ? statepcaupd[6] : -statepcaupd[6];
            Jpsi_x = statepcaupd[7];
            Jpsi_y = statepcaupd[8];
            Jpsi_z = statepcaupd[9];
          }
          else {
            // define sign of d wrt charge of tracks
            Jpsicons_d = firstplus ? statepcaupd[6] : -statepcaupd[6];
            Jpsicons_x = statepcaupd[7];
            Jpsicons_y = statepcaupd[8];
            Jpsicons_z = statepcaupd[9];
          }

          std::array<ROOT::Math::PxPyPzMVector, 2> muarr;
          std::array<Vector3d, 2> mucurvarr;
//           std::array<int, 2> muchargearr;
          
    //       std::cout << dimu_vertex->position() << std::endl;
          
          // apply the GBL fit results to the muon kinematics
          for (unsigned int id = 0; id < 2; ++id) {
            const Matrix<double, 7, 1> &refFts = refftsarr[id];
            const unsigned int trackstateidx = 3*id;

            const double qbpupd = statepcaupd[trackstateidx];
            const double lamupd = statepcaupd[trackstateidx + 1];
            const double phiupd = statepcaupd[trackstateidx + 2];
            
            const double charge = std::copysign(1., qbpupd);
            const double pupd = std::abs(1./qbpupd);
            
            const double pxupd = pupd*std::cos(lamupd)*std::cos(phiupd);
            const double pyupd = pupd*std::cos(lamupd)*std::sin(phiupd);
            const double pzupd = pupd*std::sin(lamupd);
            
            muarr[id] = ROOT::Math::PxPyPzMVector(pxupd, pyupd, pzupd, mmu);
            muchargearr[id] = charge;
            
            auto &refParms = mucurvarr[id];
            refParms << qbpupd, lamupd, phiupd;
            
    //         auto const &refFts = outparts[id]->currentState().freeTrajectoryState();
//             auto const &refFts = refftsarr[id];
//   //           auto const &jac = jacarr[id];
//             unsigned int trackstateidx = trackstateidxarr[id];
//             
// //             JacobianCurvilinearToCartesian curv2cart(refFts.parameters());
// //             const AlgebraicMatrix65& jac = curv2cart.jacobian();
// //             const Matrix<double, 6, 5> jac = curv2cartJacobianAlt(refFts);
//             const AlgebraicVector6 glob = refFts.parameters().vector();
//             
//             const Matrix<double, 3, 1> posupd = Map<const Matrix<double, 6, 1>>(glob.Array()).head<3>() + dxfull.head<3>();
//             
// //             const Matrix<double, 3, 1> momupd = Map<const Matrix<double, 6, 1>>(glob.Array()).tail<3>() + Map<const Matrix<double, 6, 5, RowMajor>>(jac.Array()).bottomLeftCorner<3, 3>()*dxfull.segment<3>(trackstateidx);
// //             const Matrix<double, 3, 1> momupd = Map<const Matrix<double, 6, 1>>(glob.Array()).tail<3>() + jac.bottomLeftCorner<3, 3>()*dxfull.segment<3>(trackstateidx);
//             
//             const GlobalPoint pos(posupd[0], posupd[1], posupd[2]);
// //             const GlobalVector mom(momupd[0], momupd[1], momupd[2]);
// //             const double charge = std::copysign(1., refFts.charge()/refFts.momentum().mag() + dxfull[trackstateidx]);
//       //         std::cout << "before update: reffts:" << std::endl;
//       //         std::cout << refFts.parameters().vector() << std::endl;
//       //         std::cout << "charge " << refFts.charge() << std::endl;
//   //           updFts = FreeTrajectoryState(pos, mom, charge, field);
// 
//             
//             const CurvilinearTrajectoryParameters curv(refFts.position(), refFts.momentum(), refFts.charge());
//               
//             const double qbpupd = curv.Qbp() + dxfull(trackstateidx);
//             const double lamupd = curv.lambda() + dxfull(trackstateidx + 1);
//             const double phiupd = curv.phi() + dxfull(trackstateidx + 2);
//             
//             const double charge = std::copysign(1., qbpupd);
//             const double pupd = std::abs(1./qbpupd);
//             
//             const double pxupd = pupd*std::cos(lamupd)*std::cos(phiupd);
//             const double pyupd = pupd*std::cos(lamupd)*std::sin(phiupd);
//             const double pzupd = pupd*std::sin(lamupd);
//             
//             const GlobalVector mom(pxupd, pyupd, pzupd);
//             
//             muarr[id] = ROOT::Math::PxPyPzMVector(pxupd, pyupd, pzupd, mmu);
//             muchargearr[id] = charge;
//             
// //             std::cout << "delta eta final = " << muarr[id].eta() - refFts.momentum().eta() << std::endl;
//                     
//             auto &refParms = mucurvarr[id];
//   //           CurvilinearTrajectoryParameters curvparms(refFts.position(), refFts.momentum(), refFts.charge());
//             CurvilinearTrajectoryParameters curvparms(pos, mom, charge);
//   //           refParms << curvparms.Qbp(), curvparms.lambda(), curvparms.phi(), curvparms.xT(), curvparms.yT();
//             refParms << curvparms.Qbp(), curvparms.lambda(), curvparms.phi();
//   //           refParms += dxcurv;

          }
          
          // *TODO* better handling of this case?
          if ( (muchargearr[0] + muchargearr[1]) != 0) {
            valid = false;
            break;
          }

          const unsigned int idxplus = muchargearr[0] > 0 ? 0 : 1;
          const unsigned int idxminus = muchargearr[0] > 0 ? 1 : 0;
          
          if (icons == 0) {
          
            Muplus_pt = muarr[idxplus].pt();
            Muplus_eta = muarr[idxplus].eta();
            Muplus_phi = muarr[idxplus].phi();
            
            Muminus_pt = muarr[idxminus].pt();
            Muminus_eta = muarr[idxminus].eta();
            Muminus_phi = muarr[idxminus].phi();
            
            Mupluskin_pt = outparts[idxplus]->currentState().globalMomentum().perp();
            Mupluskin_eta = outparts[idxplus]->currentState().globalMomentum().eta();
            Mupluskin_phi = outparts[idxplus]->currentState().globalMomentum().phi();
            
            Muminuskin_pt = outparts[idxminus]->currentState().globalMomentum().perp();
            Muminuskin_eta = outparts[idxminus]->currentState().globalMomentum().eta();
            Muminuskin_phi = outparts[idxminus]->currentState().globalMomentum().phi();
          }
          else {
            Mupluscons_pt = muarr[idxplus].pt();
            Mupluscons_eta = muarr[idxplus].eta();
            Mupluscons_phi = muarr[idxplus].phi();
            
            Muminuscons_pt = muarr[idxminus].pt();
            Muminuscons_eta = muarr[idxminus].eta();
            Muminuscons_phi = muarr[idxminus].phi();
            
            Mupluskincons_pt = outparts[idxplus]->currentState().globalMomentum().perp();
            Mupluskincons_eta = outparts[idxplus]->currentState().globalMomentum().eta();
            Mupluskincons_phi = outparts[idxplus]->currentState().globalMomentum().phi();
            
            Muminuskincons_pt = outparts[idxminus]->currentState().globalMomentum().perp();
            Muminuskincons_eta = outparts[idxminus]->currentState().globalMomentum().eta();
            Muminuskincons_phi = outparts[idxminus]->currentState().globalMomentum().phi();
          }
          
//           std::cout << "Muplus pt, eta, phi = " << Muplus_pt << ", " << Muplus_eta << ", " << Muplus_phi << std::endl;
//           std::cout << "Muminus pt, eta, phi = " << Muminus_pt << ", " << Muminus_eta << ", " << Muminus_phi << std::endl;
          
          
          MatrixXd covstate =  2.*Cinvd.solve(MatrixXd::Identity(nstateparms,nstateparms));
          
//           Matrix<double, 6, 6> covrefmom;
          covrefmom  = Matrix<double, 6, 6>::Zero();
          
          constexpr std::array<unsigned int, 2> localidxs = {{ 0, 3 }};
          const std::array<unsigned int, 2> globalidxs = {{ trackstateidxarr[0], trackstateidxarr[1] }};
          
          for (unsigned int iidx = 0; iidx < localidxs.size(); ++iidx) {
            for (unsigned int jidx = 0; jidx < localidxs.size(); ++jidx) {
              covrefmom.block<3, 3>(localidxs[iidx], localidxs[jidx]) = covstate.block<3, 3>(globalidxs[iidx], globalidxs[jidx]);
            } 
          }
          
          if (icons == 0) {
          
            Map<Matrix<float, 3, 1>>(Muplus_refParms.data()) = mucurvarr[idxplus].cast<float>();
            Map<Matrix<float, 3, 1>>(Muminus_refParms.data()) = mucurvarr[idxminus].cast<float>();
            
    //         std::cout << "nstateparms = " << nstateparms << std::endl;
    //         std::cout << "dxdparms " << dxdparms.rows() << " " << dxdparms.cols() << std::endl;
            
//             Muplus_jacRef.resize(3*npars);
//             Map<Matrix<float, 3, Dynamic, RowMajor>>(Muplus_jacRef.data(), 3, npars) = dxdparms.block(0, trackstateidxarr[idxplus], npars, 3).transpose().cast<float>();
//             
//             Muminus_jacRef.resize(3*npars);
//             Map<Matrix<float, 3, Dynamic, RowMajor>>(Muminus_jacRef.data(), 3, npars) = dxdparms.block(0, trackstateidxarr[idxminus], npars, 3).transpose().cast<float>();
            
            
//             MatrixXd covstate =  2.*Cinvd.solve(MatrixXd::Identity(nstateparms,nstateparms));
            
//             Matrix<double, 6, 6> covrefmom;
//             
//             constexpr std::array<unsigned int, 2> localidxs = {{ 0, 3 }};
//             const std::array<unsigned int, 2> globalidxs = {{ trackstateidxarr[0], trackstateidxarr[1] }};
//             
//             for (unsigned int iidx = 0; iidx < localidxs.size(); ++iidx) {
//               for (unsigned int jidx = 0; jidx < localidxs.size(); ++jidx) {
//                 covrefmom.block<3, 3>(localidxs[iidx], localidxs[jidx]) = covstate.block<3, 3>(globalidxs[iidx], globalidxs[jidx]);
//               } 
//             }
            
            
            
            const Matrix<double, 1, 6> mjacalt = massJacobianAltD(refftsarr[0], refftsarr[1], mmu);

            
            Jpsi_sigmamass = std::sqrt((mjacalt*covrefmom*mjacalt.transpose())[0]);
            
//             std::cout << "covrefmom" << std::endl;
//             std::cout << covrefmom << std::endl;
//             std::cout << "Jpsi_sigmamass = " << Jpsi_sigmamass << std::endl;
          
          }
  //         
  //         (jacarr[idxplus].topLeftCorner(5, nstateparms)*dxdparms.transpose() + jacarr[idxplus].topRightCorner(5, npars)).cast<float>();
  //         
  //         Muminus_jacRef.resize(3*npars);
  //         Map<Matrix<float, 3, Dynamic, RowMajor>>(Muminus_jacRef.data(), 3, npars) = (jacarr[idxminus].topLeftCorner(5, nstateparms)*dxdparms.transpose() + jacarr[idxminus].topRightCorner(5, npars)).cast<float>();
          
          //TODO fix this
  //         Muplus_jacRef.resize(5*npars);
  //         Map<Matrix<float, 5, Dynamic, RowMajor>>(Muplus_jacRef.data(), 5, npars) = (jacarr[idxplus].topLeftCorner(5, nstateparms)*dxdparms.transpose() + jacarr[idxplus].topRightCorner(5, npars)).cast<float>();
  //         
  //         Muminus_jacRef.resize(5*npars);
  //         Map<Matrix<float, 5, Dynamic, RowMajor>>(Muminus_jacRef.data(), 5, npars) = (jacarr[idxminus].topLeftCorner(5, nstateparms)*dxdparms.transpose() + jacarr[idxminus].topRightCorner(5, npars)).cast<float>();
          
          auto const jpsimom = muarr[0] + muarr[1];
          
          if (icons == 0) {
            Jpsi_pt = jpsimom.pt();
            Jpsi_eta = jpsimom.eta();
            Jpsi_phi = jpsimom.phi();
            Jpsi_mass = jpsimom.mass();
          }
          else {
            Jpsicons_pt = jpsimom.pt();
            Jpsicons_eta = jpsimom.eta();
            Jpsicons_phi = jpsimom.phi();
            Jpsicons_mass = jpsimom.mass();
          }
          
          Muplus_nhits = nhitsarr[idxplus];
          Muplus_nvalid = nvalidarr[idxplus];
          Muplus_nvalidpixel = nvalidpixelarr[idxplus];
          Muplus_nmatchedvalid = nmatchedvalidarr[idxplus];
          Muplus_nambiguousmatchedvalid = nambiguousmatchedvalidarr[idxplus];
          
          Muminus_nhits = nhitsarr[idxminus];
          Muminus_nvalid = nvalidarr[idxminus];
          Muminus_nvalidpixel = nvalidpixelarr[idxminus];
          Muminus_nmatchedvalid = nmatchedvalidarr[idxminus];
          Muminus_nambiguousmatchedvalid = nambiguousmatchedvalidarr[idxminus];
          
          const ROOT::Math::PxPyPzMVector mompluskin(outparts[idxplus]->currentState().globalMomentum().x(),
                                                            outparts[idxplus]->currentState().globalMomentum().y(),
                                                            outparts[idxplus]->currentState().globalMomentum().z(),
                                                            mmu);
          
          const ROOT::Math::PxPyPzMVector momminuskin(outparts[idxminus]->currentState().globalMomentum().x(),
                                                            outparts[idxminus]->currentState().globalMomentum().y(),
                                                            outparts[idxminus]->currentState().globalMomentum().z(),
                                                            mmu);
          
          auto const jpsimomkin = mompluskin + momminuskin;
          
          if (icons == 0) {
            Jpsikin_pt = jpsimomkin.pt();
            Jpsikin_eta = jpsimomkin.eta();
            Jpsikin_phi = jpsimomkin.phi();
            Jpsikin_mass = jpsimomkin.mass();
          }
          else {
            Jpsikincons_pt = jpsimomkin.pt();
            Jpsikincons_eta = jpsimomkin.eta();
            Jpsikincons_phi = jpsimomkin.phi();
            Jpsikincons_mass = jpsimomkin.mass();            
          }
          
          const reco::Candidate *muplusgen = nullptr;
          const reco::Candidate *muminusgen = nullptr;
          
          if (doGen_) {
            for (auto const &genpart : *genPartCollection) {
              if (genpart.status() != 1) {
                continue;
              }
              if (std::abs(genpart.pdgId()) != 13) {
                continue;
              }
              
//               float dRplus = deltaR(genpart.phi(), muarr[idxplus].phi(), genpart.eta(), muarr[idxplus].eta());
              float dRplus = deltaR(genpart, muarr[idxplus]);
              if (dRplus < 0.1 && genpart.charge() > 0) {
                muplusgen = &genpart;
              }
              
//               float dRminus = deltaR(genpart.phi(), muarr[idxminus].phi(), genpart.eta(), muarr[idxminus].eta());
              float dRminus = deltaR(genpart, muarr[idxminus]);
              if (dRminus < 0.1 && genpart.charge() < 0) {
                muminusgen = &genpart;
              }
            }
          }
          
          if (muplusgen != nullptr) {
            Muplusgen_pt = muplusgen->pt();
            Muplusgen_eta = muplusgen->eta();
            Muplusgen_phi = muplusgen->phi();
          }
          else {
            Muplusgen_pt = -99.;
            Muplusgen_eta = -99.;
            Muplusgen_phi = -99.;
          }
          
          if (muminusgen != nullptr) {
            Muminusgen_pt = muminusgen->pt();
            Muminusgen_eta = muminusgen->eta();
            Muminusgen_phi = muminusgen->phi();
          }
          else {
            Muminusgen_pt = -99.;
            Muminusgen_eta = -99.;
            Muminusgen_phi = -99.;
          }
          
          if (muplusgen != nullptr && muminusgen != nullptr) {
            auto const jpsigen = ROOT::Math::PtEtaPhiMVector(muplusgen->pt(), muplusgen->eta(), muplusgen->phi(), mmu) +
                                ROOT::Math::PtEtaPhiMVector(muminusgen->pt(), muminusgen->eta(), muminusgen->phi(), mmu);
            
            Jpsigen_pt = jpsigen.pt();
            Jpsigen_eta = jpsigen.eta();
            Jpsigen_phi = jpsigen.phi();
            Jpsigen_mass = jpsigen.mass();
            
            Jpsigen_x = muplusgen->vx();
            Jpsigen_y = muplusgen->vy();
            Jpsigen_z = muplusgen->vz();
          }
          else {
            Jpsigen_pt = -99.;
            Jpsigen_eta = -99.;
            Jpsigen_phi = -99.;
            Jpsigen_mass = -99.;
            
            Jpsigen_x = -99.;
            Jpsigen_y = -99.;
            Jpsigen_z = -99.;
          }


          Muplus_isMuon = false;
          Muplus_muonLoose = false;
          Muplus_muonMedium = false;
          Muplus_muonTight = false;
          Muplus_muonIsPF = false;
          Muplus_muonIsTracker = false;
          Muplus_muonIsGlobal = false;
          Muplus_muonIsStandalone = false;
          Muplus_muonInnerTrackBest = false;

          Muminus_isMuon = false;
          Muminus_muonLoose = false;
          Muminus_muonMedium = false;
          Muminus_muonTight = false;
          Muminus_muonIsPF = false;
          Muminus_muonIsTracker = false;
          Muminus_muonIsGlobal = false;
          Muminus_muonIsStandalone = false;
          Muminus_muonInnerTrackBest = false;

          if (doMuons_) {
            const reco::Muon *matchedmuonplus = idxplus == 0 ? matchedmuon0 : matchedmuon1;
            const reco::Muon *matchedmuonminus = idxminus == 0 ? matchedmuon0 : matchedmuon1;

            if (matchedmuonplus != nullptr) {
              Muplus_isMuon = true;
              Muplus_muonLoose = matchedmuonplus->passed(reco::Muon::CutBasedIdLoose);
              Muplus_muonMedium = matchedmuonplus->passed(reco::Muon::CutBasedIdMedium);
              Muplus_muonTight = matchedmuonplus->passed(reco::Muon::CutBasedIdTight);
              Muplus_muonIsPF = matchedmuonplus->isPFMuon();
              Muplus_muonIsTracker = matchedmuonplus->isTrackerMuon();
              Muplus_muonIsGlobal = matchedmuonplus->isGlobalMuon();
              Muplus_muonIsStandalone = matchedmuonplus->isStandAloneMuon();
              Muplus_muonInnerTrackBest = matchedmuonplus->muonBestTrackType() == reco::Muon::InnerTrack;
            }

            if (matchedmuonminus != nullptr) {
              Muminus_isMuon = true;
              Muminus_muonLoose = matchedmuonminus->passed(reco::Muon::CutBasedIdLoose);
              Muminus_muonMedium = matchedmuonminus->passed(reco::Muon::CutBasedIdMedium);
              Muminus_muonTight = matchedmuonminus->passed(reco::Muon::CutBasedIdTight);
              Muminus_muonIsPF = matchedmuonminus->isPFMuon();
              Muminus_muonIsTracker = matchedmuonminus->isTrackerMuon();
              Muminus_muonIsGlobal = matchedmuonminus->isGlobalMuon();
              Muminus_muonIsStandalone = matchedmuonminus->isStandAloneMuon();
              Muminus_muonInnerTrackBest = matchedmuonminus->muonBestTrackType() == reco::Muon::InnerTrack;
            }

          }
          
          
      //     const Vector5d dxRef = dxfull.head<5>();
      //     const Matrix5d Cinner = Cinvd.solve(MatrixXd::Identity(nstateparms,nstateparms)).topLeftCorner<5,5>();


          niter = iiter + 1;
          edmval = -deltachisq[0];

          if (std::isnan(edmval) || std::isinf(edmval)) {
            std::cout << "WARNING: invalid parameter update!!!" << " edmval = " << edmval << " deltachisqval = " << deltachisqval << std::endl;
            valid = false;
            break;
          }
          
          if (icons == 0) {
            edmval_cons0 = edmval;
            niter_cons0 = niter;
          }
          
//           std::cout << "icons = " << icons << " iiter = " << iiter << " edmval = " << edmval << " deltachisqval = " << deltachisqval << " chisqval = " << chisqval << std::endl;
//           std::cout << "dxvtx" << std::endl;
          
//           std::cout << "pt0 = " << refftsarr[0].momentum().perp() << " eta0 = " << refftsarr[0].momentum().eta() << " charge0 = " << refftsarr[0].charge() <<  " pt1 = " << refftsarr[1].momentum().perp() << " eta1 = " << refftsarr[1].momentum().eta() << " charge1 = " << refftsarr[1].charge() << std::endl;
//           std::cout << "dxref" << std::endl;
          
//           std::cout << "dxvtx:" << std::endl;
//           std::cout << dxfull.head<3>() << std::endl;
//           std::cout << "dxmom0" << std::endl;
//           std::cout << dxfull.segment<3>(trackstateidxarr[0]) << std::endl;
//           std::cout << "dxmom1" << std::endl;
//           std::cout << dxfull.segment<3>(trackstateidxarr[1]) << std::endl;
//           std::cout << "qop0 = " << refftsarr[0].signedInverseMomentum() + dxfull[trackstateidxarr[0]] << std::endl;
//           std::cout << "qop1 = " << refftsarr[1].signedInverseMomentum() + dxfull[trackstateidxarr[1]] << std::endl;
          
          
//           std::cout << "dx0" << std::endl;
//           std::cout << dxfull.segment(trackstateidxarr[0], trackstateidxarr[1]-trackstateidxarr[0]) << std::endl;
//           std::cout << "dx1" << std::endl;
//           std::cout << dxfull.segment(trackstateidxarr[1], nstateparms - trackstateidxarr[1]) << std::endl;
//           std::cout << dxfull.segment<3>(trackstateidxarr[0]) << std::endl;
//           std::cout << dxfull.segment<3>(trackstateidxarr[1]) << std::endl;
         
//           if (std::abs(deltachisqval)<1e-2) {
//             break;
//           }
          
//           if (iiter > 0 && edmval < 1e-5) {
//             break;
//           }
          
          if (iiter > 0 && dolocalupdate && edmval < 1e-5) {
            break;
          }
          else if (iiter > 0 && !dolocalupdate && std::fabs(deltachisqval)<1e-5) {
            break;
          }
          
//           if (iiter > 1 && std::abs(deltachisq[0])<1e-3) {
//             break;
//           }
      
        }
      
        if (!valid) {
          break;
        }

        auto const& dchisqdx = gradfull.head(nstateparms);
        auto const& dchisqdparms = gradfull.tail(npars);

        auto const& d2chisqdx2 = hessfull.topLeftCorner(nstateparms, nstateparms);
        auto const& d2chisqdxdparms = hessfull.topRightCorner(nstateparms, npars);
        auto const& d2chisqdparms2 = hessfull.bottomRightCorner(npars, npars);

        std::unordered_map<unsigned int, unsigned int> idxmap;

        globalidxvfinal.clear();
        globalidxvfinal.reserve(npars);
        idxmap.reserve(npars);

        for (unsigned int idx : globalidxv) {
          if (!idxmap.count(idx)) {
            idxmap[idx] = globalidxvfinal.size();
            globalidxvfinal.push_back(idx);
          }
        }

        const unsigned int nparsfinal = globalidxvfinal.size();

        VectorXd dchisqdparmsfinal = VectorXd::Zero(nparsfinal);
        MatrixXd d2chisqdxdparmsfinal = MatrixXd::Zero(nstateparms, nparsfinal);
        MatrixXd d2chisqdparms2final = MatrixXd::Zero(nparsfinal, nparsfinal);

        for (unsigned int i = 0; i < npars; ++i) {
          const unsigned int iidx = idxmap.at(globalidxv[i]);
          dchisqdparmsfinal[iidx] += dchisqdparms[i];
          d2chisqdxdparmsfinal.col(iidx) += d2chisqdxdparms.col(i);
          for (unsigned int j = 0; j < npars; ++j) {
            const unsigned int jidx = idxmap.at(globalidxv[j]);
            d2chisqdparms2final(iidx, jidx) += d2chisqdparms2(i, j);
          }
        }

        dxdparms = -Cinvd.solve(d2chisqdxdparmsfinal).transpose();

    //     grad = dchisqdparmsfinal + dxdparms*dchisqdx;
        grad = dchisqdparmsfinal + d2chisqdxdparmsfinal.transpose()*dxfull;
        hess = d2chisqdparms2final + dxdparms*d2chisqdxdparmsfinal;

        if (icons == 0) {
          const unsigned int idxplus = muchargearr[0] > 0 ? 0 : 1;
          const unsigned int idxminus = muchargearr[0] > 0 ? 1 : 0;

          Muplus_jacRef.resize(3*nparsfinal);
          Map<Matrix<float, 3, Dynamic, RowMajor>>(Muplus_jacRef.data(), 3, nparsfinal) = dxdparms.block(0, trackstateidxarr[idxplus], nparsfinal, 3).transpose().cast<float>();

          Muminus_jacRef.resize(3*nparsfinal);
          Map<Matrix<float, 3, Dynamic, RowMajor>>(Muminus_jacRef.data(), 3, nparsfinal) = dxdparms.block(0, trackstateidxarr[idxminus], nparsfinal, 3).transpose().cast<float>();
        }
      
      }
      
      if (!valid) {
        continue;
      }
      
//       std::cout << "gradfull rows cols " << gradfull.rows() << "  " << gradfull.cols() << "nstateparms = " << nstateparms << std::endl;
    
//       auto const& dchisqdx = gradfull.head(nstateparms);
//       auto const& dchisqdparms = gradfull.tail(npars);
//
//       auto const& d2chisqdx2 = hessfull.topLeftCorner(nstateparms, nstateparms);
//       auto const& d2chisqdxdparms = hessfull.topRightCorner(nstateparms, npars);
//       auto const& d2chisqdparms2 = hessfull.bottomRightCorner(npars, npars);
//
//       std::unordered_map<unsigned int, unsigned int> idxmap;
//
//       globalidxvfinal.clear();
//       globalidxvfinal.reserve(npars);
//       idxmap.reserve(npars);
//
//       for (unsigned int idx : globalidxv) {
//         if (!idxmap.count(idx)) {
//           idxmap[idx] = globalidxvfinal.size();
//           globalidxvfinal.push_back(idx);
//         }
//       }
//
      const unsigned int nparsfinal = globalidxvfinal.size();
//
//       VectorXd dchisqdparmsfinal = VectorXd::Zero(nparsfinal);
//       MatrixXd d2chisqdxdparmsfinal = MatrixXd::Zero(nstateparms, nparsfinal);
//       MatrixXd d2chisqdparms2final = MatrixXd::Zero(nparsfinal, nparsfinal);
//
//       for (unsigned int i = 0; i < npars; ++i) {
//         const unsigned int iidx = idxmap.at(globalidxv[i]);
//         dchisqdparmsfinal[iidx] += dchisqdparms[i];
//         d2chisqdxdparmsfinal.col(iidx) += d2chisqdxdparms.col(i);
//         for (unsigned int j = 0; j < npars; ++j) {
//           const unsigned int jidx = idxmap.at(globalidxv[j]);
//           d2chisqdparms2final(iidx, jidx) += d2chisqdparms2(i, j);
//         }
//       }
//
//       dxdparms = -Cinvd.solve(d2chisqdxdparmsfinal).transpose();
//
//   //     grad = dchisqdparmsfinal + dxdparms*dchisqdx;
//       grad = dchisqdparmsfinal + d2chisqdxdparmsfinal.transpose()*dxfull;
//       hess = d2chisqdparms2final + dxdparms*d2chisqdxdparmsfinal;
      
  //     if (debugprintout_) {
  //       std::cout << "dxrefdparms" << std::endl;
  //       std::cout << dxdparms.leftCols<5>() << std::endl;
  //     }
      
//       grad = dchisqdparms + dxdparms*dchisqdx;
      
//       std::cout << "dchisqdparms" << std::endl;
//       std::cout << dchisqdparms.transpose() << std::endl;
//       std::cout << "dxdparms*dchisqdx" << std::endl;
//       std::cout << (dxdparms*dchisqdx).transpose() << std::endl;
//       std::cout << "grad" << std::endl;
//       std::cout << grad.transpose() << std::endl;
      //TODO check the simplification
  //     hess = d2chisqdparms2 + 2.*dxdparms*d2chisqdxdparms + dxdparms*d2chisqdx2*dxdparms.transpose();
//       hess = d2chisqdparms2 + dxdparms*d2chisqdxdparms;
  
//       for (unsigned int iparm = 0; iparm < npars; ++iparm) {
//         if (detidparmsrev[globalidxv[iparm]].first != 7) {
//           hess.row(iparm) *= 0.;
//           hess.col(iparm) *= 0.;
//           hess(iparm, iparm) = 1e6;
//         }
//       }
      
//       SelfAdjointEigenSolver<MatrixXd> es(hess, EigenvaluesOnly);
//       const double condition = es.eigenvalues()[nstateparms-1]/es.eigenvalues()[0];
//       std::cout << "hess eigenvalues:" << std::endl;
//       std::cout << es.eigenvalues().transpose() << std::endl;
//       std::cout << "condition: " << condition << std::endl;
      
//       std::cout << "hess diagonal:" << std::endl;
//       std::cout << hess.diagonal().transpose() << std::endl;
//       
//       assert(es.eigenvalues()[0] > -1e-5);
//       assert(hess.diagonal().minCoeff() > 0.);
      
      nParms = nparsfinal;

      gradv.clear();
      gradv.resize(nparsfinal,0.);
      
      if (fillTrackTree_ && fillGrads_) {
        tree->SetBranchAddress("gradv", gradv.data());
      }
      
      //eigen representation of the underlying vector storage
      Map<VectorXf> gradout(gradv.data(), nparsfinal);

      gradout = grad.cast<float>();
      
        

      gradmax = 0.;
      for (unsigned int i=0; i<nparsfinal; ++i) {
        const float absval = std::abs(grad[i]);
        if (absval>gradmax) {
          gradmax = absval;
        }      
      }
      
      
      hessmax = 0.;
      for (unsigned int i=0; i<nparsfinal; ++i) {
        for (unsigned int j=i; j<nparsfinal; ++j) {
          const unsigned int iidx = globalidxvfinal[i];
          const unsigned int jidx = globalidxvfinal[j];
          
          const float absval = std::abs(hess(i,j));
          if (absval>hessmax) {
            hessmax = absval;
          }
          
        }
        
      }
        
      //fill packed hessian and indices
      const unsigned int nsym = nparsfinal*(1+nparsfinal)/2;
      hesspackedv.clear();    
      hesspackedv.resize(nsym, 0.);
      
      nSym = nsym;
      if (fillTrackTree_ && fillGrads_) {
        tree->SetBranchAddress("hesspackedv", hesspackedv.data());
      }
      
      Map<VectorXf> hesspacked(hesspackedv.data(), nsym);
      const Map<const VectorXu> globalidx(globalidxvfinal.data(), nparsfinal);

      unsigned int packedidx = 0;
      for (unsigned int ipar = 0; ipar < nparsfinal; ++ipar) {
        const unsigned int segmentsize = nparsfinal - ipar;
        hesspacked.segment(packedidx, segmentsize) = hess.block<1, Dynamic>(ipar, ipar, 1, segmentsize).cast<float>();
        packedidx += segmentsize;
      }
      
      
//       assert(globalidxvfinal.size() == (2*Muplus_nhits + 2*Muminus_nhits + 2*Muplus_nvalid + 2*Muminus_nvalid + Muplus_nvalidpixel + Muminus_nvalidpixel));

//       hessv.resize(nparsfinal*nparsfinal);
//       Map<Matrix<float, Dynamic, Dynamic, RowMajor>>(hessv.data(), nparsfinal, nparsfinal) = hess.cast<float>();
      
      if (fillTrackTree_) {
        tree->Fill();
      }
      

      

      
//       const Matrix3d covvtx = Cinvd.solve(MatrixXd::Identity(nstateparms,nstateparms)).topLeftCorner<3,3>();
//       
//       const double covqop0 = Cinvd.solve(MatrixXd::Identity(nstateparms,nstateparms))(trackstateidxarr[0], trackstateidxarr[0]);
//       const double covqop1 = Cinvd.solve(MatrixXd::Identity(nstateparms,nstateparms))(trackstateidxarr[1], trackstateidxarr[1]);
//       
//       const double covqop0kin = outparts[0]->currentState().freeTrajectoryState().curvilinearError().matrix()(0,0);
//       const double covqop1kin = outparts[1]->currentState().freeTrajectoryState().curvilinearError().matrix()(0,0);
//       
// //       Matrix<double, 1, 1> covmass = 2.*massjac*Cinvd.solve(MatrixXd::Identity(nstateparms,nstateparms))*massjac.transpose();
//       
// //       const VectorXd cinvrow0 = Cinvd.solve(MatrixXd::Identity(nstateparms,nstateparms)).row(0).head(nstateparms);
// //       
//       std::cout << "kinfit covariance:" << std::endl;
//       std::cout << dimu_vertex->error().matrix() << std::endl;
//       
//       std::cout << "GBL covariance:" << std::endl;
//       std::cout << 2.*covvtx << std::endl;
//       
//       std::cout << "kinfit qop0 covariance:" << std::endl;
//       std::cout << covqop0kin << std::endl;
//       
//       std::cout << "GBL qop0 covariance:" << std::endl;
//       std::cout << 2.*covqop0 << std::endl;
//       
//       std::cout << "kinfit qop1 covariance:" << std::endl;
//       std::cout << covqop1kin << std::endl;
//       
//       std::cout << "GBL qop1 covariance:" << std::endl;
//       std::cout << 2.*covqop1 << std::endl;
      
//       std::cout << "dqop0 beamline" << std::endl;
//       std::cout << dxfull[trackstateidxarr[0]] << std::endl;
//       std::cout << "dqop0 first layer" << std::endl;
//       std::cout << dxfull[trackstateidxarr[0]+3] << std::endl;
//       std::cout << "dqop0 second layer" << std::endl;
//       std::cout << dxfull[trackstateidxarr[0]+6] << std::endl;
//       
//       std::cout << "dqop1 beamline" << std::endl;
//       std::cout << dxfull[trackstateidxarr[1]] << std::endl;
//       std::cout << "dqop1 first layer" << std::endl;
//       std::cout << dxfull[trackstateidxarr[1]+3] << std::endl;
//       std::cout << "dqop1 second layer" << std::endl;
//       std::cout << dxfull[trackstateidxarr[1]+6] << std::endl;
//       
//       std::cout << "sigmam" << std::endl;
//       std::cout << std::sqrt(covmass[0]) << std::endl;
      
//       
//       std::cout << "cinvrow0" << std::endl;
//       std::cout << cinvrow0 << std::endl;
      
      //TODO restore statejac stuff
//       dxstate = statejac*dxfull;
//       const Vector5d dxRef = dxstate.head<5>();
//       const Matrix5d Cinner = (statejac*Cinvd.solve(MatrixXd::Identity(nstateparms,nstateparms))*statejac.transpose()).topLeftCorner<5,5>();
      
      //TODO fill outputs
      
    }
  }
  
}


DEFINE_FWK_MODULE(ResidualGlobalCorrectionMakerTwoTrackG4e);
