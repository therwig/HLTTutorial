// -*- C++ -*-
//
// Package:    HLTAnalysis/TriggerAnalyzerRAWMiniAOD
// Class:      TriggerAnalyzerRAWMiniAOD
// 
/**\class TriggerAnalyzerRAWMiniAOD TriggerAnalyzerRAWMiniAOD.cc HLTAnalysis/TriggerAnalyzerRAWMiniAOD/plugins/TriggerAnalyzerRAWMiniAOD.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Laurent Thomas
//         Created:  Fri, 24 Mar 2017 04:09:55 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/Common/interface/AssociationMap.h"

#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TLorentzVector.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class TriggerAnalyzerRAWMiniAOD : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit TriggerAnalyzerRAWMiniAOD(const edm::ParameterSet&);
      ~TriggerAnalyzerRAWMiniAOD();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
  bool PassOfflineMuonSelection(const pat::Muon *mu, reco::Vertex::Point PV);

      // ----------member data ---------------------------


  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> trigobjectsMINIAODToken_;
  edm::EDGetTokenT<edm::TriggerResults> trgresultsORIGToken_;
  edm::EDGetTokenT<trigger::TriggerEvent> trigobjectsRAWToken_;
  edm::EDGetTokenT<edm::TriggerResults>  trgresultsHLT2Token_;

  edm::EDGetTokenT<std::vector<pat::Jet> > jet_token;
  edm::EDGetTokenT<std::vector<pat::Muon> > muon_token;
  edm::EDGetTokenT<std::vector<pat::MET> > met_token;
  edm::EDGetTokenT<std::vector<pat::Electron> > electron_token;
  edm::EDGetTokenT<std::vector<reco::Vertex> > PV_token;

  edm::EDGetTokenT<std::vector<reco::GenMET> > genMET_token;
  edm::EDGetTokenT<std::vector<reco::GenJet> > genJet_token;
  edm::EDGetTokenT<std::vector<reco::GenParticle> > genPart_token;

  edm::Service<TFileService> fs;
  
  // histos
  TH1F* h_counts;

  TH1F* h_HLT_PFMET120_PFMHT120_IDTight_vs_met_num;
  TH1F* h_HLT_PFMET120_PFMHT120_IDTight_vs_mht_num;
  TH1F* h_HLT_DoubleMu3_DZ_PFMET50_PFMHT60_vs_met_num;
  TH1F* h_HLT_DoubleMu3_DZ_PFMET50_PFMHT60_vs_mht_num;
  TH1F* h_HLT_DoubleMu3_DZ_PFMET50_PFMHT60_vs_mu2_num;
  TH1F* h_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_vs_mu_num;
  TH1F* h_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_vs_jet_num;
  TH1F* h_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_vs_met_num;
  TH1F* h_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_vs_mht_num;
  TH1F* h_HLT_Mu15_IsoVVVL_PFHT450_vs_mu_num;
  TH1F* h_HLT_Mu15_IsoVVVL_PFHT450_vs_ht_num;
  TH1F* h_HLT_PFHT800_PFMET75_PFMHT75_IDTight_vs_met_num;
  TH1F* h_HLT_PFHT800_PFMET75_PFMHT75_IDTight_vs_mht_num;
  TH1F* h_HLT_PFHT800_PFMET75_PFMHT75_IDTight_vs_ht_num;
  TH1F* h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_mu_num;
  TH1F* h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_jet_num;
  TH1F* h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_deta_num;
  TH1F* h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_mjj_num;
  TH1F* h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_ht_num;
  TH1F* h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_met_num;
  TH1F* h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_mht_num;

  TH1F* h_HLT_PFMET120_PFMHT120_IDTight_vs_met_den;
  TH1F* h_HLT_PFMET120_PFMHT120_IDTight_vs_mht_den;
  TH1F* h_HLT_DoubleMu3_DZ_PFMET50_PFMHT60_vs_met_den;
  TH1F* h_HLT_DoubleMu3_DZ_PFMET50_PFMHT60_vs_mht_den;
  TH1F* h_HLT_DoubleMu3_DZ_PFMET50_PFMHT60_vs_mu2_den;
  TH1F* h_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_vs_mu_den;
  TH1F* h_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_vs_jet_den;
  TH1F* h_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_vs_met_den;
  TH1F* h_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_vs_mht_den;
  TH1F* h_HLT_Mu15_IsoVVVL_PFHT450_vs_mu_den;
  TH1F* h_HLT_Mu15_IsoVVVL_PFHT450_vs_ht_den;
  TH1F* h_HLT_PFHT800_PFMET75_PFMHT75_IDTight_vs_met_den;
  TH1F* h_HLT_PFHT800_PFMET75_PFMHT75_IDTight_vs_mht_den;
  TH1F* h_HLT_PFHT800_PFMET75_PFMHT75_IDTight_vs_ht_den;
  TH1F* h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_mu_den;
  TH1F* h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_jet_den;
  TH1F* h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_deta_den;
  TH1F* h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_mjj_den;
  TH1F* h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_ht_den;
  TH1F* h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_met_den;
  TH1F* h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_mht_den;

  // distributions of all quantities for triggers where effs are measured per-leg (no cuts)
  TH1F* h_met;
  TH1F* h_mht5p0;
  TH1F* h_mht2p5;
  TH1F* h_mu1pt_2p4_high;
  TH1F* h_mu1pt_1p4;
  TH1F* h_mu2pt_1p4;
  TH1F* h_mu1pt_2p4;
  TH1F* h_mu2pt_2p4;
  TH1F* h_j1pt_2p5;
  TH1F* h_j1pt_5p0;
  TH1F* h_j2pt_5p0;
  TH1F* h_ht5p0;
  TH1F* h_ht2p5;
  TH1F* h_ht5p0_high;
  TH1F* h_ht2p5_high;
  TH1F* h_mjj;
  TH1F* h_deta;

};

//
// constructors and destructor
//
TriggerAnalyzerRAWMiniAOD::TriggerAnalyzerRAWMiniAOD(const edm::ParameterSet& iConfig)

{
  trigobjectsMINIAODToken_ = consumes<pat::TriggerObjectStandAloneCollection>( edm::InputTag("slimmedPatTrigger"));
  trigobjectsRAWToken_=consumes<trigger::TriggerEvent>(edm::InputTag("hltTriggerSummaryAOD::HLT2"));  

  trgresultsORIGToken_= consumes<edm::TriggerResults>( edm::InputTag("TriggerResults::HLT") );
  trgresultsHLT2Token_= consumes<edm::TriggerResults>( edm::InputTag("TriggerResults::HLT2") );  

  jet_token = consumes< std::vector<pat::Jet> >(edm::InputTag("slimmedJets") );
  muon_token = consumes<std::vector<pat::Muon> >(edm::InputTag("slimmedMuons") );
  met_token = consumes<std::vector<pat::MET> >(edm::InputTag("slimmedMETs") );
  electron_token = consumes<std::vector<pat::Electron> >(edm::InputTag("slimmedElectrons") );
  PV_token = consumes<std::vector<reco::Vertex> > (edm::InputTag("offlineSlimmedPrimaryVertices"));
  genMET_token = consumes<std::vector<reco::GenMET> >(edm::InputTag("genMetTrue") );
  genJet_token = consumes<std::vector<reco::GenJet> >(edm::InputTag("ak4GenJetsNoNu") );
  genPart_token = consumes<std::vector<reco::GenParticle> >(edm::InputTag("genParticles") );

  //now do what ever initialization is needed
  usesResource("TFileService");

  h_counts = fs->make<TH1F>("h_counts","",20,0,20);
  h_counts->SetCanExtend(TH1::kAllAxes);
  h_counts->Fill("dummy", 1 );

  h_HLT_PFMET120_PFMHT120_IDTight_vs_met_num                                  = fs->make<TH1F>("h_HLT_PFMET120_PFMHT120_IDTight_vs_met_num", "", 80, 0, 400);
  h_HLT_PFMET120_PFMHT120_IDTight_vs_mht_num                                  = fs->make<TH1F>("h_HLT_PFMET120_PFMHT120_IDTight_vs_mht_num", "", 80, 0, 400);
  h_HLT_DoubleMu3_DZ_PFMET50_PFMHT60_vs_met_num                               = fs->make<TH1F>("h_HLT_DoubleMu3_DZ_PFMET50_PFMHT60_vs_met_num", "", 80, 0, 400);
  h_HLT_DoubleMu3_DZ_PFMET50_PFMHT60_vs_mht_num                               = fs->make<TH1F>("h_HLT_DoubleMu3_DZ_PFMET50_PFMHT60_vs_mht_num", "", 80, 0, 400);
  h_HLT_DoubleMu3_DZ_PFMET50_PFMHT60_vs_mu2_num                               = fs->make<TH1F>("h_HLT_DoubleMu3_DZ_PFMET50_PFMHT60_vs_mu2_num", "", 80, 0, 8);
  h_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_vs_mu_num              = fs->make<TH1F>("h_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_vs_mu_num", "", 80, 0, 8);
  h_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_vs_jet_num             = fs->make<TH1F>("h_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_vs_jet_num", "", 80, 0, 400);
  h_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_vs_met_num             = fs->make<TH1F>("h_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_vs_met_num", "", 80, 0, 400);
  h_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_vs_mht_num             = fs->make<TH1F>("h_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_vs_mht_num", "", 80, 0, 400);
  h_HLT_Mu15_IsoVVVL_PFHT450_vs_mu_num                                        = fs->make<TH1F>("h_HLT_Mu15_IsoVVVL_PFHT450_vs_mu_num", "", 120, 0, 30);
  h_HLT_Mu15_IsoVVVL_PFHT450_vs_ht_num                                        = fs->make<TH1F>("h_HLT_Mu15_IsoVVVL_PFHT450_vs_ht_num", "", 80, 0, 800);
  h_HLT_PFHT800_PFMET75_PFMHT75_IDTight_vs_met_num                                    = fs->make<TH1F>("h_HLT_PFHT800_PFMET75_PFMHT75_IDTight_vs_met_num", "", 80, 0, 400);
  h_HLT_PFHT800_PFMET75_PFMHT75_IDTight_vs_mht_num                                    = fs->make<TH1F>("h_HLT_PFHT800_PFMET75_PFMHT75_IDTight_vs_mht_num", "", 80, 0, 400);
  h_HLT_PFHT800_PFMET75_PFMHT75_IDTight_vs_ht_num                                     = fs->make<TH1F>("h_HLT_PFHT800_PFMET75_PFMHT75_IDTight_vs_ht_num", "", 120, 0, 1200);
  h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_mu_num   = fs->make<TH1F>("h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_mu_num", "", 80, 0, 8);
  h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_jet_num  = fs->make<TH1F>("h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_jet_num", "", 80, 0, 400);
  h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_deta_num = fs->make<TH1F>("h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_deta_num", "", 50, 0, 5);
  h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_mjj_num  = fs->make<TH1F>("h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_mjj_num", "", 120, 0, 1200);
  h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_ht_num   = fs->make<TH1F>("h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_ht_num", "", 120, 0, 600);
  h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_met_num  = fs->make<TH1F>("h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_met_num", "", 80, 0, 400);
  h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_mht_num  = fs->make<TH1F>("h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_mht_num", "", 80, 0, 400);

  h_HLT_PFMET120_PFMHT120_IDTight_vs_met_den                                  = fs->make<TH1F>("h_HLT_PFMET120_PFMHT120_IDTight_vs_met_den", "", 80, 0, 400);
  h_HLT_PFMET120_PFMHT120_IDTight_vs_mht_den                                  = fs->make<TH1F>("h_HLT_PFMET120_PFMHT120_IDTight_vs_mht_den", "", 80, 0, 400);
  h_HLT_DoubleMu3_DZ_PFMET50_PFMHT60_vs_met_den                               = fs->make<TH1F>("h_HLT_DoubleMu3_DZ_PFMET50_PFMHT60_vs_met_den", "", 80, 0, 400);
  h_HLT_DoubleMu3_DZ_PFMET50_PFMHT60_vs_mht_den                               = fs->make<TH1F>("h_HLT_DoubleMu3_DZ_PFMET50_PFMHT60_vs_mht_den", "", 80, 0, 400);
  h_HLT_DoubleMu3_DZ_PFMET50_PFMHT60_vs_mu2_den                               = fs->make<TH1F>("h_HLT_DoubleMu3_DZ_PFMET50_PFMHT60_vs_mu2_den", "", 80, 0, 8);
  h_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_vs_mu_den              = fs->make<TH1F>("h_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_vs_mu_den", "", 80, 0, 8);
  h_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_vs_jet_den             = fs->make<TH1F>("h_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_vs_jet_den", "", 80, 0, 400);
  h_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_vs_met_den             = fs->make<TH1F>("h_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_vs_met_den", "", 80, 0, 400);
  h_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_vs_mht_den             = fs->make<TH1F>("h_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_vs_mht_den", "", 80, 0, 400);
  h_HLT_Mu15_IsoVVVL_PFHT450_vs_mu_den                                        = fs->make<TH1F>("h_HLT_Mu15_IsoVVVL_PFHT450_vs_mu_den", "", 120, 0, 30);
  h_HLT_Mu15_IsoVVVL_PFHT450_vs_ht_den                                        = fs->make<TH1F>("h_HLT_Mu15_IsoVVVL_PFHT450_vs_ht_den", "", 80, 0, 800);
  h_HLT_PFHT800_PFMET75_PFMHT75_IDTight_vs_met_den                                    = fs->make<TH1F>("h_HLT_PFHT800_PFMET75_PFMHT75_IDTight_vs_met_den", "", 80, 0, 400);
  h_HLT_PFHT800_PFMET75_PFMHT75_IDTight_vs_mht_den                                    = fs->make<TH1F>("h_HLT_PFHT800_PFMET75_PFMHT75_IDTight_vs_mht_den", "", 80, 0, 400);
  h_HLT_PFHT800_PFMET75_PFMHT75_IDTight_vs_ht_den                                     = fs->make<TH1F>("h_HLT_PFHT800_PFMET75_PFMHT75_IDTight_vs_ht_den", "", 120, 0, 1200);
  h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_mu_den   = fs->make<TH1F>("h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_mu_den", "", 80, 0, 8);
  h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_jet_den  = fs->make<TH1F>("h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_jet_den", "", 80, 0, 400);
  h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_deta_den = fs->make<TH1F>("h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_deta_den", "", 50, 0, 5);
  h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_mjj_den  = fs->make<TH1F>("h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_mjj_den", "", 120, 0, 1200);
  h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_ht_den   = fs->make<TH1F>("h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_ht_den", "", 120, 0, 600);
  h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_met_den  = fs->make<TH1F>("h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_met_den", "", 80, 0, 400);
  h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_mht_den  = fs->make<TH1F>("h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_mht_den", "", 80, 0, 400);

  h_met            = fs->make<TH1F>("h_met" , "", 80, 0, 400);
  h_mht5p0         = fs->make<TH1F>("h_mht5p0" , "", 80, 0, 400);
  h_mht2p5         = fs->make<TH1F>("h_mht2p5" , "", 80, 0, 400);
  h_mu1pt_2p4_high = fs->make<TH1F>("h_mu1pt_2p4_high" , "", 50, 0, 50);
  h_mu1pt_1p4      = fs->make<TH1F>("h_mu1pt_1p4" , "", 40, 0, 10);
  h_mu2pt_1p4      = fs->make<TH1F>("h_mu2pt_1p4" , "", 40, 0, 10);
  h_mu1pt_2p4      = fs->make<TH1F>("h_mu1pt_2p4" , "", 40, 0, 10);
  h_mu2pt_2p4      = fs->make<TH1F>("h_mu2pt_2p4" , "", 40, 0, 10);
  h_j1pt_2p5       = fs->make<TH1F>("h_j1pt_2p5" , "", 80, 0, 400);
  h_j1pt_5p0       = fs->make<TH1F>("h_j1pt_5p0" , "", 80, 0, 400);
  h_j2pt_5p0       = fs->make<TH1F>("h_j2pt_5p0" , "", 80, 0, 400);
  h_ht5p0          = fs->make<TH1F>("h_ht5p0" , "", 80, 0, 800);
  h_ht2p5          = fs->make<TH1F>("h_ht2p5" , "", 80, 0, 800);
  h_ht5p0_high     = fs->make<TH1F>("h_ht5p0_high" , "", 80, 0, 1200);
  h_ht2p5_high     = fs->make<TH1F>("h_ht2p5_high" , "", 80, 0, 1200);
  h_mjj            = fs->make<TH1F>("h_mjj" , "", 120, 0, 1200);
  h_deta           = fs->make<TH1F>("h_deta", "", 50, 0, 5); 

}


TriggerAnalyzerRAWMiniAOD::~TriggerAnalyzerRAWMiniAOD()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TriggerAnalyzerRAWMiniAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   using namespace edm;
   using namespace reco;
   using namespace std;



   // // ****************Part 1. Accessing some trigger information ************* 
   // Single muon triggers
   bool passHLT_IsoMu24(false);
   // MET triggers (vs gen met)
   bool passHLT_PFMET110_PFMHT110_IDTight(false);
   bool passHLT_PFMET120_PFMHT120_IDTight(false);
   bool passHLT_PFMET130_PFMHT130_IDTight(false);
   // 2mu + MET (vs met, mu2)
   bool passHLT_DoubleMu3_DZ_PFMET50_PFMHT60(false);
   // mu + MET + VBF (vs met, ptmu, mjj, ptj2, deta, ht)
   bool passHLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60(false);
   // mu + MET + ISR (vs met, ptmu, ptj1)   
   bool passHLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight(false);
   // mu+HT (vs mu, ht)
   bool passHLT_Mu15_IsoVVVL_PFHT450(false);
   // HT + MET (vs ht, met)
   bool passHLT_PFHT800_PFMET75_PFMHT75_IDTight(false);


   // //Accessing trigger bits:
   // //This works in both RAW, AOD or MINIAOD 
   // //Here we access the decision provided by the HLT (i.e. original trigger step). 

   edm::Handle<edm::TriggerResults> trigResults;
   iEvent.getByToken(trgresultsORIGToken_, trigResults);
   if( !trigResults.failedToGet() ) {
     int N_Triggers = trigResults->size();
     const edm::TriggerNames & trigName = iEvent.triggerNames(*trigResults);

     for( int i_Trig = 0; i_Trig < N_Triggers; ++i_Trig ) {
       if (trigResults.product()->accept(i_Trig)) {
         TString TrigPath =trigName.triggerName(i_Trig);
         //	 cout << "Passed path: " << TrigPath<<endl;
         if(TrigPath.Index("HLT_IsoMu24_v") >=0) passHLT_IsoMu24=true; 
         //Notice the special syntax: since the path version can change during data taking one only looks for the string "HLT_IsoMu24_v"
       }
     }
   }
   //cout << "passHLT_IsoMu24 = " << passHLT_IsoMu24 << endl;

   //Exercise 1: 
   //Clone and *then* modify the code above in order to save the decision of your customized HLT menu in the booleans passHLT_Mu3_PFJet200DeepCSV_1p59, passHLT_Mu3_L1SingleJet180, passHLT_PFJet200DeepCSV_1p59
   //Do not directly edit the code above as you will also need the use the original HLT_IsoMu24 decision later on.

   //Solution: 
   edm::Handle<edm::TriggerResults> trigResultsHLT2;
   iEvent.getByToken(trgresultsHLT2Token_, trigResultsHLT2);
   if( !trigResultsHLT2.failedToGet() ) {
     int N_Triggers = trigResultsHLT2->size();
     const edm::TriggerNames & trigName = iEvent.triggerNames(*trigResultsHLT2);

     for( int i_Trig = 0; i_Trig < N_Triggers; ++i_Trig ) {
       if (trigResultsHLT2.product()->accept(i_Trig)) {
         TString TrigPath =trigName.triggerName(i_Trig);
         //	 cout << "Passed path: " << TrigPath<<endl;
         if(TrigPath.Index("HLT_PFMET110_PFMHT110_IDTight_v") >=0) passHLT_PFMET110_PFMHT110_IDTight=true; 
         if(TrigPath.Index("HLT_PFMET120_PFMHT120_IDTight_v") >=0) passHLT_PFMET120_PFMHT120_IDTight=true; 
         if(TrigPath.Index("HLT_PFMET130_PFMHT130_IDTight_v") >=0) passHLT_PFMET130_PFMHT130_IDTight=true; 
         // if(TrigPath.Index("HLT_Mu3_PFJet200DeepCSV_1p59_v") >=0) passHLT_Mu3_PFJet200DeepCSV_1p59=true; 
         // if(TrigPath.Index("HLT_Mu3_L1SingleJet180_v") >=0)passHLT_Mu3_L1SingleJet180=true;
         // if(TrigPath.Index("HLT_PFJet200DeepCSV_1p59_v") >=0)passHLT_PFJet200DeepCSV_1p59=true;
         //Notice the special syntax: since the path version can change during data taking one only looks for the string "HLT_IsoMu24_v"
         if(TrigPath.Index("HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60") >=0) passHLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60 = true;
         if(TrigPath.Index("HLT_DoubleMu3_DZ_PFMET50_PFMHT60") >=0)                              passHLT_DoubleMu3_DZ_PFMET50_PFMHT60 = true;
         if(TrigPath.Index("HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight") >=0)            passHLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight = true;
         if(TrigPath.Index("HLT_Mu15_IsoVVVL_PFHT450") >=0)                                      passHLT_Mu15_IsoVVVL_PFHT450 = true;
         if(TrigPath.Index("HLT_PFHT800_PFMET75_PFMHT75_IDTight") >=0)                                   passHLT_PFHT800_PFMET75_PFMHT75_IDTight = true;
       }
     }
   }


   //
   // "Offline" (gen) quantities to plot the trigger performance versus
   float met=0;
   float mu1pt_1p4=0;
   float mu2pt_1p4=0;
   float mu1pt_2p4=0;
   float mu2pt_2p4=0;
   float mjj=0;
   float deta=0;
   float ht5p0=0;
   float ht2p5=0;
   float j1pt_5p0=0;
   float j2pt_5p0=0;
   float j1pt_2p5=0;
   float mht5p0x=0;
   float mht5p0y=0;
   float mht5p0=0;
   float mht2p5x=0;
   float mht2p5y=0;
   float mht2p5=0;


   edm::Handle< std::vector<reco::GenMET> > genMETs;
   iEvent.getByToken(genMET_token, genMETs );
   met = (*genMETs).front().pt();

   edm::Handle< std::vector<reco::GenJet> > genJets;
   iEvent.getByToken(genJet_token, genJets );
   
   std::vector<const reco::GenJet*> jets5p0;
   std::vector<const reco::GenJet*> jets2p5;
   for( const auto &j : *genJets ){
       if(fabs(j.eta()) < 5.0 && j.pt()>20) jets5p0.push_back(&j);
       if(fabs(j.eta()) < 2.5 && j.pt()>20) jets2p5.push_back(&j);
   }
   std::sort(jets5p0.begin(), jets5p0.end(), [](const reco::GenJet* i, const reco::GenJet* j) {return (i->pt() > j->pt());});
   std::sort(jets2p5.begin(), jets2p5.end(), [](const reco::GenJet* i, const reco::GenJet* j) {return (i->pt() > j->pt());});
   if (jets5p0.size()) j1pt_5p0 = jets5p0.front()->pt();
   if (jets2p5.size()) j1pt_2p5 = jets2p5.front()->pt();

   // used for vbf only
   float m;
   for(unsigned int i=0; i<jets5p0.size(); i++){
       for(unsigned int j=i+1; j<jets5p0.size(); j++){
           auto j1 = jets5p0[i];
           auto j2 = jets5p0[j];
           TLorentzVector v1, v2;
           v1.SetPtEtaPhiM(j1->pt(), j1->eta(), j1->phi(), 0.);
           v2.SetPtEtaPhiM(j2->pt(), j2->eta(), j2->phi(), 0.);
           m = (v1+v2).M();
           // prefer the largest mjj pair
           if (m > mjj){
               mjj = m;
               deta = fabs(j1->eta()-j2->eta());
               j2pt_5p0 = j2->pt();
           }
           // vbf selection (a little below threshold)
           // if (j2->pt() < 60) continue;
           // if( fabs(j1->eta()-j2->eta()) < 3.0) continue;
           // if (m<500) continue;
       }
   }
   
   // no pt cuts here at the moment
   for( const auto &j : *genJets ){
       if(fabs(j.eta()) < 5.0){
           ht5p0 += j.pt();
           mht5p0x += j.px();
           mht5p0y += j.py();
       }
       if(fabs(j.eta()) < 2.5){
           ht2p5 += j.pt();
           mht2p5x += j.px();
           mht2p5y += j.py();
       }
   }
   mht5p0 = hypot(mht5p0x, mht5p0y);
   mht2p5 = hypot(mht2p5x, mht2p5y);

   edm::Handle< std::vector<reco::GenParticle> > genParts;
   iEvent.getByToken(genPart_token, genParts );
   std::vector<const reco::GenParticle*> muons2p4;
   std::vector<const reco::GenParticle*> muons1p4;
   for(const auto & p : *genParts){
       if (p.status()!=1) continue;
       if ( abs(p.pdgId())==13 ){
           if (p.pt()<3) continue;
           if (fabs(p.eta())< 2.4) muons2p4.push_back(&p);
           if (fabs(p.eta())< 1.4) muons1p4.push_back(&p);
       }
   }   
   std::sort(muons2p4.begin(), muons2p4.end(), [](const reco::GenParticle* i, const reco::GenParticle* j) {return (i->pt() > j->pt());});
   std::sort(muons1p4.begin(), muons1p4.end(), [](const reco::GenParticle* i, const reco::GenParticle* j) {return (i->pt() > j->pt());});

   // mu1 pt
   if( muons1p4.size() ) mu1pt_1p4 = muons1p4.front()->pt();
   if( muons2p4.size() ) mu1pt_2p4 = muons2p4.front()->pt();

   for(unsigned int i=0; i<muons2p4.size(); i++){
       for(unsigned int j=i+1; j<muons2p4.size(); j++){
           if( deltaR(muons2p4[i]->eta(), muons2p4[i]->phi(), muons2p4[j]->eta(), muons2p4[j]->phi()) < 0.3 ) continue;
           if( muons2p4[j]->pt() > mu2pt_2p4) mu2pt_2p4 = muons2p4[j]->pt();
       }
   }

   for(unsigned int i=0; i<muons1p4.size(); i++){
       for(unsigned int j=i+1; j<muons1p4.size(); j++){
           if( deltaR(muons1p4[i]->eta(), muons1p4[i]->phi(), muons1p4[j]->eta(), muons1p4[j]->phi()) < 0.3 ) continue;
           if( muons1p4[j]->pt() > mu2pt_1p4) mu2pt_1p4 = muons1p4[j]->pt();
       }
   }

   // triggers for histograms
   h_counts->Fill("all", 1 );
   h_counts->Fill("HLT_IsoMu24",                                                  passHLT_IsoMu24                                                  );
   h_counts->Fill("HLT_PFMET110_PFMHT110_IDTight",                                passHLT_PFMET110_PFMHT110_IDTight                                );
   h_counts->Fill("HLT_PFMET120_PFMHT120_IDTight",                                passHLT_PFMET120_PFMHT120_IDTight                                );
   h_counts->Fill("HLT_PFMET130_PFMHT130_IDTight",                                passHLT_PFMET130_PFMHT130_IDTight                                );
   h_counts->Fill("HLT_DoubleMu3_DZ_PFMET50_PFMHT60",                             passHLT_DoubleMu3_DZ_PFMET50_PFMHT60                             );
   h_counts->Fill("HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight",           passHLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight           );
   h_counts->Fill("HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60",passHLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60);
   h_counts->Fill("HLT_Mu15_IsoVVVL_PFHT450",                                     passHLT_Mu15_IsoVVVL_PFHT450                                     );
   h_counts->Fill("HLT_PFHT800_PFMET75_PFMHT75_IDTight",                                  passHLT_PFHT800_PFMET75_PFMHT75_IDTight                                  );
   

   // event distributions
   h_met->Fill(met);
   h_mht5p0->Fill(mht5p0);
   h_mht2p5->Fill(mht2p5);
   h_mu1pt_2p4_high->Fill(mu1pt_2p4);
   h_mu1pt_1p4->Fill(mu1pt_1p4);
   h_mu2pt_1p4->Fill(mu2pt_1p4);
   h_mu1pt_2p4->Fill(mu1pt_2p4);
   h_mu2pt_2p4->Fill(mu2pt_2p4);
   h_j1pt_2p5->Fill(j1pt_2p5);
   h_j1pt_5p0->Fill(j1pt_5p0);
   h_j2pt_5p0->Fill(j2pt_5p0);
   h_ht5p0->Fill(ht5p0);
   h_ht2p5->Fill(ht2p5);
   h_ht5p0_high->Fill(ht5p0);
   h_ht2p5_high->Fill(ht2p5);
   h_mjj->Fill(mjj);
   h_deta->Fill(deta);


   // met
   h_HLT_PFMET120_PFMHT120_IDTight_vs_met_den->Fill(met) ;
   if (passHLT_PFMET120_PFMHT120_IDTight) h_HLT_PFMET120_PFMHT120_IDTight_vs_met_num->Fill(met) ;
   h_HLT_PFMET120_PFMHT120_IDTight_vs_mht_den->Fill(mht2p5) ;
   if (passHLT_PFMET120_PFMHT120_IDTight) h_HLT_PFMET120_PFMHT120_IDTight_vs_mht_num->Fill(mht2p5) ;

   // 2mu + MET
   if(mu2pt_1p4 > 6){
       h_HLT_DoubleMu3_DZ_PFMET50_PFMHT60_vs_met_den->Fill(met) ;
       if (passHLT_DoubleMu3_DZ_PFMET50_PFMHT60) h_HLT_DoubleMu3_DZ_PFMET50_PFMHT60_vs_met_num->Fill(met) ;
       h_HLT_DoubleMu3_DZ_PFMET50_PFMHT60_vs_mht_den->Fill(mht2p5) ;
       if (passHLT_DoubleMu3_DZ_PFMET50_PFMHT60) h_HLT_DoubleMu3_DZ_PFMET50_PFMHT60_vs_mht_num->Fill(mht2p5) ;
   }
   if(met > 125){
       h_HLT_DoubleMu3_DZ_PFMET50_PFMHT60_vs_mu2_den->Fill(mu2pt_1p4) ;
       if (passHLT_DoubleMu3_DZ_PFMET50_PFMHT60) h_HLT_DoubleMu3_DZ_PFMET50_PFMHT60_vs_mu2_num->Fill(mu2pt_1p4) ;
   }

   // mu + ISR + MET
   if(j1pt_2p5 > 150 && met > 200 && mht2p5 > 200){
       h_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_vs_mu_den->Fill(mu1pt_1p4) ;
       if (passHLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight) h_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_vs_mu_num->Fill(mu1pt_1p4) ;
   }
   if(mu1pt_1p4 > 6 && met > 200 && mht2p5 > 200){
       h_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_vs_jet_den->Fill(j1pt_2p5) ;
       if (passHLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight) h_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_vs_jet_num->Fill(j1pt_2p5) ;
   }
   if(mu1pt_1p4 > 6 && j1pt_2p5 > 150){
       h_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_vs_met_den->Fill(met) ;
       if (passHLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight) h_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_vs_met_num->Fill(met) ;
   }
   if(mu1pt_1p4 > 6 && j1pt_2p5 > 150){
       h_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_vs_mht_den->Fill(mht2p5) ;
       if (passHLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight) h_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_vs_mht_num->Fill(mht2p5) ;
   }

   // mu + HT
   if(ht2p5>500){
       h_HLT_Mu15_IsoVVVL_PFHT450_vs_mu_den->Fill(mu1pt_1p4) ;
       if (passHLT_Mu15_IsoVVVL_PFHT450) h_HLT_Mu15_IsoVVVL_PFHT450_vs_mu_num->Fill(mu1pt_1p4) ;
   }
   if(mu1pt_1p4 > 17){
       h_HLT_Mu15_IsoVVVL_PFHT450_vs_ht_den->Fill(ht2p5) ;
       if (passHLT_Mu15_IsoVVVL_PFHT450) h_HLT_Mu15_IsoVVVL_PFHT450_vs_ht_num->Fill(ht2p5) ;
   }

   // HT + MET
   if(ht2p5 > 1000){
       h_HLT_PFHT800_PFMET75_PFMHT75_IDTight_vs_met_den->Fill(met) ;
       if (passHLT_PFHT800_PFMET75_PFMHT75_IDTight) h_HLT_PFHT800_PFMET75_PFMHT75_IDTight_vs_met_num->Fill(met) ;
       h_HLT_PFHT800_PFMET75_PFMHT75_IDTight_vs_mht_den->Fill(mht2p5) ;
       if (passHLT_PFHT800_PFMET75_PFMHT75_IDTight) h_HLT_PFHT800_PFMET75_PFMHT75_IDTight_vs_mht_num->Fill(mht2p5) ;
   }
   if(met > 150){
       h_HLT_PFHT800_PFMET75_PFMHT75_IDTight_vs_ht_den->Fill(ht2p5) ;
       if (passHLT_PFHT800_PFMET75_PFMHT75_IDTight) h_HLT_PFHT800_PFMET75_PFMHT75_IDTight_vs_ht_num->Fill(ht2p5) ;
   }

   // mu + VBF + MET
   if(j2pt_5p0 > 120 && deta > 3.6 && mjj > 1000 && ht5p0 > 400 && (met > 125 || mht5p0 > 125)){
       h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_mu_den->Fill(mu1pt_1p4) ;
       if (passHLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60) h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_mu_num->Fill(mu1pt_1p4) ;
   }
   if(mu1pt_1p4 > 6 && deta > 3.6 && mjj > 1000 && ht5p0 > 400 && (met > 125 || mht5p0 > 125)){
       h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_jet_den->Fill(j2pt_5p0) ;
       if (passHLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60) h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_jet_num->Fill(j2pt_5p0) ;
   }
   if(mu1pt_1p4 > 6 && j2pt_5p0 > 120 && mjj > 1000 && ht5p0 > 400 && (met > 125 || mht5p0 > 125)){
       h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_deta_den->Fill(deta) ;
       if (passHLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60) h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_deta_num->Fill(deta) ;
   }
   if(mu1pt_1p4 > 6 && j2pt_5p0 > 120 && deta > 3.6 && ht5p0 > 400 && (met > 125 || mht5p0 > 125)){
       h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_mjj_den->Fill(mjj) ;
       if (passHLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60) h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_mjj_num->Fill(mjj) ;
   }
   if(mu1pt_1p4 > 6 && j2pt_5p0 > 120 && deta > 3.6 && mjj > 1000 && (met > 125 || mht5p0 > 125)){
       h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_ht_den->Fill(ht5p0) ;
       if (passHLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60) h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_ht_num->Fill(ht5p0) ;
   }
   if(mu1pt_1p4 > 6 && j2pt_5p0 > 120 && deta > 3.6 && mjj > 1000 && ht5p0 > 400){
       h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_met_den->Fill(met) ;
       if (passHLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60) h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_met_num->Fill(met) ;
   }
   if(mu1pt_1p4 > 6 && j2pt_5p0 > 120 && deta > 3.6 && mjj > 1000 && ht5p0 > 400){
       h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_mht_den->Fill(mht5p0) ;
       if (passHLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60) h_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_vs_mht_num->Fill(mht5p0) ;
   }

   
   // h_DoubleMu3_DZ_PFMET50_PFMHT60_vs_met_den->Fill( met );
   // //if(nMuPlus5 && nMuMinus5) h_DoubleMu3_DZ_PFMET50_PFMHT60_vs_met_num->Fill( met );

   // h_met110_vs_met_den->Fill( met );
   // if(passHLT_PFMET110_PFMHT110_IDTight) h_met110_vs_met_num->Fill( met );
   // h_met120_vs_met_den->Fill( met );
   // if(passHLT_PFMET120_PFMHT120_IDTight) h_met120_vs_met_num->Fill( met );
   // h_met130_vs_met_den->Fill( met );
   // if(passHLT_PFMET130_PFMHT130_IDTight) h_met130_vs_met_num->Fill( met );

}


// ------------ method called once each job just before starting event loop  ------------
void 
TriggerAnalyzerRAWMiniAOD::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TriggerAnalyzerRAWMiniAOD::endJob() 
{
  h_counts->LabelsDeflate("X");
  h_counts->LabelsOption("v");
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TriggerAnalyzerRAWMiniAOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

bool TriggerAnalyzerRAWMiniAOD::PassOfflineMuonSelection(const pat::Muon *mu, reco::Vertex::Point PV){
  if ( !(mu->isGlobalMuon() || mu->isTrackerMuon() )) return false;
  if ( !(mu->isPFMuon()) ) return false;
  const reco::TrackRef innerTrack = mu->innerTrack();
  if( innerTrack.isNull() )return false;
  
  bool goodGlb =  mu->isGlobalMuon() &&  mu->globalTrack()->normalizedChi2() < 3
    &&  mu->combinedQuality().chi2LocalPosition < 12  && mu->combinedQuality().trkKink < 20;
  bool good =  mu->innerTrack()->validFraction() >= 0.8 &&  mu->segmentCompatibility() >= (goodGlb ? 0.303 : 0.451)  ;
  
  if(!good) return false;
  if(TMath::Abs(innerTrack->dxy(PV)) >0.1 ) return false;
  if(TMath::Abs(innerTrack->dz(PV)) >0.1 ) return false;  
  
  double chargedHadronIso = mu->pfIsolationR03().sumChargedHadronPt;
  double neutralHadronIso = mu->pfIsolationR03().sumNeutralHadronEt;
  double photonIso = mu->pfIsolationR03().sumPhotonEt;
  
  double beta = mu->pfIsolationR03().sumPUPt;
  double pfRelIsoMu  = ( chargedHadronIso + TMath::Max ( 0.0 ,neutralHadronIso + photonIso - 0.5 * beta ) )/mu->pt() ;
  
  if(pfRelIsoMu >0.4) return false;
  return true;
}









// bool TriggerAnalyzerRAWMiniAOD::PassOfflineElectronSelection(const pat::Electron * ele, reco::Vertex::Point PV){
//   const reco::GsfTrackRef gsfTrack = ele->gsfTrack();
//   if (!gsfTrack.isNonnull()) return false;
//   if( TMath::Abs(gsfTrack->dxy(PV)) > 0.05  )  return false;
//   if( TMath::Abs(gsfTrack->dz(PV)) > 0.2  )  return false;
//   if(TMath::Abs(ele->superCluster()->eta()) >2.5) return false;
//   else if( TMath::Abs(ele->superCluster()->eta()) < 1.479  ) {
//     if( TMath::Abs(ele->full5x5_sigmaIetaIeta()) > 0.0103 ) return  false;
//     if( TMath::Abs(ele->deltaEtaSuperClusterTrackAtVtx()) > 0.0105 ) return  false;
//     if( TMath::Abs(ele->deltaPhiSuperClusterTrackAtVtx()) > 0.115  ) return  false;
//     if( TMath::Abs(ele->hadronicOverEm())  > 0.104  ) return  false;
//     if( TMath::Abs(1.0/ele->ecalEnergy() - ele->eSuperClusterOverP()/ele->ecalEnergy() )>0.102 ) return  false;
//     if( TMath::Abs(gsfTrack->dxy(PV)) > 0.0261) return false;
//     if( TMath::Abs(gsfTrack->dz(PV)) > 0.41) return false;
//     if(gsfTrack->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS)>2) return  false;
//   }
//   else {
//     if( TMath::Abs(ele->full5x5_sigmaIetaIeta()) > 0.0301 ) return  false;
//     if( TMath::Abs(ele->deltaEtaSuperClusterTrackAtVtx()) > 0.00814 ) return  false;
//     if( TMath::Abs(ele->deltaPhiSuperClusterTrackAtVtx()) > 0.182  ) return  false;
//     if( TMath::Abs(ele->hadronicOverEm())  > 0.0897  ) return  false;
//     if( TMath::Abs(1.0/ele->ecalEnergy() - ele->eSuperClusterOverP()/ele->ecalEnergy() )>0.126 ) return  false;
//     if( TMath::Abs(gsfTrack->dxy(PV)) > 0.0118) return false;
//     if( TMath::Abs(gsfTrack->dz(PV)) > 0.822) return false;
//     if(gsfTrack->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS)>1) return  false;


//   }

//   double iso = (ele->pfIsolationVariables().sumChargedHadronPt 
// 		+ TMath::Max(0.0, ele->pfIsolationVariables().sumNeutralHadronEt + ele->pfIsolationVariables().sumPhotonEt - 0.5*ele->pfIsolationVariables().sumPUPt ) 
// 		) /ele->pt() ; 
//   if(iso>0.2) return false; 
//   return true;



// }


// bool TriggerAnalyzerRAWMiniAOD::RecoHLTMatching(const edm::Event& iEvent, double recoeta, double recophi, std::string filtername, double dRmatching){
//   //In the next few lines one loops over all the trigger objects (corresponding to a given filter) and check whether one of them matches the reco object under study                                       
//   edm::Handle<edm::TriggerResults> trigResults;
//   iEvent.getByToken(trgresultsORIGToken_, trigResults);

//   edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
//   iEvent.getByToken(trigobjectsMINIAODToken_, triggerObjects);

//   const edm::TriggerNames &names = iEvent.triggerNames(*trigResults);
//   for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
//     obj.unpackFilterLabels(iEvent,*trigResults);
//     obj.unpackPathNames(names);
//     for (unsigned h = 0; h < obj.filterLabels().size(); ++h){
//       std::string myfillabl=obj.filterLabels()[h];
//       if( myfillabl.find(filtername)!=std::string::npos   && deltaR(recoeta,recophi, obj.eta(),obj.phi())<dRmatching ) return true;
//     }
//   }

//   return false;
// }




// double TriggerAnalyzerRAWMiniAOD::VarStudied( const edm::Event& iEvent, double recoeta, double recophi,
// 					      edm::EDGetTokenT<edm::AssociationMap<edm::OneToValue<std::vector<reco::RecoEcalCandidate>, float > > > varToken_,  edm::EDGetTokenT<trigger::TriggerFilterObjectWithRefs> candToken_,   bool  dividebyE, bool dividebyEt, double dRmatching ){

//   double thevar = 0.;

//   //Inspired from http://cmslxr.fnal.gov/source/HLTrigger/Egamma/src/HLTGenericFilter.cc        
//   edm::Handle<trigger::TriggerFilterObjectWithRefs> PrevFilterOutput;
//   iEvent.getByToken (candToken_, PrevFilterOutput);

//   edm::Handle<edm::AssociationMap<edm::OneToValue<std::vector<reco::RecoEcalCandidate>, float > > > depMap;
//   iEvent.getByToken (varToken_,depMap);

//   std::vector<edm::Ref<std::vector<reco::RecoEcalCandidate> > > recoCands;


//   if(PrevFilterOutput.isValid()&&  depMap.isValid() ){

//     PrevFilterOutput->getObjects(trigger::TriggerCluster, recoCands);
//     if(recoCands.empty())PrevFilterOutput->getObjects(trigger::TriggerPhoton, recoCands);

//     double dRmin = dRmatching;
//     for (unsigned int i=0; i<recoCands.size(); i++) {
//       edm::Ref<std::vector<reco::RecoEcalCandidate> > ref = recoCands[i];
//       typename edm::AssociationMap<edm::OneToValue<std::vector<reco::RecoEcalCandidate>, float > >::const_iterator mapi = (*depMap).find( ref );
//       float vali = mapi->val;
//       float EtaSC = ref->eta();
//       float PhiSC = ref->phi();

//       if(deltaR(recoeta,recophi,EtaSC,PhiSC ) > dRmin )continue;
//       dRmin = deltaR(recoeta,recophi,EtaSC,PhiSC ) ;
//       float energy = ref->superCluster()->energy();
//       float et = ref->superCluster()->energy() * sin (2*atan(exp(-ref->eta())));
//       thevar = (double) vali;
//       if(dividebyE)thevar = (double)vali/energy;
//       if(dividebyEt)thevar =(double) vali/et;

//     }
//   }
//   return thevar;
// }



//define this as a plug-in
DEFINE_FWK_MODULE(TriggerAnalyzerRAWMiniAOD);
