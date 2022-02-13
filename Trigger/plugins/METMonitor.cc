#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMOffline/Trigger/plugins/TriggerDQMBase.h"
#include "CommonTools/TriggerUtils/interface/GenericTriggerEventFlag.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaPhi.h"

class METMonitor : public DQMEDAnalyzer, public TriggerDQMBase {
public:
  typedef dqm::reco::MonitorElement MonitorElement;
  typedef dqm::reco::DQMStore DQMStore;

  METMonitor(const edm::ParameterSet&);
  ~METMonitor() throw() override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

protected:
  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&) override;
  void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) override;

private:
  const std::string folderName_;

  const bool requireValidHLTPaths_;
  bool hltPathsAreValid_;
  const bool enableFullMonitoring_;//=== added the flag==

  edm::InputTag metInputTag_;
  edm::InputTag jetInputTag_;
  edm::InputTag eleInputTag_;
  edm::InputTag muoInputTag_;
  edm::InputTag vtxInputTag_;

  edm::EDGetTokenT<reco::PFMETCollection> metToken_;
  edm::EDGetTokenT<reco::PFJetCollection> jetToken_;
  edm::EDGetTokenT<reco::GsfElectronCollection> eleToken_;
  edm::EDGetTokenT<reco::MuonCollection> muoToken_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;

  std::vector<double> met_variable_binning_;
  std::vector<double> jetptBinning_; //==== added jetPt binning
  MEbinning met_binning_;
  MEbinning ls_binning_;

  ObjME metME_;
  ObjME metME_variableBinning_;
  ObjME metVsLS_;
  ObjME metPhiME_;
  ObjME deltaphimetj1ME_;
  ObjME deltaphij1j2ME_;
  //========= added =========
  ObjME metEtaME_;  
  ObjME phij1ME_;
  ObjME phij2ME_;
  ObjME jetPhiME_;
  ObjME etaj1ME_;
  ObjME etaj2ME_; 
  ObjME jetEtaME_; 
  ObjME jetpt1ME_;
  ObjME jetpt2ME_;
  ObjME jetPtME_;   
  ObjME nJetsME_; 
  //=============================

  std::unique_ptr<GenericTriggerEventFlag> num_genTriggerEventFlag_;
  std::unique_ptr<GenericTriggerEventFlag> den_genTriggerEventFlag_;

  StringCutObjectSelector<reco::MET, true> metSelection_;
  StringCutObjectSelector<reco::PFJet, true> jetSelection_;
  StringCutObjectSelector<reco::GsfElectron, true> eleSelection_;
  StringCutObjectSelector<reco::Muon, true> muoSelection_;

  unsigned njets_;
  unsigned nelectrons_;
  unsigned nmuons_;

  static constexpr double MAX_PHI = 3.2;
  static constexpr int N_PHI = 64;
  static constexpr MEbinning phi_binning_{N_PHI, -MAX_PHI, MAX_PHI};
  //============ added ===============
  // Define Eta Bining
  static constexpr double MAX_ETA = 3.0;
  static constexpr int N_ETA = 30;
  static constexpr MEbinning eta_binning{N_ETA, -MAX_ETA, MAX_ETA}; 
  static constexpr MEbinning jet_eta_binning{50, -5, 5};  //eta binning for inclusive jet eta
  //Define nJets Binning general selection Pt>0
  static constexpr int MIN_NJETS = 0;
  static constexpr int MAX_NJETS = 200;
  static constexpr int N_BIN_NJETS = 200;
  static constexpr MEbinning nJets_binning{N_BIN_NJETS, MIN_NJETS, MAX_NJETS};
  
  ///==============================

  std::vector<bool> warningPrinted4token_;
};

METMonitor::METMonitor(const edm::ParameterSet& iConfig)
    : folderName_(iConfig.getParameter<std::string>("FolderName")),
      requireValidHLTPaths_(iConfig.getParameter<bool>("requireValidHLTPaths")),
      hltPathsAreValid_(false),
      enableFullMonitoring_(iConfig.getParameter<bool>("enableFullMonitoring")), //===added=====
      metInputTag_(iConfig.getParameter<edm::InputTag>("met")),
      jetInputTag_(iConfig.getParameter<edm::InputTag>("jets")),
      eleInputTag_(iConfig.getParameter<edm::InputTag>("electrons")),
      muoInputTag_(iConfig.getParameter<edm::InputTag>("muons")),
      vtxInputTag_(iConfig.getParameter<edm::InputTag>("vertices")),
      metToken_(consumes<reco::PFMETCollection>(metInputTag_)),
      jetToken_(mayConsume<reco::PFJetCollection>(jetInputTag_)),
      eleToken_(mayConsume<reco::GsfElectronCollection>(eleInputTag_)),
      muoToken_(mayConsume<reco::MuonCollection>(muoInputTag_)),
      vtxToken_(mayConsume<reco::VertexCollection>(vtxInputTag_)),
      met_variable_binning_(
          iConfig.getParameter<edm::ParameterSet>("histoPSet").getParameter<std::vector<double> >("metBinning")),
      jetptBinning_(
          iConfig.getParameter<edm::ParameterSet>("histoPSet").getParameter<std::vector<double> >("jetptBinning")), //==== added (jet pt binning)
      met_binning_(getHistoPSet(
          iConfig.getParameter<edm::ParameterSet>("histoPSet").getParameter<edm::ParameterSet>("metPSet"))),
      ls_binning_(
          getHistoPSet(iConfig.getParameter<edm::ParameterSet>("histoPSet").getParameter<edm::ParameterSet>("lsPSet"))),
      num_genTriggerEventFlag_(new GenericTriggerEventFlag(
          iConfig.getParameter<edm::ParameterSet>("numGenericTriggerEventPSet"), consumesCollector(), *this)),
      den_genTriggerEventFlag_(new GenericTriggerEventFlag(
          iConfig.getParameter<edm::ParameterSet>("denGenericTriggerEventPSet"), consumesCollector(), *this)),
      metSelection_(iConfig.getParameter<std::string>("metSelection")),
      jetSelection_(iConfig.getParameter<std::string>("jetSelection")),
      eleSelection_(iConfig.getParameter<std::string>("eleSelection")),
      muoSelection_(iConfig.getParameter<std::string>("muoSelection")),
      njets_(iConfig.getParameter<unsigned>("njets")),
      nelectrons_(iConfig.getParameter<unsigned>("nelectrons")),
      nmuons_(iConfig.getParameter<unsigned>("nmuons")) {
  // this vector has to be alligned to the the number of Tokens accessed by this module
  warningPrinted4token_.push_back(false);  // PFMETCollection
  warningPrinted4token_.push_back(false);  // JetCollection
  warningPrinted4token_.push_back(false);  // GsfElectronCollection
  warningPrinted4token_.push_back(false);  // MuonCollection
  warningPrinted4token_.push_back(false);  // VertexCollection
}

METMonitor::~METMonitor() throw() {
  if (num_genTriggerEventFlag_) {
    num_genTriggerEventFlag_.reset();
  }
  if (den_genTriggerEventFlag_) {
    den_genTriggerEventFlag_.reset();
  }
}

void METMonitor::bookHistograms(DQMStore::IBooker& ibooker, edm::Run const& iRun, edm::EventSetup const& iSetup) {
  // Initialize the GenericTriggerEventFlag
  if (num_genTriggerEventFlag_ && num_genTriggerEventFlag_->on()) {
    num_genTriggerEventFlag_->initRun(iRun, iSetup);
  }
  if (den_genTriggerEventFlag_ && den_genTriggerEventFlag_->on()) {
    den_genTriggerEventFlag_->initRun(iRun, iSetup);
  }

  // check if every HLT path specified in numerator and denominator has a valid match in the HLT Menu
  hltPathsAreValid_ = (num_genTriggerEventFlag_ && den_genTriggerEventFlag_ && num_genTriggerEventFlag_->on() &&
                       den_genTriggerEventFlag_->on() && num_genTriggerEventFlag_->allHLTPathsAreValid() &&
                       den_genTriggerEventFlag_->allHLTPathsAreValid());

  // if valid HLT paths are required,
  // create DQM outputs only if all paths are valid
  if (requireValidHLTPaths_ and (not hltPathsAreValid_)) {
    return;
  }

  std::string histname, histtitle;

  std::string currentFolder = folderName_;
  ibooker.setCurrentFolder(currentFolder);

  histname = "deltaphi_metjet1";
  histtitle = "DPHI_METJ1";
  bookME(ibooker, deltaphimetj1ME_, histname, histtitle, phi_binning_.nbins, phi_binning_.xmin, phi_binning_.xmax);
  setMETitle(deltaphimetj1ME_, "delta phi (met, j1)", "events / 0.1 rad");

  histname = "deltaphi_jet1jet2";
  histtitle = "DPHI_J1J2";
  bookME(ibooker, deltaphij1j2ME_, histname, histtitle, phi_binning_.nbins, phi_binning_.xmin, phi_binning_.xmax);
  setMETitle(deltaphij1j2ME_, "delta phi (j1, j2)", "events / 0.1 rad");

  histname = "met";
  histtitle = "PFMET";
  bookME(ibooker, metME_, histname, histtitle, met_binning_.nbins, met_binning_.xmin, met_binning_.xmax);
  setMETitle(metME_, "PF MET [GeV]", "events / [GeV]");

  histname = "met_variable";
  histtitle = "PFMET";
  bookME(ibooker, metME_variableBinning_, histname, histtitle, met_variable_binning_);
  setMETitle(metME_variableBinning_, "PF MET [GeV]", "events / [GeV]");

  histname = "metVsLS";
  histtitle = "PFMET vs LS";
  bookME(ibooker,
         metVsLS_,
         histname,
         histtitle,
         ls_binning_.nbins,
         ls_binning_.xmin,
         ls_binning_.xmax,
         met_binning_.xmin,
         met_binning_.xmax);
  setMETitle(metVsLS_, "LS", "PF MET [GeV]");

  histname = "metPhi";
  histtitle = "PFMET phi";
  bookME(ibooker, metPhiME_, histname, histtitle, phi_binning_.nbins, phi_binning_.xmin, phi_binning_.xmax);
  setMETitle(metPhiME_, "PF MET #phi", "events / 0.1 rad");

  //=================== added ================
  histname = "metEta";
  histtitle = "PFMET eta";
  bookME(ibooker,
     metEtaME_,
     histname,
     histtitle,
     eta_binning.nbins,
     eta_binning.xmin,
     eta_binning.xmax);
  setMETitle(metEtaME_, "PFMET #eta", "events");

  histname = "jetEta";
  histtitle = "Jet eta";
  bookME(ibooker,
     jetEtaME_,
     histname,
     histtitle,
     jet_eta_binning.nbins,
     jet_eta_binning.xmin,
     jet_eta_binning.xmax);
  setMETitle(jetEtaME_, "Jet #eta", "events");

  histname = "jetPhi";
  histtitle = "Jet phi";
  bookME(ibooker,
     jetPhiME_,
     histname,
     histtitle,
     phi_binning_.nbins,
     phi_binning_.xmin,
     phi_binning_.xmax);
  setMETitle(jetPhiME_, "Jet #phi", "events / 0.1 rad");

  histname = "jetPt";
  histtitle = "Jet Pt";
  bookME(ibooker, jetPtME_, histname, histtitle, jetptBinning_);
  setMETitle(jetPtME_, "jet Pt[GeV]", "events");
           
  histname = "NJets";
  histtitle = "number of Jets";
  bookME(ibooker,
     nJetsME_,
     histname,
     histtitle,
     nJets_binning.nbins,
     nJets_binning.xmin,
     nJets_binning.xmax);
  setMETitle(nJetsME_, "number of Jets", "events"); 

  //check the flag 
  if (!enableFullMonitoring_) {
        return;
   }
 
  histname = "jetEta1";
  histtitle = "leading Jet eta";
  bookME(ibooker,
     etaj1ME_,
     histname,
     histtitle,
     eta_binning.nbins,
     eta_binning.xmin,
     eta_binning.xmax);
  setMETitle(etaj1ME_, "Jet_1 #eta", "events");

      histname = "jetEta2";
      histtitle = "subleading Jet eta";
      bookME(ibooker,
         etaj2ME_,
         histname,
         histtitle,
         eta_binning.nbins,
         eta_binning.xmin,
         eta_binning.xmax);
      setMETitle(etaj2ME_, "Jet_2 #eta", "events");
   
      histname = "jetPhi_1";
      histtitle = "leading Jet phi";
      bookME(ibooker,
         phij1ME_,
         histname,
         histtitle,
         phi_binning_.nbins,
         phi_binning_.xmin,
         phi_binning_.xmax);
      setMETitle(phij1ME_, "Jet_1 #phi", "events / 0.1 rad");

      histname = "jetPhi_2";
      histtitle = "subleading Jet phi";
      bookME(ibooker,
         phij2ME_,
         histname,
         histtitle,
         phi_binning_.nbins,
         phi_binning_.xmin,
         phi_binning_.xmax);
      setMETitle(phij2ME_, "Jet_2 #phi", "events / 0.1 rad");  
  
      histname = "jetPt_1";
      histtitle = "leading Jet Pt";
      bookME(ibooker, jetpt1ME_, histname, histtitle, jetptBinning_);
      setMETitle(jetpt1ME_, "Pt_1 [GeV]", "events");

      histname = "jetPt_2";
      histtitle = "second leading Jet Pt";
      bookME(ibooker, jetpt2ME_, histname, histtitle, jetptBinning_);
      setMETitle(jetpt2ME_, "Pt_2 [GeV]", "events");
      //====================================================================== 

}

void METMonitor::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {
  // if valid HLT paths are required,
  // analyze event only if all paths are valid
  if (requireValidHLTPaths_ and (not hltPathsAreValid_)) {
    return;
  }

  // Filter out events if Trigger Filtering is requested
  if (den_genTriggerEventFlag_->on() && !den_genTriggerEventFlag_->accept(iEvent, iSetup))
    return;

  edm::Handle<reco::PFMETCollection> metHandle;
  iEvent.getByToken(metToken_, metHandle);
  if (!metHandle.isValid()) {
    if (!warningPrinted4token_[0]) {
      edm::LogWarning("METMonitor") << "skipping events because the collection " << metInputTag_.label().c_str()
                                    << " is not available";
      warningPrinted4token_[0] = true;
    }
    return;
  }
  reco::PFMET pfmet = metHandle->front();
  if (!metSelection_(pfmet))
    return;

  float met = pfmet.pt();
  float phi = pfmet.phi();
  float eta = pfmet.eta();

  std::vector<reco::PFJet> jets;
  edm::Handle<reco::PFJetCollection> jetHandle;
  iEvent.getByToken(jetToken_, jetHandle);
  if (jetHandle.isValid()) {
    if (jetHandle->size() < njets_)
      return;
    for (auto const& j : *jetHandle) {
      if (jetSelection_(j)) {
        jets.push_back(j);
      }
    }
  } else {
    if (!warningPrinted4token_[1]) {
      if (jetInputTag_.label().empty())
        edm::LogWarning("METMonitor") << "JetCollection not set";
      else
        edm::LogWarning("METMonitor") << "skipping events because the collection " << jetInputTag_.label().c_str()
                                      << " is not available";
      warningPrinted4token_[1] = true;
    }
    // if Handle is not valid, because the InputTag has been mis-configured, then skip the event
    if (!jetInputTag_.label().empty())
      return;
  }
  float deltaPhi_met_j1 = 10.0;
  float deltaPhi_j1_j2 = 10.0;
  //========= added ===========
  int nJets = jets.size();   
  float Pt_J1 = 0.0;
  float Pt_J2 = 0.0;
  float Phi_J1 = 10.0;
  float Phi_J2 = 10.0;
  float Eta_J1 = 10.0;
  float Eta_J2 = 10.0;  
 
 //==========================
  if (!jets.empty()){
    deltaPhi_met_j1 = fabs(deltaPhi(pfmet.phi(), jets[0].phi()));
    //====== added==============
    if(enableFullMonitoring_){  //check the flag
       Pt_J1  = jets[0].pt();
       Phi_J1 = jets[0].phi();
       Eta_J1 = jets[0].p4().eta(); 
     } 
    //======================= 
   }     
  if (jets.size() >= 2){
    deltaPhi_j1_j2 = fabs(deltaPhi(jets[0].phi(), jets[1].phi()));
    //========= added =========
    if(enableFullMonitoring_){  //check the flag
       Pt_J2 = jets[1].pt();
       Phi_J2 = jets[1].phi();
       Eta_J2 = jets[1].p4().eta();  
     }
    //==================== 
   }
  std::vector<reco::GsfElectron> electrons;
  edm::Handle<reco::GsfElectronCollection> eleHandle;
  iEvent.getByToken(eleToken_, eleHandle);
  if (eleHandle.isValid()) {
    if (eleHandle->size() < nelectrons_)
      return;
    for (auto const& e : *eleHandle) {
      if (eleSelection_(e))
        electrons.push_back(e);
    }
    if (electrons.size() < nelectrons_)
      return;
  } else {
    if (!warningPrinted4token_[2]) {
      warningPrinted4token_[2] = true;
      if (eleInputTag_.label().empty())
        edm::LogWarning("METMonitor") << "GsfElectronCollection not set";
      else
        edm::LogWarning("METMonitor") << "skipping events because the collection " << eleInputTag_.label().c_str()
                                      << " is not available";
    }
    if (!eleInputTag_.label().empty())
      return;
  }

  reco::Vertex vtx;
  edm::Handle<reco::VertexCollection> vtxHandle;
  iEvent.getByToken(vtxToken_, vtxHandle);
  if (vtxHandle.isValid()) {
    for (auto const& v : *vtxHandle) {
      bool isFake = v.isFake();

      if (!isFake) {
        vtx = v;
        break;
      }
    }
  } else {
    if (!warningPrinted4token_[3]) {
      warningPrinted4token_[3] = true;
      if (vtxInputTag_.label().empty())
        edm::LogWarning("METMonitor") << "VertexCollection is not set";
      else
        edm::LogWarning("METMonitor") << "skipping events because the collection " << vtxInputTag_.label().c_str()
                                      << " is not available";
    }
    if (!vtxInputTag_.label().empty())
      return;
  }

  std::vector<reco::Muon> muons;
  edm::Handle<reco::MuonCollection> muoHandle;
  iEvent.getByToken(muoToken_, muoHandle);
  if (muoHandle.isValid()) {
    if (muoHandle->size() < nmuons_)
      return;
    for (auto const& m : *muoHandle) {
      bool pass = m.isGlobalMuon() && m.isPFMuon() && m.globalTrack()->normalizedChi2() < 10. &&
                  m.globalTrack()->hitPattern().numberOfValidMuonHits() > 0 && m.numberOfMatchedStations() > 1 &&
                  fabs(m.muonBestTrack()->dxy(vtx.position())) < 0.2 &&
                  fabs(m.muonBestTrack()->dz(vtx.position())) < 0.5 &&
                  m.innerTrack()->hitPattern().numberOfValidPixelHits() > 0 &&
                  m.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5;
      if (muoSelection_(m) && pass)
        muons.push_back(m);
    }
    if (muons.size() < nmuons_)
      return;
  } else {
    if (!warningPrinted4token_[4]) {
      warningPrinted4token_[4] = true;
      if (muoInputTag_.label().empty())
        edm::LogWarning("METMonitor") << "MuonCollection not set";
      else
        edm::LogWarning("METMonitor") << "skipping events because the collection " << muoInputTag_.label().c_str()
                                      << " is not available";
    }
    if (!muoInputTag_.label().empty())
      return;
  }

  // filling histograms (denominator)
  metME_.denominator->Fill(met);
  metME_variableBinning_.denominator->Fill(met);
  metPhiME_.denominator->Fill(phi);
  deltaphimetj1ME_.denominator->Fill(deltaPhi_met_j1);
  deltaphij1j2ME_.denominator->Fill(deltaPhi_j1_j2); 

  const int ls = iEvent.id().luminosityBlock();
  metVsLS_.denominator->Fill(ls, met);
 
  //=========== added ============
  metEtaME_.denominator->Fill(eta);
  nJetsME_.denominator->Fill(nJets);
  if (!jets.empty()){
    for(int ii=0; ii<nJets; ii++){
       jetPtME_.denominator->Fill(jets[ii].pt());
       jetPhiME_.denominator->Fill(jets[ii].phi());
       jetEtaME_.denominator->Fill(jets[ii].p4().eta());
    }
  }
 if(enableFullMonitoring_){  //===check the flag
      jetpt1ME_.denominator->Fill(Pt_J1);
      jetpt2ME_.denominator->Fill(Pt_J2);       
      phij1ME_.denominator->Fill(Phi_J1);
      phij2ME_.denominator->Fill(Phi_J2);
      etaj1ME_.denominator->Fill(Eta_J1);
      etaj2ME_.denominator->Fill(Eta_J2);  
   }
  //========================= 
 

  
  // applying selection for numerator
  if (num_genTriggerEventFlag_->on() && !num_genTriggerEventFlag_->accept(iEvent, iSetup))
    return;

  // filling histograms (num_genTriggerEventFlag_)
  metME_.numerator->Fill(met);
  metME_variableBinning_.numerator->Fill(met);
  metPhiME_.numerator->Fill(phi);
  metVsLS_.numerator->Fill(ls, met);
  deltaphimetj1ME_.numerator->Fill(deltaPhi_met_j1);
  deltaphij1j2ME_.numerator->Fill(deltaPhi_j1_j2);

   //============ added ===========
   metEtaME_.numerator->Fill(eta); 
   nJetsME_.numerator->Fill(nJets);    
   if (!jets.empty()){
     for(int kk=0; kk<nJets; kk++){
       jetPtME_.numerator->Fill(jets[kk].pt());
       jetPhiME_.numerator->Fill(jets[kk].phi());
       jetEtaME_.numerator->Fill(jets[kk].p4().eta());
     }
    }
   if(enableFullMonitoring_){  //===check the flag
      jetpt1ME_.numerator->Fill(Pt_J1);
      jetpt2ME_.numerator->Fill(Pt_J2);       
      phij1ME_.numerator->Fill(Phi_J1);
      phij2ME_.numerator->Fill(Phi_J2);
      etaj1ME_.numerator->Fill(Eta_J1);
      etaj2ME_.numerator->Fill(Eta_J2);       
     }
    //======================== 
}

void METMonitor::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("FolderName", "HLT/MET");
  desc.add<bool>("requireValidHLTPaths", true);
  //========== add a flag to enable more or less outputs in each trigger =====
  desc.add<bool>("enableFullMonitoring", true);
  //==============================================

  desc.add<edm::InputTag>("met", edm::InputTag("pfMet"));
  desc.add<edm::InputTag>("jets", edm::InputTag("ak4PFJetsCHS"));
  desc.add<edm::InputTag>("electrons", edm::InputTag("gedGsfElectrons"));
  desc.add<edm::InputTag>("muons", edm::InputTag("muons"));
  desc.add<edm::InputTag>("vertices", edm::InputTag("offlinePrimaryVertices"));
  desc.add<std::string>("metSelection", "pt > 0");
  desc.add<std::string>("jetSelection", "pt > 0");
  desc.add<std::string>("eleSelection", "pt > 0");
  desc.add<std::string>("muoSelection", "pt > 0");
  desc.add<unsigned>("njets", 0);
  desc.add<unsigned>("nelectrons", 0);
  desc.add<unsigned>("nmuons", 0);

  edm::ParameterSetDescription genericTriggerEventPSet;
  genericTriggerEventPSet.add<bool>("andOr");
  genericTriggerEventPSet.add<edm::InputTag>("dcsInputTag", edm::InputTag("scalersRawToDigi"));
  genericTriggerEventPSet.add<std::vector<int> >("dcsPartitions", {});
  genericTriggerEventPSet.add<bool>("andOrDcs", false);
  genericTriggerEventPSet.add<bool>("errorReplyDcs", true);
  genericTriggerEventPSet.add<std::string>("dbLabel", "");
  genericTriggerEventPSet.add<bool>("andOrHlt", true);
  genericTriggerEventPSet.add<edm::InputTag>("hltInputTag", edm::InputTag("TriggerResults::HLT"));
  genericTriggerEventPSet.add<std::vector<std::string> >("hltPaths", {});
  genericTriggerEventPSet.add<std::string>("hltDBKey", "");
  genericTriggerEventPSet.add<bool>("errorReplyHlt", false);
  genericTriggerEventPSet.add<unsigned int>("verbosityLevel", 1);

  desc.add<edm::ParameterSetDescription>("numGenericTriggerEventPSet", genericTriggerEventPSet);
  desc.add<edm::ParameterSetDescription>("denGenericTriggerEventPSet", genericTriggerEventPSet);

  edm::ParameterSetDescription histoPSet;
  edm::ParameterSetDescription metPSet;
  fillHistoPSetDescription(metPSet);
  histoPSet.add<edm::ParameterSetDescription>("metPSet", metPSet);
  std::vector<double> bins = {0.,   20.,  40.,  60.,  80.,  90.,  100., 110., 120., 130., 140., 150., 160.,
                              170., 180., 190., 200., 220., 240., 260., 280., 300., 350., 400., 450., 1000.};
  histoPSet.add<std::vector<double> >("metBinning", bins);
  //======= added
  std::vector<double> bins_ = {
      0.,   20.,  40.,  60.,  80.,  90.,  100., 110., 120., 130., 140., 150., 160.,
      170., 180., 190., 200., 220., 240., 260., 280., 300., 350., 400., 450., 1000.};  // Jet pT Binning
  histoPSet.add<std::vector<double> >("jetptBinning", bins_);
  //=======

  edm::ParameterSetDescription lsPSet;
  fillHistoLSPSetDescription(lsPSet);
  histoPSet.add<edm::ParameterSetDescription>("lsPSet", lsPSet);

  desc.add<edm::ParameterSetDescription>("histoPSet", histoPSet);

  descriptions.add("metMonitoring", desc);
}

DEFINE_FWK_MODULE(METMonitor);
