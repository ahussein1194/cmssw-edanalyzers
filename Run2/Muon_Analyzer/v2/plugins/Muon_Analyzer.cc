// -*- C++ -*-
//
// Package:    my_EDAnalyzers/Muon_Analyzer
// Class:      Muon_Analyzer
//
/**\class Muon_Analyzer Muon_Analyzer.cc my_EDAnalyzers/Muon_Analyzer/plugins/Muon_Analyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author: Ahmed Hussein
//         Created:  Thu, 15 Dec 2022 08:40:07 GMT
////////////////////////////////////////////////////////////////////////////////
#include "Muon_Analyzer.h"

// Constructor.
Muon_Analyzer::Muon_Analyzer(const edm::ParameterSet& iConfig):
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  testBool_(iConfig.getParameter<bool>("testBool"))

{
   // Now do what ever initialization is needed.
   edm::Service<TFileService> fs;
   mtree = fs->make<TTree>("Events", "Events");

   // Create the TTree TBranches.
   mtree->Branch("numbermuon",&num_muon);
   mtree->GetBranch("numbermuon")->SetTitle("number of muons");
   mtree->Branch("muon_e",&muon_e);
   mtree->GetBranch("muon_e")->SetTitle("muon energy");
   mtree->Branch("muon_pt",&muon_pt);
   mtree->GetBranch("muon_pt")->SetTitle("muon transverse momentum");
   mtree->Branch("muon_px",&muon_px);
   mtree->GetBranch("muon_px")->SetTitle("muon momentum x-component");
   mtree->Branch("muon_py",&muon_py);
   mtree->GetBranch("muon_py")->SetTitle("muon momentum y-component");
   mtree->Branch("muon_pz",&muon_pz);
   mtree->GetBranch("muon_pz")->SetTitle("muon momentum z-component");
   mtree->Branch("muon_eta",&muon_eta);
   mtree->GetBranch("muon_eta")->SetTitle("muon pseudorapidity");
   mtree->Branch("muon_phi",&muon_phi);
   mtree->GetBranch("muon_phi")->SetTitle("muon polar angle");
   mtree->Branch("muon_ch",&muon_ch);
   mtree->GetBranch("muon_ch")->SetTitle("muon charge");
   mtree->Branch("muon_isLoose",&muon_isLoose);
   mtree->GetBranch("muon_isLoose")->SetTitle("muon tagged loose");
   mtree->Branch("muon_isMedium",&muon_isMedium);
   mtree->GetBranch("muon_isMedium")->SetTitle("muon tagged medium");
   mtree->Branch("muon_isTight",&muon_isTight);
   mtree->GetBranch("muon_isTight")->SetTitle("muon tagged tight");
   mtree->Branch("muon_isSoft",&muon_isSoft);
   mtree->GetBranch("muon_isSoft")->SetTitle("muon tagged soft");
   mtree->Branch("muon_isHighPt",&muon_isHighPt);
   mtree->GetBranch("muon_isHighPt")->SetTitle("muon tagged high pt");
   mtree->Branch("muon_dxy",&muon_dxy);
   mtree->GetBranch("muon_dxy")->SetTitle("muon transverse plane impact parameter (mm)");
   mtree->Branch("muon_dz",&muon_dz);
   mtree->GetBranch("muon_dz")->SetTitle("muon longitudinal impact parameter (mm)");
   mtree->Branch("muon_dxyError",&muon_dxyError);
   mtree->GetBranch("muon_dxyError")->SetTitle("muon transverse impact parameter uncertainty (mm)");
   mtree->Branch("muon_dzError",&muon_dzError);
   mtree->GetBranch("muon_dzError")->SetTitle("muon longitudinal impact parameter uncertainty (mm)");
   mtree->Branch("muon_pfreliso03all",&muon_pfreliso03all);
   mtree->GetBranch("muon_pfreliso03all")->SetTitle("muon particle flow relative isolation cone 03");
   mtree->Branch("muon_pfreliso04all",&muon_pfreliso04all);
   mtree->GetBranch("muon_pfreliso04all")->SetTitle("muon particle flow relative isolation cone 04");
   mtree->Branch("muon_pfreliso04DBCorr", &muon_pfreliso04DBCorr);
   mtree->GetBranch("muon_pfreliso04DBCorr")->SetTitle("muon particle flow relative isolation cone 04 DB Corrected");
   mtree->Branch("muon_TkIso03",&muon_TkIso03);
   mtree->GetBranch("muon_TkIso03")->SetTitle("muon tracker based isolation with cone size 03");
   mtree->Branch("muon_jetidx",&muon_jetidx);
   mtree->GetBranch("muon_jetidx")->SetTitle("index of the associated jet");
   mtree->Branch("muon_genpartidx",&muon_genpartidx);
   mtree->GetBranch("muon_genpartidx")->SetTitle("index into genParticle list for MC matching to status==1 muons");
   mtree->Branch("muon_ip3d", &muon_ip3d);
   mtree->GetBranch("muon_ip3d")->SetTitle("muon impact parameter in 3d");
   mtree->Branch("muon_sip3d", &muon_sip3d);
   mtree->GetBranch("muon_sip3d")->SetTitle("muon significance on impact parameter in 3d");
}

// Destructor.
Muon_Analyzer::~Muon_Analyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}

// Member functions.
// Method called for each event.
void
Muon_Analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

   Handle<pat::MuonCollection> muons;
   iEvent.getByToken(muonToken_, muons);
   // The result of getToken() method is an object  called "muons",
   //which is a collection of all the muon objects.
   // Collection classes are generally constructed as vectors.

   Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);
   // Skip the event if no PV is found.
   if(vertices->empty()) return;
   const reco::Vertex& PV = vertices->front();

   num_muon = 0;
   muon_e.clear();
   muon_pt.clear();
   muon_px.clear();
   muon_py.clear();
   muon_pz.clear();
   muon_eta.clear();
   muon_phi.clear();
   muon_ch.clear();
   muon_isLoose.clear();
   muon_isMedium.clear();
   muon_isTight.clear();
   muon_isSoft.clear();
   muon_isHighPt.clear();
   muon_dxy.clear();
   muon_dz.clear();
   muon_dxyError.clear();
   muon_dzError.clear();
   muon_pfreliso03all.clear();
   muon_pfreliso04all.clear();
   muon_pfreliso04DBCorr.clear();
   muon_TkIso03.clear();
   muon_jetidx.clear();
   muon_genpartidx.clear();
   muon_ip3d.clear();
   muon_sip3d.clear();

   // Muons loop.
   for(const pat::Muon& mu : *muons){
     muon_e.push_back(mu.energy());
     muon_pt.push_back(mu.pt());
     muon_px.push_back(mu.px());
     muon_py.push_back(mu.py());
     muon_pz.push_back(mu.pz());
     muon_eta.push_back(mu.eta());
     muon_phi.push_back(mu.phi());
     muon_ch.push_back(mu.charge());
     muon_isLoose.push_back(mu.isLooseMuon());
     muon_isMedium.push_back(mu.isMediumMuon());
     muon_isTight.push_back(mu.isTightMuon(PV));
     muon_isSoft.push_back(mu.isSoftMuon(PV));
     muon_isHighPt.push_back(mu.isHighPtMuon(PV));
     muon_dxy.push_back(mu.muonBestTrack()->dxy(PV.position()));
     muon_dz.push_back(mu.muonBestTrack()->dz(PV.position()));
     muon_dxyError.push_back(mu.muonBestTrack()->d0Error());
     muon_dzError.push_back(mu.muonBestTrack()->dzError());
     auto iso03 = mu.pfIsolationR03();
     muon_pfreliso03all.push_back((iso03.sumChargedHadronPt + iso03.sumNeutralHadronEt + iso03.sumPhotonEt)/mu.pt());
     auto iso04 = mu.pfIsolationR04();
     muon_pfreliso04all.push_back((iso04.sumChargedHadronPt + iso04.sumNeutralHadronEt + iso04.sumPhotonEt)/mu.pt());
     muon_pfreliso04DBCorr.push_back((iso04.sumChargedHadronPt + max(0., iso04.sumNeutralHadronEt + iso04.sumPhotonEt - 0.5*iso04.sumPUPt))/mu.pt());
     auto TkIso03 = mu.isolationR03();
     muon_TkIso03.push_back(TkIso03.sumPt/mu.pt());
     muon_genpartidx.push_back(-1);
     muon_jetidx.push_back(-1);

     // Get impact parameter in 3D.
     // https://github.com/cms-sw/cmssw/blob/CMSSW_7_6_X/PhysicsTools/PatAlgos/plugins/PATElectronProducer.cc
     // This is needed by the IPTools methods from the tracking group.
     edm::ESHandle<TransientTrackBuilder> trackBuilder;
     iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilder);
     reco::TransientTrack tt = trackBuilder->build(mu.muonBestTrack());
     std::pair<bool,Measurement1D> ip3dv = IPTools::absoluteImpactParameter3D(tt, PV);
     muon_ip3d.push_back(ip3dv.second.value());
     muon_sip3d.push_back(ip3dv.second.significance());

     num_muon++;
   } // End of Muons loop.
   mtree->Fill();
   return;
}


// Method called once each job just before starting event loop.
void
Muon_Analyzer::beginJob()
{
  std::cout << "Test Bool: " << testBool_ << std::endl;
}

// Method called once each job just after ending the event loop.
void
Muon_Analyzer::endJob()
{
}

// Method fills 'descriptions' with the allowed parameters for the module.
void
Muon_Analyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("muons", edm::InputTag("slimmedMuons"));
  desc.add<edm::InputTag>("vertices", edm::InputTag("offlineSlimmedPrimaryVertices"));
  desc.add<bool>("testBool", true);
  descriptions.add("muonAnalyzer", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Muon_Analyzer);
