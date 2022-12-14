// -*- C++ -*-
//
// Package:    my_EDAnalyzers/Electron_Analyzer
// Class:      Electron_Analyzer
//
/**\class Electron_Analyzer Electron_Analyzer.cc my_EDAnalyzers/Electron_Analyzer/plugins/Electron_Analyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:
//         Created:  Fri, 09 Dec 2022 12:09:43 GMT
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

// additional includes.
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/ValueMap.h"

// classes to output root files.
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

// classes to extract electron information.
// Hints (miniAOD - Run2):
// -https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2015#High_level_physics_objects
// -https://cms-opendata-guide.web.cern.ch/analysis/selection/objects/objects/#dataformats
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

// classes to save data.
#include "TTree.h"
#include "TFile.h"
#include <vector>
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class Electron_Analyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
//class Electron_Analyzer : public edm::one::EDAnalyzer<>  {

   public:
      explicit Electron_Analyzer(const edm::ParameterSet&);
      ~Electron_Analyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      /// Registering for data access.
      //(1): Declare the token for ElectronCollection and VertexCollection.
      edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      // or
      //edm::EDGetTokenT<std::vector<pat::Electron>> electronToken_, electronToken2_;
      //edm::EDGetTokenT<std::vector<reco::Vertex>> vtxToken_;

      // How to know <pat::ElectronCollection>?
      // Go to the class you included above for electrons:
      // "DataFormats/PatCandidates/interface/Electron.h" on LXR or Doxygen.
      // Find the type defined for "vector<pat::Electron>" type for "slimmedElectrons"
      // Same for ["DataFormats/VertexReco/interface/VertexFwd.h", "DataFormats/VertexReco/interface/Vertex.h"]

      // ----------member data ---------------------------

      TTree* mtree;
      int num_electron; //number of electrons in the event
      // vectors to store electrons' properties.
      // Check the properties here: https://cmssdt.cern.ch/lxr/source/DataFormats/Candidate/interface/Particle.h?v=CMSSW_7_6_7
      std::vector<float> electron_e;
      std::vector<float> electron_pt;
      std::vector<float> electron_px;
      std::vector<float> electron_py;
      std::vector<float> electron_pz;
      std::vector<float> electron_eta;
      std::vector<float> electron_phi;
      std::vector<float> electron_ch;
      std::vector<float> electron_iso;
      std::vector<bool> electron_veto;
      std::vector<bool> electron_isLoose;
      std::vector<bool> electron_isMedium;
      std::vector<bool> electron_isTight;
      std::vector<float> electron_dxy;
      std::vector<float> electron_dz;
      std::vector<float> electron_dxyError;
      std::vector<float> electron_dzError;
      std::vector<int> electron_ismvaLoose;
      std::vector<int> electron_ismvaTight;

      bool testBool_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Electron_Analyzer::Electron_Analyzer(const edm::ParameterSet& iConfig):
 /// Registering for data access.
 //(2): Pass the Input Tag to the token defined above.
 electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
 vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
 testBool_(iConfig.getParameter<bool>("testBool"))
 // or
 //electronToken_(consumes<std::vector<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
 //vtxToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices")))
{
   //now do what ever initialization is needed
   //usesResource("TFileService");
   edm::Service<TFileService> fs;
   mtree = fs->make<TTree>("Events", "Events");
   //or (if you did not define it before in member data)
   //TTree* mtree = fs->make<TTree>("Events", "Events");

  mtree->Branch("numberelectron", &num_electron);
  mtree->GetBranch("numberelectron")->SetTitle("number of electrons");
  mtree->Branch("electron_e",&electron_e);
  mtree->GetBranch("electron_e")->SetTitle("electron energy");
  mtree->Branch("electron_pt",&electron_pt);
  mtree->GetBranch("electron_pt")->SetTitle("electron transverse momentum");
  mtree->Branch("electron_px",&electron_px);
  mtree->GetBranch("electron_px")->SetTitle("electron momentum x-component");
  mtree->Branch("electron_py",&electron_py);
  mtree->GetBranch("electron_py")->SetTitle("electron momentum y-component");
  mtree->Branch("electron_pz",&electron_pz);
  mtree->GetBranch("electron_pz")->SetTitle("electron momentum z-component");
  mtree->Branch("electron_eta",&electron_eta);
  mtree->GetBranch("electron_eta")->SetTitle("electron pseudorapidity");
  mtree->Branch("electron_phi",&electron_phi);
  mtree->GetBranch("electron_phi")->SetTitle("electron polar angle");
  mtree->Branch("electron_ch",&electron_ch);
  mtree->GetBranch("electron_ch")->SetTitle("electron charge");
  mtree->Branch("electron_iso",&electron_iso);
  mtree->GetBranch("electron_iso")->SetTitle("electron isolation");
  mtree->Branch("electron_veto",&electron_veto);//
  mtree->GetBranch("electron_veto")->SetTitle("electron veto");//
  mtree->Branch("electron_isLoose",&electron_isLoose);
  mtree->GetBranch("electron_isLoose")->SetTitle("electron tagged loose");
  mtree->Branch("electron_isMedium",&electron_isMedium);
  mtree->GetBranch("electron_isMedium")->SetTitle("electron tagged medium");
  mtree->Branch("electron_isTight",&electron_isTight);
  mtree->GetBranch("electron_isTight")->SetTitle("electron tagged tight");
  mtree->Branch("electron_dxy",&electron_dxy);
  mtree->GetBranch("electron_dxy")->SetTitle("electron transverse plane impact parameter (mm)");
  mtree->Branch("electron_dz",&electron_dz);
  mtree->GetBranch("electron_dz")->SetTitle("electron longitudinal impact parameter (mm)");
  mtree->Branch("electron_dxyError",&electron_dxyError);
  mtree->GetBranch("electron_dxyError")->SetTitle("electron transverse impact parameter uncertainty (mm)");
  mtree->Branch("electron_dzError",&electron_dzError);
  mtree->GetBranch("electron_dzError")->SetTitle("electron longitudinal impact parameter uncertainty (mm)");
  mtree->Branch("electron_ismvaLoose",&electron_ismvaLoose);
  mtree->GetBranch("electron_ismvaLoose")->SetTitle("electron mva Loose");
  mtree->Branch("electron_ismvaTight",&electron_ismvaTight);
  mtree->GetBranch("electron_ismvaTight")->SetTitle("electron mva Tight");
}

// Destructor.
Electron_Analyzer::~Electron_Analyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Electron_Analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   /// Registering for data access.
   //(3): Now, use the token.
   Handle<pat::ElectronCollection> electrons;
   iEvent.getByToken(electronToken_, electrons);

   Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);
   math::XYZPoint pv(vertices->begin()->position());

   num_electron = 0;
   electron_e.clear();
   electron_pt.clear();
   electron_px.clear();
   electron_py.clear();
   electron_pz.clear();
   electron_eta.clear();
   electron_phi.clear();
   electron_ch.clear();
   electron_iso.clear();
   electron_veto.clear();//
   electron_isLoose.clear();
   electron_isMedium.clear();
   electron_isTight.clear();
   electron_dxy.clear();
   electron_dz.clear();
   electron_dxyError.clear();
   electron_dzError.clear();
   electron_ismvaLoose.clear();
   electron_ismvaTight.clear();

   // "pat::ElectronCollection" is a vector of "pat::Electron" objects,
   // Look at: https://cmssdt.cern.ch/lxr/source/DataFormats/PatCandidates/interface/Electron.h?v=CMSSW_7_6_7
   // "electrons" is a pointer to a vector of "pat::Electron" objects?
   for(const pat::Electron& el : *electrons){
     electron_e.push_back(el.energy());
     electron_pt.push_back(el.pt());
     electron_px.push_back(el.px());
      electron_py.push_back(el.py());
      electron_pz.push_back(el.pz());
      electron_eta.push_back(el.eta());
      electron_phi.push_back(el.phi());
      electron_ch.push_back(el.charge());
      electron_iso.push_back(el.ecalPFClusterIso());
      electron_veto.push_back(el.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-veto"));//
      electron_isLoose.push_back(el.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-loose"));
      electron_isMedium.push_back(el.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-medium"));
      electron_isTight.push_back(el.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-tight"));
      electron_dxy.push_back(el.gsfTrack()->dxy(pv));
      electron_dz.push_back(el.gsfTrack()->dz(pv));
      electron_dxyError.push_back(el.gsfTrack()->d0Error());
      electron_dzError.push_back(el.gsfTrack()->dzError());
      electron_ismvaLoose.push_back(el.electronID("mvaEleID-Spring15-25ns-nonTrig-V1-wp90"));
      electron_ismvaTight.push_back(el.electronID("mvaEleID-Spring15-25ns-nonTrig-V1-wp80"));
      num_electron++;
   }

   mtree->Fill();
   //return;

/*#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif*/
}


// ------------ method called once each job just before starting event loop  ------------
void
Electron_Analyzer::beginJob()
{
  std::cout<< "testBool: " << testBool_ << std::endl;
}

// ------------ method called once each job just after ending the event loop  ------------
void
Electron_Analyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Electron_Analyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Electron_Analyzer);
