// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// Additional includes.
#include "FWCore/Utilities/interface/InputTag.h"

// Classes to output ROOT files.
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

// Classes to extract Muon information.
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

// Classes to save data.
#include "TTree.h"
#include "TFile.h"
#include <vector>
////////////////////////////////////////////////////////////////////////////////
// class declaration
class Muon_Analyzer : public edm::one::EDAnalyzer<>  {
   public:
      explicit Muon_Analyzer(const edm::ParameterSet&);
      ~Muon_Analyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // Declare your tokens.
      edm::EDGetTokenT<pat::MuonCollection> muonToken_;
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;

      // Other member data.
      TTree* mtree;
      int num_muon; //number of Muons in the event
      // vectors to store electrons' properties.
      // Check the properties here: https://cmssdt.cern.ch/lxr/source/DataFormats/Candidate/interface/Particle.h?v=CMSSW_7_6_7
      std::vector<float> muon_e;
      std::vector<float> muon_pt;
      std::vector<float> muon_px;
      std::vector<float> muon_py;
      std::vector<float> muon_pz;
      std::vector<float> muon_eta;
      std::vector<float> muon_phi;
      std::vector<float> muon_ch;
      std::vector<int> muon_isLoose;
      std::vector<int> muon_isMedium;
      std::vector<int> muon_isTight;
      std::vector<int> muon_isSoft;
      std::vector<int> muon_isHighPt;
      std::vector<float> muon_dxy;
      std::vector<float> muon_dz;
      std::vector<float> muon_dxyError;
      std::vector<float> muon_dzError;
      std::vector<float> muon_pfreliso03all;
      std::vector<float> muon_pfreliso04all;
      std::vector<float> muon_TkIso03;
      std::vector<float> muon_genpartidx;
      std::vector<float> muon_jetidx;

      bool testBool_;

};
