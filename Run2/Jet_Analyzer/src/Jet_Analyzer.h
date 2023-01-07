// system include files
#include <memory>
#include <TMath.h>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Ref.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "math.h"
// classes to extract jet information.
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetResolutionObject.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
// Classes to output ROOT files.
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
// Classes to save data.
#include "TTree.h"
#include "TFile.h"
#include <vector>

// Other classes.
#include "TRandom3.h"
////////////////////////////////////////////////////////////////////////////////

// class declaration
class Jet_Analyzer : public edm::one::EDAnalyzer<>  {
   public:
      explicit Jet_Analyzer(const edm::ParameterSet&);
      ~Jet_Analyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob();
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob();

      virtual double getBtagEfficiency(double pt);
      virtual double getCtagEfficiency(double pt);
      virtual double getLFtagEfficiency(double pt);
      virtual double getBorCtagSF(double pt, double eta);
      virtual double getLFtagSF(double pt, double eta);
      virtual double uncertaintyForBTagSF( double pt, double eta);
      virtual double uncertaintyForLFTagSF( double pt, double eta);
      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override; //Original
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // Declare the input tag for PFJetCollection.
      edm::EDGetTokenT<pat::JetCollection> jetToken_;
      edm::EDGetTokenT<double> rhoToken_;

      // Other member data.
      // Jet variables.
      std::string jetJECUncName_;
      std::string jetResName_;
      std::string sfName_;
      boost::shared_ptr<JetCorrectionUncertainty> jecUnc_;
      bool isData;

      JME::JetResolution resolution;
      JME::JetResolutionScaleFactor resolution_sf;

      int num_jets; //number of jets in the event
      TTree *mtree; // tree to store jet properties
      std::vector<float> jet_e;
      std::vector<float> jet_pt;
      std::vector<float> jet_eta;
      std::vector<float> jet_phi;
      std::vector<float> jet_ch;
      std::vector<float> jet_mass;
      std::vector<double> jet_btag;
      std::vector<int>   jet_hflav;
      std::vector<float> jet_corrpt;
      std::vector<float> jet_corrptUp;
      std::vector<float> jet_corrptDown;
      std::vector<float> jet_corrptSmearUp;
      std::vector<float> jet_corrptSmearDown;
      std::vector<float> jet_corrmass;
      std::vector<float> jet_corre;
      std::vector<float> jet_corrpx;
      std::vector<float> jet_corrpy;
      std::vector<float> jet_corrpz;
      float btagWeight;
      float btagWeightUp;
      float btagWeightDn;
};
