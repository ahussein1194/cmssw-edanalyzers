// -*- C++ -*-
//
// Package:    my_packages/muon_e
// Class:      muon_e
// 
/**\class muon_e muon_e.cc my_packages/muon_e/plugins/muon_e.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Wed, 07 Dec 2022 14:09:57 GMT
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

// classes to extract muon info.
#include "DataFormats/PatCandidates/interface/Muon.h"

// standard C++ vector library.
#include <vector>

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class muon_e : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit muon_e(const edm::ParameterSet&);
      ~muon_e();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      std::vector<float> muonE; //energy values for muons in the event.
      edm::EDGetTokenT<pat::MuonCollection> muonToken_;
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
muon_e::muon_e(const edm::ParameterSet& iConfig):
	muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))) //here we are reading the 'Collection' variable from configuration.

{
   //now do what ever initialization is needed
   usesResource("TFileService");
   
   // Hardcoding the tag.
   //muonToken_ = consumes<pat::MuonCollection>(edm::InputTag("slimmedMuons")); // "slimmedMuons" used from (edmDumpEventContent) then checking the 'Label'
   
   // We are making the selection of this tag configurable.
   //muonToken_ = consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"));
}


muon_e::~muon_e()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
muon_e::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   // Clean the container.
   muonE.clear();

   // Define the handler, a token and get the information by toke.
   Handle<pat::MuonCollection> mymuons;
   iEvent.getByToken(muonToken_, mymuons);

   // If the collection is valid, loop over muons in the event.
   if(mymuons.isValid()) {
      for(const pat::Muon& itmuon : *mymuons) {
         // "pat::Muon": Look at class definitions.
         muonE.push_back(itmuon.energy());
         // "pat::Muon::energy": Look at class definitions
      }
   }

   // Print the vector of muons for each event.
   for(unsigned int i=0; i < muonE.size(); i++) {
   std::cout<< "Muon # " << i << " with E = " << muonE.at(i) << " GeV." << std::endl; 
   }



#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
muon_e::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
muon_e::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
muon_e::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(muon_e);
