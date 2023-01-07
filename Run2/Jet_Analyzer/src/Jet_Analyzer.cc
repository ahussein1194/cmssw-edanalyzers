// -*- C++ -*-
//
// Package:    my_EDAnalyzers/Jet_Analyzer
// Class:      Jet_Analyzer
//
/**\class Jet_Analyzer Jet_Analyzer.cc my_EDAnalyzers/Jet_Analyzer/plugins/Jet_Analyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author: Ahmed Hussein
//         Created:  Wed, 21 Dec 2022 21:52:54 GMT
////////////////////////////////////////////////////////////////////////////////
#include "Jet_Analyzer.h"

// Constructor.
Jet_Analyzer::Jet_Analyzer(const edm::ParameterSet& iConfig):
  jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
  jetJECUncName_(iConfig.getParameter<edm::FileInPath>("jetJECUncName").fullPath()),
  jetResName_(iConfig.getParameter<edm::FileInPath>("jerResName").fullPath()),
  sfName_(iConfig.getParameter<edm::FileInPath>("jerSFName").fullPath()),
  isData(iConfig.getParameter<bool>("isData"))
{
   //now do what ever initialization is needed
   edm::Service<TFileService> fs;
   mtree = fs->make<TTree>("Events", "Events");

   jecUnc_ = boost::shared_ptr<JetCorrectionUncertainty>(new JetCorrectionUncertainty(jetJECUncName_));
   resolution = JME::JetResolution(jetResName_);
   resolution_sf = JME::JetResolutionScaleFactor(sfName_);

   edm::InputTag rhotag("fixedGridRhoFastjetAll");
   rhoToken_ = consumes<double>(rhotag);

   // Create the TTree TBranches.
   mtree->Branch("numberjet", &num_jets);
   mtree->GetBranch("numberjet")->SetTitle("Number of Jets");
   mtree->Branch("jet_e",&jet_e);
   mtree->GetBranch("jet_e")->SetTitle("Uncorrected Jet energy");
   mtree->Branch("jet_pt",&jet_pt);
   mtree->GetBranch("jet_pt")->SetTitle("Uncorrected Transverse Jet Momentum");
   mtree->Branch("jet_eta",&jet_eta);
   mtree->GetBranch("jet_eta")->SetTitle("Jet Eta");
   mtree->Branch("jet_phi",&jet_phi);
   mtree->GetBranch("jet_phi")->SetTitle("Jet Phi");
   mtree->Branch("jet_ch",&jet_ch);
   mtree->GetBranch("jet_ch")->SetTitle("Jet Charge");
   mtree->Branch("jet_mass",&jet_mass);
   mtree->GetBranch("jet_mass")->SetTitle("Jet Mass");
   mtree->Branch("jet_btag", &jet_btag);
   mtree->GetBranch("jet_btag")->SetTitle("Jet Btagging Discriminant (CSV)");
   mtree->Branch("jet_hflav",&jet_hflav);
   mtree->GetBranch("jet_hflav")->SetTitle("Jet hadron-based flavour (in MC)");
   mtree->Branch("jet_corrpt",&jet_corrpt);
   mtree->GetBranch("jet_corrpt")->SetTitle("Corrected Transverse Jet Momentum");
   mtree->Branch("jet_corrptUp",&jet_corrptUp);
   mtree->GetBranch("jet_corrptUp")->SetTitle("Corrected Transverse Jet Momentum (JEC Shifted Up)");
   mtree->Branch("jet_corrptDown",&jet_corrptDown);
   mtree->GetBranch("jet_corrptDown")->SetTitle("Corrected Transverse Jet Momentum (JEC Shifted Down)");
   mtree->Branch("jet_corrptSmearUp",&jet_corrptSmearUp);
   mtree->GetBranch("jet_corrptSmearUp")->SetTitle("Corrected Transverse Jet Momentum (JER Shifted Up)");
   mtree->Branch("jet_corrptSmearDown",&jet_corrptSmearDown);
   mtree->GetBranch("jet_corrptSmearDown")->SetTitle("Corrected Transverse Jet Momentum (JER Shifted Down)");
   mtree->Branch("jet_corrmass",&jet_corrmass);
   mtree->GetBranch("jet_corrmass")->SetTitle("Corrected Jet Mass");
   mtree->Branch("jet_corre",&jet_corre);
   mtree->GetBranch("jet_corre")->SetTitle("Corrected Jet Energy");
   mtree->Branch("jet_corrpx",&jet_corrpx);
   mtree->GetBranch("jet_corrpx")->SetTitle("Corrected X-Component of Jet Momentum");
   mtree->Branch("jet_corrpy",&jet_corrpy);
   mtree->GetBranch("jet_corrpy")->SetTitle("Corrected Y-Component of Jet Momentum");
   mtree->Branch("jet_corrpz",&jet_corrpz);
   mtree->GetBranch("jet_corrpz")->SetTitle("Corrected Z-Component of Jet Momentum");
   mtree->Branch("btag_Weight", &btagWeight);
   mtree->GetBranch("btag_Weight")->SetTitle("B-Tag Event Weight");
   mtree->Branch("btag_WeightUp", &btagWeightUp);
   mtree->GetBranch("btag_WeightUp")->SetTitle("B-Tag Up Event Weight");
   mtree->Branch("btag_WeightDn", &btagWeightDn);
   mtree->GetBranch("btag_WeightDn")->SetTitle("B-Tag Down Event Weight");
}

// Destructor.
Jet_Analyzer::~Jet_Analyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}

// Member functions
//
// Function that gets the efficiences for B tag jets.
//returns hard coded values from data analyzed from out root -l plotBEff.C
double
Jet_Analyzer::getBtagEfficiency(double pt){
  // Efficiencies from TTJets madgraph sample miniAODv2
  // ******* REPLACE ME BY USING THE BTAGGING MODULE!!! ***********
  if(pt < 20) return 0.228679;
  else if(pt < 40) return 0.472444;
  else if(pt < 60) return 0.572408;
  else if(pt < 80) return 0.608838;
  else if(pt < 100) return 0.629398;
  else if(pt < 120) return 0.634980;
  else if(pt < 140) return 0.645666;
  else if(pt < 2150) return 0.706965 - 0.000328*pt;
  else return 0;
}

// Function that gets the efficiences for C tag jets.
//returns hard coded values from data analyzed from out root -l plotBEff.C
double
Jet_Analyzer::getCtagEfficiency(double pt){
  // Hand-waving assumption that C efficiency is 20% of the B efficiency
  // ******* REPLACE ME BY USING THE BTAGGING MODULE!!! ***********

  return 0.2*getBtagEfficiency(pt);
}

// Function that gets the efficiences for LF tag jets.
//returns hard coded values from data analyzed from out root -l plotBEff.C
double
Jet_Analyzer::getLFtagEfficiency(double pt){
  // Mistag rates from TTJets madgraph sample miniAODv2.
  // ******* REPLACE ME BY USING THE BTAGGING MODULE!!! ***********
  if(pt < 20) return 0.003402;
  else if(pt < 40) return 0.008067;
  else if(pt < 60) return 0.006980;
  else if(pt < 80) return 0.006316;
  else if(pt < 100) return 0.006713;
  else if(pt < 115) return 0.006598;
  else if(pt < 400) return -0.00242504 + 9.15452e-05*pt - 9.63553e-08*pt*pt; // functional form is overkill, but oh well.
  else return 0.0134038 + 1.24358e-05*pt; // same here, you can just use a binned value
}

// Function that gets the SF for B or C tag jets since their equations from the CSV.csv files are the same.
//returns the SF of the B or C tag jet depending on the pt of the jet
double
Jet_Analyzer::getBorCtagSF(double pt, double eta){
  if(pt > 670) pt = 670;
  if(fabs(eta) > 2.4 or pt < 30.) return 1.0;

  return 0.934588*((1.+(0.00678184*pt))/(1.+(0.00627144*pt)));
}

// Function that gets the SF for lf tag jets since there equations from the CSV.csv files are the same.
//returns the SF of the lf jet depending on the eta of the jet
double
Jet_Analyzer::getLFtagSF(double pt, double eta){
  if (pt > 1000.) pt = 1000;
  if(fabs(eta) > 2.4 or pt < 20.) return 1.0;

  if(fabs(eta) < 0.8) return ((0.994351+(0.000250077*pt))+(9.24801e-07*(pt*pt)))+(-8.73293e-10*(pt*(pt*pt)));
  else if(fabs(eta) < 1.6) return ((1.00939+(0.000461283*pt))+(-6.30306e-07*(pt*pt)))+(3.53075e-10*(pt*(pt*pt)));
  else return ((0.955798+(0.00146058*pt))+(-3.76689e-06*(pt*pt)))+(2.39196e-09*(pt*(pt*pt)));
}

// Function that gets the uncertainty for B tag jets.
//returns the uncertainty from the CSVV2.csv file depending on the pt of the jet.
//This will also be used to find the uncertainty of C tag jets but will be multiplied by 2 since the uncertainty of C tag jets are 2 times that of B tag jets.
double
Jet_Analyzer::uncertaintyForBTagSF(double pt, double eta){
  if(fabs(eta) > 2.4 or pt<20.) return 0;

  if(pt < 30) return 0.018076473847031593;
  else if(pt < 50) return 0.024799736216664314;
  else if(pt < 70) return 0.024073265492916107;
  else if(pt < 100) return 0.020040607079863548;
  else if(pt < 140) return 0.016540588811039925;
  else if(pt < 200) return 0.025977084413170815;
  else return 0.027120551094412804;
}

// Function that gets the uncertainty for LF tag jets.
//returns the uncertainty from the CSVV2.csv file depending on the eta of the jet.
double
Jet_Analyzer::uncertaintyForLFTagSF(double pt, double eta){
  if(pt > 1000.)  pt = 1000;
  if(fabs(eta) > 2.4 or pt<20.)  return 0;
  if(fabs(eta) < 0.8)
    return 0.5*((((1.03928+(0.000857422*pt))+(-4.02756e-07*(pt*pt)))+(-8.45836e-11*(pt*(pt*pt)))) - (((0.949401+(-0.000356232*pt))+(2.24887e-06*(pt*pt)))+(-1.66011e-09*(pt*(pt*pt)))));
  else if (fabs(eta) < 1.6)
    return 0.5*((((1.05392+(0.000944135*pt))+(-1.73386e-06*(pt*pt)))+(1.04242e-09*(pt*(pt*pt)))) - (((0.964857+(-2.19898e-05*pt))+(4.74117e-07*(pt*pt)))+(-3.36548e-10*(pt*(pt*pt)))));
  else
    return 0.5*((((1.00151+(0.00175547*pt))+(-4.50251e-06*(pt*pt)))+(2.91473e-09*(pt*(pt*pt)))) - (((0.910086+(0.00116371*pt))+(-3.02747e-06*(pt*pt)))+(1.86906e-09*(pt*(pt*pt)))));
}

// ------------ method called for each event  ------------
void
Jet_Analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

   Handle<pat::JetCollection> jets;
   iEvent.getByToken(jetToken_, jets);
   Handle<double> rhoHandle;
   iEvent.getByToken(rhoToken_, rhoHandle);

   // Clear the varibles before beginning the event loop.
   num_jets = 0;
   jet_e.clear();
   jet_pt.clear();
   jet_eta.clear();
   jet_phi.clear();
   jet_ch.clear();
   jet_mass.clear();
   jet_btag.clear();
   jet_hflav.clear();
   jet_corrpt.clear();
   jet_corrptUp.clear();
   jet_corrptDown.clear();
   jet_corrptSmearUp.clear();
   jet_corrptSmearDown.clear();
   jet_corrmass.clear();
   jet_corre.clear();
   jet_corrpx.clear();
   jet_corrpy.clear();
   jet_corrpz.clear();

   int hadronFlavour;
   double SF, SFu, SFd, eff, corrpt;
   double corrUp, corrDown;
   float ptscale, ptscale_down, ptscale_up;
   double MC = 1;
   btagWeight = 1;
   btagWeightUp = 1;
   btagWeightDn = 1;
   int min_pt = 20;

   if(jets.isValid()){
     // Jets loop.
     for(const pat::Jet& jet : *jets){
       pat::Jet uncorrJet = jet.correctedJet(0);

       corrpt = jet.pt();
       corrUp = 1.0;
       corrDown = 1.0;

       if(fabs(jet.eta()) < 5) jecUnc_->setJetEta(jet.eta());
       else jecUnc_->setJetEta(4.99);
       jecUnc_->setJetPt(corrpt);
       corrDown = (1 - fabs(jecUnc_->getUncertainty(-1)));

       ptscale = 1;
       ptscale_down = 1;
       ptscale_up = 1;

       if(!isData){
         JME::JetParameters JERparameters = {{JME::Binning::JetPt, corrpt}, {JME::Binning::JetEta, jet.eta()}, {JME::Binning::Rho, *(rhoHandle.product())}};
	       float res = resolution.getResolution(JERparameters);
	       float sf = resolution_sf.getScaleFactor(JERparameters);
	       float sf_up = resolution_sf.getScaleFactor(JERparameters, Variation::UP);
	       float sf_down = resolution_sf.getScaleFactor(JERparameters, Variation::DOWN);

	       const reco::GenJet *genJet = jet.genJet();
	       bool smeared = false;
         if(genJet){
           double deltaPt = fabs(genJet->pt() - corrpt);
	         double deltaR = reco::deltaR(genJet->p4(),jet.p4());
           if((deltaR < 0.2) && deltaPt <= 3*corrpt*res){
             ptscale = max(0.0, 1 + (sf - 1.0)*(corrpt - genJet->pt())/corrpt);
	           ptscale_down = max(0.0, 1 + (sf_down - 1.0)*(corrpt - genJet->pt())/corrpt);
	           ptscale_up = max(0.0, 1 + (sf_up - 1.0)*(corrpt - genJet->pt())/corrpt);
	           smeared = true;
           }
         }
         if(!smeared){
           TRandom3 JERrand;

	         JERrand.SetSeed(abs(static_cast<int>(jet.phi()*1e4)));
	         ptscale = max(0.0, 1.0 + JERrand.Gaus(0, res)*sqrt(max(0.0, sf*sf - 1.0)));

	         JERrand.SetSeed(abs(static_cast<int>(jet.phi()*1e4)));
	         ptscale_down = max(0.0, 1.0 + JERrand.Gaus(0, res)*sqrt(max(0.0, sf_down*sf_down - 1.0)));

	         JERrand.SetSeed(abs(static_cast<int>(jet.phi()*1e4)));
	         ptscale_up = max(0.0, 1.0 + JERrand.Gaus(0, res)*sqrt(max(0.0, sf_up*sf_up - 1.0)));
         }
       }
       if(ptscale*corrpt <= min_pt) continue;

       pat::Jet smearedjet = jet;
       smearedjet.scaleEnergy(ptscale);

       jet_e.push_back(uncorrJet.energy());
       jet_pt.push_back(uncorrJet.pt());
       jet_eta.push_back(uncorrJet.eta());
       jet_phi.push_back(uncorrJet.phi());
       jet_ch.push_back(uncorrJet.charge());
       jet_mass.push_back(uncorrJet.mass());
       jet_btag.push_back(smearedjet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
       jet_hflav.push_back(smearedjet.hadronFlavour());
       jet_corrpt.push_back(smearedjet.pt());
       jet_corrptUp.push_back(corrUp*smearedjet.pt());
       jet_corrptDown.push_back(corrDown*smearedjet.pt());
       jet_corrptSmearUp.push_back(ptscale_up*smearedjet.pt()/ptscale);
       jet_corrptSmearDown.push_back(ptscale_down*smearedjet.pt()/ptscale);
       jet_corrmass.push_back(smearedjet.mass());
       jet_corre.push_back(smearedjet.energy());
       jet_corrpx.push_back(smearedjet.px());
       jet_corrpy.push_back(smearedjet.py());
       jet_corrpz.push_back(smearedjet.pz());

       if(!isData){
         SF = 1;
         SFu = 1;
         SFd = 1;
         eff = 1;
         hadronFlavour = smearedjet.hadronFlavour();
         corrpt = jet_corrpt.at(num_jets);

         if(jet_btag.at(num_jets) > 0.800){ // MEDIUM working point
           if(abs(hadronFlavour) == 5){
             eff = getBtagEfficiency(corrpt);
	           SF = getBorCtagSF(corrpt, jet_eta.at(num_jets));
	           SFu = SF + uncertaintyForBTagSF(corrpt, jet_eta.at(num_jets));
	           SFd = SF - uncertaintyForBTagSF(corrpt, jet_eta.at(num_jets));
           }
           else if(abs(hadronFlavour) == 4){
             eff = getCtagEfficiency(corrpt);
	           SF = getBorCtagSF(corrpt, jet_eta.at(num_jets));
	           SFu = SF + (2 * uncertaintyForBTagSF(corrpt, jet_eta.at(num_jets)));
	           SFd = SF - (2 * uncertaintyForBTagSF(corrpt, jet_eta.at(num_jets)));
           }
           else{
             eff = getLFtagEfficiency(corrpt);
	           SF = getLFtagSF(corrpt, jet_eta.at(num_jets));
	           SFu = SF + ( uncertaintyForLFTagSF(corrpt, jet_eta.at(num_jets)));
	           SFd = SF - ( uncertaintyForLFTagSF(corrpt, jet_eta.at(num_jets)));
           }
           MC *= eff;
           btagWeight *= SF * eff;
           btagWeightUp *= SFu * eff;
	         btagWeightDn *= SFd * eff;
         }
         else {
           if(abs(hadronFlavour) == 5){
             eff = getBtagEfficiency(corrpt);
	           SF = getBorCtagSF(corrpt, jet_eta.at(num_jets));
	           SFu = SF + uncertaintyForBTagSF(corrpt, jet_eta.at(num_jets));
	           SFd = SF - uncertaintyForBTagSF(corrpt, jet_eta.at(num_jets));
           }
           else if(abs(hadronFlavour) == 4){
             eff = getCtagEfficiency(corrpt);
	           SF = getBorCtagSF(corrpt, jet_eta.at(num_jets));
	           SFu = SF + (2 * uncertaintyForBTagSF(corrpt, jet_eta.at(num_jets)));
	           SFd = SF - (2 * uncertaintyForBTagSF(corrpt, jet_eta.at(num_jets)));
           }
           else{
             eff = getLFtagEfficiency(corrpt);
	           SF = getLFtagSF(corrpt, jet_eta.at(num_jets));
	           SFu = SF + ( uncertaintyForLFTagSF(corrpt, jet_eta.at(num_jets)));
	           SFd = SF - ( uncertaintyForLFTagSF(corrpt, jet_eta.at(num_jets)));
           }
           MC *= (1 - eff);
           btagWeight *= (1 - ( SF * eff));
	         btagWeightUp *= (1 - (SFu * eff));
	         btagWeightDn *= (1 -  (SFd * eff));
         }
       }
       ++num_jets;
     }
     btagWeight = (btagWeight/MC);
     btagWeightUp = (btagWeightUp/MC);
     btagWeightDn = (btagWeightDn/MC);
   } // End of Jets loop.
   mtree->Fill();
   return;
}

// ------------ method called once each job just before starting event loop  ------------
void
Jet_Analyzer::beginJob()
{}

// ------------ method called once each job just after ending the event loop  ------------
void
Jet_Analyzer::endJob()
{}

// ------------ method called when starting to processes a run  ------------
void
Jet_Analyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a run  ------------
void
Jet_Analyzer::endRun(edm::Run const&, edm::EventSetup const&)
{}

// ------------ method called when starting to processes a luminosity block  ------------
void
Jet_Analyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a luminosity block  ------------
void
Jet_Analyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Jet_Analyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("jets", edm::InputTag("slimmedJetsNewJEC"));
  desc.add<bool>("isData");
  desc.add<edm::FileInPath>("jetJECUncName", edm::FileInPath("my_EDAnalyzers/Jet_Analyzer/JEC/Fall15_25nsV2_MC_Uncertainty_AK4PFchs.txt"));
  desc.add<edm::FileInPath>("jerResName", edm::FileInPath("my_EDAnalyzers/Jet_Analyzer/JEC/Fall15_25nsV2_MC_PtResolution_AK4PFchs.txt"));
  desc.add<edm::FileInPath>("jerSFName", edm::FileInPath("my_EDAnalyzers/Jet_Analyzer/JEC/Fall15_25nsV2_MC_SF_AK4PFchs.txt"));
  descriptions.add("jetAnalyzer", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Jet_Analyzer);
