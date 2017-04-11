#pragma once

#include "BTagCalibrationStandalone.h"
#include "BaseOperator.h"

#include "TH1.h"
#include "TH2.h"
#include "TRandom3.h"
#include "TString.h"

template <class EventClass> class TriggerStudiesOperator : public BaseOperator<EventClass> {

  public:

    std::vector<std::string> or_paths_; 
    std::vector<std::string> weights_;

    // get short disc name
    std::map < std::string, std::string> s_name_map;
    // get flavour from hadronFlavour
    std::map < int , BTagEntry::JetFlavor > flavour_map;
    // systematics to take into account per flavour
    std::map< BTagEntry::JetFlavor, std::vector<std::string>> syst_map;
    // Load BTagSF
    BTagCalibration btcalib;
    // BTagSF options
    std::map<std::string, BTagCalibrationReader> cr_map;
    std::string sf_mode = "iterativefit";

    const static int nsel   = 3;
    const static int npaths = 2;
    TString sselection[nsel  ]={"_All","_HLTbits","_HLTemulation"};
    TString spaths        [npaths+1]={"_QuadJet45","_DoubleJet90_Double30","_QuadJet_OR_DoubleJet"};
    TH1F* h_jet3_pt [nsel][npaths+1];
    TH1F* h_jet3_CSV[nsel][npaths+1];
    TH1F* h_HT      [nsel][npaths+1];
    TH2F h2_triggerEmulation {"h2_triggerEmulation","h2_triggerEmulation",npaths+1,0,npaths+1,npaths+1,0,npaths+1};
    TH2F h2_triggerBits      {"h2_triggerBits"     ,"h2_triggerBits"     ,npaths+1,0,npaths+1,npaths+1,0,npaths+1};

    TRandom3 r;

    TriggerStudiesOperator(std::vector<std::string> or_paths, const std::vector<std::string> & weights = {}) :
      or_paths_(or_paths),
      weights_(weights)
      {
//      disc_("pfCombinedMVAV2BJetTags"), // FIXME HACK: do not know how to pass it from the python.
      flavour_map = {{5, BTagEntry::FLAV_B},
                     {4, BTagEntry::FLAV_C},
                     {0, BTagEntry::FLAV_UDSG}};
    }
    virtual ~TriggerStudiesOperator() {}

    virtual void init(TDirectory * tdir) {
      h2_triggerBits     .SetDirectory(tdir);
      h2_triggerEmulation.SetDirectory(tdir);

      for (unsigned int isel=0; isel<nsel; ++isel) {
	  for (unsigned int ipath=0; ipath<npaths+1; ++ipath) {
	      if (isel==0 && ipath) continue;
	      h_jet3_pt [isel][ipath] = new TH1F("h_jet3_pt"  +(isel ? spaths[ipath] : TString(""))+sselection[isel], "jet3 pt"  +(isel ? spaths[ipath] : TString(""))+sselection[isel], 100, 0., 300. );
	      h_jet3_CSV[isel][ipath] = new TH1F("h_jet3_CSV" +(isel ? spaths[ipath] : TString(""))+sselection[isel], "jet3 CSV" +(isel ? spaths[ipath] : TString(""))+sselection[isel], 100, 0., 1. );
	      h_HT      [isel][ipath] = new TH1F("h_HT"       +(isel ? spaths[ipath] : TString(""))+sselection[isel], "HT"       +(isel ? spaths[ipath] : TString(""))+sselection[isel], 100, 0., 3000.);
	      h_jet3_pt [isel][ipath]->Sumw2();
	      h_jet3_CSV[isel][ipath]->Sumw2();
	      h_HT      [isel][ipath]->Sumw2();
	      h_jet3_pt [isel][ipath]->SetDirectory(tdir);
	      h_jet3_CSV[isel][ipath]->SetDirectory(tdir);
	      h_HT      [isel][ipath]->SetDirectory(tdir);
	  }
      }
    }

    virtual bool process( EventClass & ev ) {

      float w = 1.0;
      w*=ev.eventInfo_.eventWeight(weights_);
//      cout << endl;
      vector <float> vjet_pt;
      vector <float> vjet_eta;
      // order by pt
      std::string disc_("pt");
      order_jets_by_disc(ev.jets_, disc_);
      double sumpt=0;
      for (auto & jet : ev.jets_) {
	  auto jet_eta = jet.eta();
	  auto jet_pt = jet.pt();
	  if (abs(jet_eta)<2.4 && jet_pt>30) {
	      sumpt+=jet_pt;
	  }
	  vjet_pt .push_back(jet_pt );
	  vjet_eta.push_back(jet_eta);
      }
//      cout << "sumpt= " << sumpt << endl;

      // order by discriminator
      disc_="pfCombinedInclusiveSecondaryVertexV2BJetTags";
      vector <float> vjet_disc;
      order_jets_by_disc(ev.jets_, disc_);
      for (auto & jet : ev.jets_) {
	  vjet_disc.push_back(jet.disc(disc_));
      }

      vector <float> triggerEmulation(npaths,0);
      if (vjet_pt.size()<4) {
//	  cout << "*************NOT ENOUGH JETS***********" << endl;
	  return 0;
      }
      else {
	  triggerEmulation[0] = TurnOnQuad  (sumpt,vjet_pt.at(1),vjet_pt.at(3),vjet_disc.at(0),vjet_disc.at(1),vjet_disc.at(2),vjet_disc.at(3));
	  triggerEmulation[1] = TurnOnDouble(sumpt,vjet_pt.at(1),vjet_pt.at(3),vjet_disc.at(0),vjet_disc.at(1),vjet_disc.at(2),vjet_disc.at(3));
      }

      disc_="pfCombinedMVAV2BJetTags"; // FIXME HACK: do not know how to pass it from the python.
//      vector <float> vjet_disc;
      order_jets_by_disc(ev.jets_, disc_);
      vjet_disc.clear();
      for (auto & jet : ev.jets_) {
	  vjet_disc.push_back(jet.disc(disc_));
      }
      if (vjet_disc.at(2)<0.4432) return 0; // at least 4 medium CMSVA


      vector <bool> triggerBits      (npaths,0);
      vector <bool> triggerEmulations(npaths,0);
      if (or_paths_.size() != npaths) {
	  cout << "ERROR: or_paths_.size() != npaths" << endl;
	  return 0;
      }
      for (unsigned int ipath=0; ipath<npaths; ++ipath) {
	  triggerBits      [ipath]=(ev.eventInfo_.getFilter(or_paths_.at(ipath)));
	  float rand = r.Uniform(1);
	  triggerEmulations[ipath]=(rand<triggerEmulation[ipath]);
	  for (unsigned int isel=0; isel<nsel; ++isel) {
	      if (isel==0 && ipath) continue;
	      if (isel==0 || (isel==1&&triggerBits[ipath]) || (isel==2&&triggerEmulations[ipath])) {
		  h_jet3_pt [isel][ipath]->Fill(vjet_pt  .at(3));
		  h_jet3_CSV[isel][ipath]->Fill(vjet_disc.at(3));
		  h_HT      [isel][ipath]->Fill(sumpt          );
	      }
	  }
      }

      // OR of the two paths


      h2_triggerBits     .Fill(triggerBits      .at(0),triggerBits      .at(1),w);
      h2_triggerEmulation.Fill(triggerEmulations.at(0),triggerEmulations.at(1),w);

      return (triggerBits.at(0) || triggerBits.at(1));
    }

    virtual std::string get_name() {
      auto name = std::string{};
      name+= "trigger";
      for (const auto & or_path : or_paths_) name+="_OR"+or_path; 
      return name;
    }
    virtual bool output( TFile * tfile) {

      return true;
    }

    #include "trig.h"
};
