/*
 * TriggerEffStudiesPlotter.C
 *
 *  Created on: Mar 29, 2017
 *      Author: rossin
 */


#include <TH2.h>
#include <TGraph.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPad.h>
#include <TString.h>
#include <TFile.h>
#include <TList.h>
#include <TCollection.h>
#include <TObject.h>
#include <vector>
#include <vector>

#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <algorithm>



void TriggerEffStudiesPlotter(TString sfile="HHTo4B_SM.root") {
  TString sdir("/data/rossin/hh2016/");
  TFile * fIn = new TFile(sdir+sfile);

  const static int nvars  = 3;
  const static int nsel   = 3;
  const static int npaths = 2;
  TString sselection[nsel  ]={"_All","_HLTbits","_HLTemulation"};
  TString spaths    [npaths]={"_QuadJet45","_DoubleJet90_Double30"};
//  TString vars      [nvars ]={"jet3_pt","jet3_CSV","HT"};
  TH1F* h_jet3_pt [nsel][npaths];
  TH1F* h_jet3_CSV[nsel][npaths];
  TH1F* h_HT      [nsel][npaths];
  TH2F* h2_triggerEmulation;
  TH2F* h2_triggerBits     ;

  TIter next(gDirectory->GetListOfKeys());
  TH1 *h;
  TObject* obj;
  while((obj= (TObject*)next())){
      TObject* o = fIn->FindObjectAny(obj->GetName());
      if (o->InheritsFrom(TH1::Class())==false) continue;
      if (TString(o->GetName()).Contains("h2_triggerBits"     )) h2_triggerBits      = (TH2F*) o;
      if (TString(o->GetName()).Contains("h2_triggerEmulation")) h2_triggerEmulation = (TH2F*) o;
      for (unsigned int isel=0; isel<nsel; ++isel) {
	  for (unsigned int ipath=0; ipath<npaths; ++ipath) {
	      if (isel==0 && ipath) continue;
	      if (TString(o->GetName()).Contains(TString("h_jet3_pt" )+(isel ? spaths[ipath] : TString(""))+sselection[isel])) h_jet3_pt [isel][ipath] = (TH1F*) o;
	      if (TString(o->GetName()).Contains(TString("h_jet3_CSV")+(isel ? spaths[ipath] : TString(""))+sselection[isel])) h_jet3_CSV[isel][ipath] = (TH1F*) o;
	      if (TString(o->GetName()).Contains(TString("h_HT"      )+(isel ? spaths[ipath] : TString(""))+sselection[isel])) h_HT      [isel][ipath] = (TH1F*) o;
	  }
      }
  }

  TCanvas* cEff_jet3_pt = new TCanvas("cjet3_pt" ,"cjet3_pt" ,0,0,900,1200);
  cEff_jet3_pt->Divide(1,2);
  TCanvas* cEff_jet3_CSV= new TCanvas("cjet3_CSV","cjet3_CSV",20,20,900,1200);
  cEff_jet3_CSV->Divide(1,2);
  TCanvas* cEff_HT      = new TCanvas("cHT"      ,"cHT"      ,40,40,900,1200);
  cEff_HT->Divide(1,2);
  unsigned int ipad=1;
  for (unsigned int ipath=0; ipath<npaths; ++ipath) {
      bool isFirst=true;
      for (unsigned int isel=1; isel<nsel; ++isel) {
	      if (isel==0 && ipath) continue;
	      h_jet3_pt [isel][ipath]->Divide(h_jet3_pt [0][0]);
	      h_jet3_pt [isel][ipath]->SetLineColor(isel);
	      h_jet3_pt [isel][ipath]->SetMarkerStyle(20);
	      h_jet3_pt [isel][ipath]->SetMarkerColor(isel);
	      h_jet3_CSV[isel][ipath]->Divide(h_jet3_CSV[0][0]);
	      h_jet3_CSV[isel][ipath]->SetLineColor(isel);
	      h_jet3_CSV[isel][ipath]->SetMarkerStyle(20);
	      h_jet3_CSV[isel][ipath]->SetMarkerColor(isel);
	      h_HT      [isel][ipath]->Divide(h_HT      [0][0]);
	      h_HT      [isel][ipath]->SetLineColor(isel);
	      h_HT      [isel][ipath]->SetMarkerStyle(20);
	      h_HT      [isel][ipath]->SetMarkerColor(isel);
	      if (isFirst) {
		  isFirst=false;
		  cEff_jet3_pt->cd(1+ipath);
		  h_jet3_pt [isel][ipath]->Draw();
		  cEff_jet3_CSV->cd(1+ipath);
		  h_jet3_CSV[isel][ipath]->Draw();
		  cEff_HT      ->cd(1+ipath);
		  h_HT      [isel][ipath]->Draw();
	      }
	      else {
		  cEff_jet3_pt->cd(1+ipath);
		  h_jet3_pt [isel][ipath]->Draw("same");
		  cEff_jet3_CSV->cd(1+ipath);
		  h_jet3_CSV[isel][ipath]->Draw("same");
		  cEff_HT      ->cd(1+ipath);
		  h_HT      [isel][ipath]->Draw("same");
	      }
	  }
  }
  cEff_jet3_pt ->SaveAs(sdir+cEff_jet3_pt ->GetName()+TString(".png"));
  cEff_jet3_CSV->SaveAs(sdir+cEff_jet3_CSV->GetName()+TString(".png"));
  cEff_HT      ->SaveAs(sdir+cEff_HT      ->GetName()+TString(".png"));

//  if (h2_triggerBits!=0) h2_triggerBits->Draw("colz0");

  return;
}
