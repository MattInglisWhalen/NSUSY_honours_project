
//This version plots stop mass vs. gluino mass
// with neutralino mass as the z dimension

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>

#include <TFile.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TH1.h>

#include <TMath.h>
#include <Math/DistFunc.h>

//holds search data
#include "MIW_search_data.h"
//holds global variables and many helper functions
#include "MIW_plot_helper.h"

#define PRINT_EXCL_V_BR 0
#define STOP_PAIR_PRODUCTION 1
#define GLUINO_PAIR_PRODUCTION 2

#define GO_TTXBBXN1 1
#define GO_T1TX 2
#define T1_TN1_BX1 3

#define PRINT_CONTOUR 1
#define PRINT_CONTOUR_TMODE 0
#define PRINT_CONTOUR_BMODE 0
#define PRINT_CONTOUR_BOTH 0
#define PRINT_CONTOUR_GEOM 0
#define PRINT_POINTS 0

#define EXACT 0
#define TMODE 1
#define BMODE 2
#define BOTH  3
#define GEOM  4

using namespace std;

void plot_excl_br_curves_GOMED_T1_MU(char* stordir,char* ansystag,char* ztag=""){

  cout << "Processing " << ansystag << endl;

  get_search_data(ansystag);

  read_data(stordir,ansystag,ztag,0,GLUINO,STOP1,NEUT1);
  read_data(stordir,ansystag,ztag,1,GLUINO,STOP1,NEUT1);
  read_data(stordir,ansystag,ztag,2,GLUINO,STOP1,NEUT1);

  sort_xyz();

  Double_t sumn1=0.,meann1,sumx1=0.,meanx1;
  for(Int_t i=0;i<(signed)GP.size();++i){
    sumn1+=GP[i].mpart(NEUT1);
    sumx1+=GP[i].mpart(CHRG1);
  }
  meann1=sumn1/GP.size();
  meanx1=sumx1/GP.size();

  refine();
  refine();
//  refine();  //this one might be too much for ROOT to handle

  sort_zyx();

// ------------------------------------------------ making plots ---------------------------------------------------- //

  string* exclfilename=new string(ansystag);  exclfilename->append("_exclusion.root");
  TFile* exclfile = new TFile(exclfilename->c_str(),"RECREATE");

  Int_t contour_colours[10]={1,2,8,4,6,7,3,9,12,46};  //See TAttMarker on the ROOT website

  string* canvname=new string("GOMED_T1_MU");
  string* framename=new string("pp#rightarrow#tilde{g}#tilde{g} , #tilde{g}#rightarrow#tilde{t}_{1}#bar{t}");
  framename->append(" , #tilde{#chi}_{1}^{0}~#tilde{H}^{0} , #tilde{#chi}_{1}^{#pm}~#tilde{H}^{#pm}");
  framename->append(" , m(#tilde{#chi}_{1}^{0}) #approx ");
  framename->append( mass_string( floor(meann1/10+0.5)*10 ) );
  framename->append(", m(#tilde{#chi}_{1}^{0})-m(#tilde{#chi}_{1}^{#pm}) < 5 GeV ");


  //create the canvas needed for this particular neutralino mass
  TCanvas* canv = new TCanvas(canvname->c_str(),"Exclusion Plot",300,325);
  canv->SetTickx();
  canv->SetTicky();
 // canv->SetGridx();
 // canv->SetGridy();

  string* pdfname=new string(ztag); pdfname->append("_");
  pdfname->append(ansystag); pdfname->append("_GOMED_T1_MU.pdf");

  TH1F* frame=new TH1F("frame",framename->c_str(),1000,200,1200);
  frame->SetStats(0);
  frame->SetMinimum(0);
  frame->SetMaximum(1300);
  frame->GetXaxis()->SetTitle("#tilde{g}   /   GeV");
  frame->GetXaxis()->SetTitleSize(0.035);
  frame->GetXaxis()->SetTickLength(0.02);
  frame->GetXaxis()->SetLabelSize(0.025);
  frame->GetXaxis()->SetLimits(700,1500); //minx/maxx limits
  frame->GetYaxis()->SetTitle("#lower[-0.3]{#tilde{t_{1}}   /   GeV}");
  frame->GetYaxis()->SetTitleSize(0.035);
  frame->GetYaxis()->SetTickLength(0.02);
  frame->GetYaxis()->SetLabelSize(0.025);
  frame->Draw(" ");

  //make the lines denoting kinematically disallowed regions
  TGraph* kinemat_stop = new TGraph(2);
  kinemat_stop->SetPoint(0,700,700-TOPMASS);
  kinemat_stop->SetPoint(1,1500,1500-TOPMASS);
  kinemat_stop->SetLineWidth(3);
  kinemat_stop->SetLineColor(13);
  kinemat_stop->SetLineStyle(9);
  kinemat_stop->Draw("L");

  TLatex* go_tt1 = new TLatex(800,700,"#color[13]{#font[22]{#tilde{g} #rightarrow t#tilde{t}_{1}   forbidden}}");
  go_tt1->SetTextSize(0.03);
  go_tt1->SetTextAngle(35);
  go_tt1->Draw();

  //make the lines denoting kinematically disallowed regions
  TGraph* kinemat_tn1 = new TGraph(2);
  kinemat_tn1->SetPoint(0,700,TOPMASS+meann1);
  kinemat_tn1->SetPoint(1,1500,TOPMASS+meann1);
  kinemat_tn1->SetLineWidth(2);
  kinemat_tn1->SetLineColor(13);
  kinemat_tn1->SetLineStyle(9);
  kinemat_tn1->Draw("L");
  
  TLatex* t1_tn1 = new TLatex(1260,TOPMASS+meann1-50,"#color[13]{#font[22]{#tilde{t}_{1} #rightarrow t#tilde{#chi}_{1}^{0}   forbidden}}");
  t1_tn1->SetTextSize(0.02);
  t1_tn1->SetTextAngle(0);
  t1_tn1->Draw();

  //make the lines denoting kinematically disallowed regions
  TGraph* kinemat_bx1 = new TGraph(2);
  kinemat_bx1->SetPoint(0,700,BOTMASS+meanx1);
  kinemat_bx1->SetPoint(1,1500,BOTMASS+meanx1);
  kinemat_bx1->SetLineWidth(2);
  kinemat_bx1->SetLineColor(13);
  kinemat_bx1->SetLineStyle(9);
  kinemat_bx1->Draw("L");

  if(BOTMASS+meanx1 > 50 ){
    TLatex* t1_bx1 = new TLatex(850,BOTMASS+meanx1-50,"#color[13]{#font[22]{#tilde{t}_{1} #rightarrow b#tilde{#chi}_{1}^{+}   forbidden}}");
    t1_bx1->SetTextSize(0.02);
    t1_bx1->SetTextAngle(0);
    t1_bx1->Draw();
  }

  //make a legend please (for the BR line colours)
  TLegend* leg = new TLegend(0.1,0.7,0.3,0.9);
  leg->SetHeader("#lower[0.2]{BR( #tilde{t}_{1}#rightarrowt#tilde{#chi}_{1}^{0} )}");

  Int_t itt=0,legend_on;  //for the BR contour
  //for each branching ratio
  for(Double_t iBR=plot_opts[0];(iBR<=1.)&&(iBR>=0.);iBR=iBR+plot_opts[1]){

    legend_on=0;

    exclGP.clear(); exclGP_tmode.clear(); exclGP_bmode.clear(); exclGP_both.clear(); exclGP_geom.clear();
    ordGP.clear(); ordGP_tmode.clear(); ordGP_bmode.clear(); ordGP_both.clear(); ordGP_geom.clear();

    find_excluded_GPs(iBR,GLUINO_PAIR_PRODUCTION,GO_T1TX,T1_TN1_BX1,EXACT);

    string* BRtitle=new string(BR_string(iBR));
    cout << "BR_T " << iBR << ":  # Excluded--> " << exclGP.size() << "   # Not Excluded--> " << GP.size()-exclGP.size() << endl;

    if(PRINT_CONTOUR){
      ordGP=order_excluded_GPs(exclGP);
      TGraph* gcontour=new TGraph(ordGP.size());
      for(Int_t j=0;j<(signed)ordGP.size();++j){
        gcontour->SetPoint(j,ordGP[j].x(),ordGP[j].y());
      }
      if( ordGP.size()>0){
        gcontour->SetName("Contour");
        gcontour->SetLineStyle(contour_colours[itt]);
        gcontour->SetLineColor(contour_colours[itt]);
        gcontour->SetLineWidth(2+2*legend_on);
        gcontour->Draw("L");
        if(!legend_on){
          leg->AddEntry(gcontour,BRtitle->c_str(),"L");
          legend_on=1;
        }
      }
    }
    if(PRINT_CONTOUR_TMODE){
      find_excluded_GPs(iBR,GLUINO_PAIR_PRODUCTION,GO_T1TX,T1_TN1_BX1,TMODE);
      ordGP_tmode=order_excluded_GPs(exclGP_tmode);
      TGraph* gcontour_tmode=new TGraph(ordGP_tmode.size());
      for(Int_t j=0;j<(signed)ordGP_tmode.size();++j){
        gcontour_tmode->SetPoint(j,ordGP_tmode[j].x(),ordGP_tmode[j].y());
      }
      if( ordGP_tmode.size()>0){
        gcontour_tmode->SetName("Contour");
        gcontour_tmode->SetLineStyle(contour_colours[itt]);
        gcontour_tmode->SetLineColor(contour_colours[itt]);
        gcontour_tmode->SetLineWidth(2+2*legend_on);
        gcontour_tmode->Draw("L");
        if(!legend_on){
          leg->AddEntry(gcontour_tmode,BRtitle->c_str(),"L");
          legend_on=1;
        }
      }
    }
    if(PRINT_CONTOUR_BMODE){
      find_excluded_GPs(iBR,GLUINO_PAIR_PRODUCTION,GO_T1TX,T1_TN1_BX1,BMODE);
      ordGP_bmode=order_excluded_GPs(exclGP_bmode);
      TGraph* gcontour_bmode=new TGraph(ordGP_bmode.size());
      for(Int_t j=0;j<(signed)ordGP_bmode.size();++j){
        gcontour_bmode->SetPoint(j,ordGP_bmode[j].x(),ordGP_bmode[j].y());
      }
      if( ordGP_bmode.size()>0){
        gcontour_bmode->SetName("Contour");
        gcontour_bmode->SetLineStyle(contour_colours[itt]);
        gcontour_bmode->SetLineColor(contour_colours[itt]);
        gcontour_bmode->SetLineWidth(2+2*legend_on);
        gcontour_bmode->Draw("L");
        if(!legend_on){
          leg->AddEntry(gcontour_bmode,BRtitle->c_str(),"L");
          legend_on=1;
        }
      }
    }
    if(PRINT_CONTOUR_BOTH){
      find_excluded_GPs(iBR,GLUINO_PAIR_PRODUCTION,GO_T1TX,T1_TN1_BX1,BOTH);
      ordGP_both=order_excluded_GPs(exclGP_both);
      TGraph* gcontour_both=new TGraph(ordGP_both.size());
      for(Int_t j=0;j<(signed)ordGP_both.size();++j){
        gcontour_both->SetPoint(j,ordGP_both[j].x(),ordGP_both[j].y());
      }
      if( ordGP_both.size()>0){
        gcontour_both->SetName("Contour Both");
        gcontour_both->SetLineStyle(contour_colours[itt]);
        gcontour_both->SetLineColor(contour_colours[itt]);
        gcontour_both->SetLineWidth(2+2*legend_on);
        gcontour_both->Draw("L");
      }
    }
    if(PRINT_CONTOUR_GEOM){
      find_excluded_GPs(iBR,GLUINO_PAIR_PRODUCTION,GO_T1TX,T1_TN1_BX1,GEOM);
      ordGP_geom=order_excluded_GPs(exclGP_geom);
      TGraph* gcontour_geom=new TGraph(ordGP_geom.size());
      for(Int_t j=0;j<(signed)ordGP_geom.size();++j){
        gcontour_geom->SetPoint(j,ordGP_geom[j].x(),ordGP_geom[j].y());
      }
      if( ordGP.size()>0){
        gcontour_geom->SetName("Contour");
        gcontour_geom->SetLineStyle(contour_colours[itt]);
        gcontour_geom->SetLineColor(contour_colours[itt]);
        gcontour_geom->SetLineWidth(2+2*legend_on);
        gcontour_geom->Draw("L");
        if(!legend_on){
          leg->AddEntry(gcontour_geom,BRtitle->c_str(),"L");
          legend_on=1;
        }
      }
    }
    if(PRINT_POINTS){
      TGraph* gexcluded=new TGraph(exclGP.size());
      for(Int_t j=0;j<(signed)exclGP.size();++j){
        gexcluded->SetPoint(j,exclGP[j].x(),exclGP[j].y());
      }

      if( exclGP.size()>0){
        gexcluded->SetName("Excluded");
        gexcluded->SetMarkerColor(contour_colours[itt]);  
        gexcluded->SetMarkerSize(1.5);   
        gexcluded->SetMarkerStyle(7);   // cross
        gexcluded->SetDrawOption("P");
        gexcluded->Draw("P");
        leg->AddEntry(gexcluded,BRtitle->c_str(),"P");
      }
    }

    ++itt;

    delete BRtitle;
  } //for various branching ratios
  leg->Draw();

  canv->SaveAs(pdfname->c_str());
  canv->Write();
  exclfile->Write();
  canv->Clear();

  delete exclfilename;
  delete exclfile;
  delete framename;
  delete pdfname;

};
