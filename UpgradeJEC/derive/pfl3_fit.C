
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH1F.h>

const int np=18;
double l3rsp[np][4]={
  //{17.3488,	0.749741,	0.010157,	0.0029907},
  //{22.3616,	0.776692,	0.011830,	0.00318862},
  {17.3488,	0.769741,	0.010157,	0.0029907},
  {22.3616,	0.786692,	0.011830,	0.00318862},
  {27.4012,	0.797747,	0.013360,	0.00325891},
  {32.4417,	0.818449,	0.014591,	0.00342929},
  {37.4059,	0.828972,	0.016138,	0.00362341},
  {42.4640,	0.836770,	0.017684,	0.00361193},
  {50.7286,	0.854519,	0.030483,	0.00258806},
  {64.0698,	0.863980,	0.039042,	0.00235095},

  {80.5561,	0.875852,	0.049890,	0.00228613},
  {103.947,	0.880087,	0.076824,	0.00187454},
  {134.290,       0.881353,	0.091183,	0.00195049},

  // {80.5561,	0.871852,	0.049890,	0.00228613},
  // {103.947,	0.879087,	0.076824,	0.00187454},
  // {134.290,     0.88153,	0.091183,	0.00195049},

  {173.306,	0.884511,	0.139963,	0.001627},
  {244.484,	0.889452,	0.257155,	0.00128926},
  {345.610,	0.893639,	0.341065,	0.00133698},
  {467.259,	0.897780,	0.554281,	0.00133149},
  {637.724,	0.904255,	0.879846,	0.00135506},
  {859.393,	0.910190,	1.382010,	0.00142224},
  {1194.04,	0.918904,	3.164200,	0.00131612}
};

double l3cor[np][4]={
  {16.9381,	1.33379,	0.072004,	0.00532047},
  {21.9662,	1.28751,	0.0988086,	0.00528573},
  {27.0328,	1.25353,	0.123967,	0.00512085},
  {32.1951,	1.22182,	0.145947,	0.00511942},
  {37.0755,	1.20631,	0.169713,	0.00527276},
  {41.7204,	1.19507,	0.188521,	0.00515854},
  {49.1930,	1.17025,	0.142932,	0.0035443 },
  {60.5918,	1.15743,	0.154297,	0.00314947},
  {74.9858,	1.14175,	0.174525,	0.00298016},
  {95.2348,	1.13625,	0.180012,	0.00242016},
  {121.523,	1.13462,	0.228079,	0.00251099},
  {155.644,	1.1357 ,	0.252660,	0.00209855},
  {219.157,	1.12429,	0.325363,	0.00162965},
  {310.532,	1.11902,	0.478065,	0.00167418},
  {420.691,	1.11386,	0.702770,	0.00165195},
  {577.366,	1.10588,	1.066490,	0.00165721},
  {782.972,	1.09867,	1.673930,	0.00171676},
  {1095.42,	1.08825,	3.449980,	0.00155867}
};


void pfl3_fit()
{
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);

  TGraphErrors *grsp = new TGraphErrors(np);
  grsp->SetName("grspL3Response");
  grsp->SetMarkerStyle(20);
  grsp->SetMarkerSize(1.2);
  grsp->SetMarkerColor(1);
  grsp->SetLineColor(1);

  TGraphErrors *gcor = new TGraphErrors(np);
  gcor->SetName("gcorL3Response");
  gcor->SetMarkerStyle(20);
  gcor->SetMarkerSize(1.2);
  gcor->SetMarkerColor(1);
  gcor->SetLineColor(1);

  TGraphErrors *gcormod = new TGraphErrors(np);
  gcormod->SetName("gcorL3Response_mod");
  gcormod->SetMarkerStyle(20);
  gcormod->SetMarkerSize(1.2);
  gcormod->SetMarkerColor(4);
  gcormod->SetLineColor(4);

  for(int i=0; i<np; i++){
    grsp->SetPoint(i, l3rsp[i][0], l3rsp[i][1]);
    grsp->SetPointError(i, l3rsp[i][2], l3rsp[i][3]);

    // l3cor[i][1] =1.0/l3rsp[i][1];
    // l3cor[i][3] =l3rsp[i][1]*l3rsp[i][1]*l3rsp[i][3];

    gcor->SetPoint(i, l3cor[i][0], l3cor[i][1]);
    gcor->SetPointError(i, l3cor[i][2], l3cor[i][3]);
    
    gcormod->SetPoint(i, l3cor[i][0], 1.0/l3rsp[i][1]);
    gcormod->SetPointError(i, l3cor[i][2], l3rsp[i][1]*l3rsp[i][1]*l3rsp[i][3]);
  }

  TF1 *fitrsp = new TF1("fitrsp","[0]-[1]/(pow(log10(x),2)+[2])-[3]*exp((-[4]*(log10(x)-[5])*(log10(x)-[5]))+([6]*(log10(x)-[5])))",
			//1.0,grsp->GetX()[grsp->GetN()-1]); 
			grsp->GetX()[0],grsp->GetX()[grsp->GetN()-1]);
  fitrsp->SetParameter(0,0.96);
  fitrsp->SetParameter(1,0.033);
  fitrsp->SetParameter(2,-0.7);
  fitrsp->SetParameter(3,0.02);
  fitrsp->SetParameter(4,1.02);
  fitrsp->SetParameter(5,2.7);
  fitrsp->SetParameter(6,0.016);
  fitrsp->SetLineWidth(2);
  grsp->Fit(fitrsp,"QR");

  TCanvas* crsp = new TCanvas("crsp","crsp",0,0,700,600);
  crsp->cd();
  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.05);
  gPad->SetBottomMargin(0.15);
  //if (logx) gPad->SetLogx();
  //if (logy) gPad->SetLogy();

  grsp->Draw("AP");
  TH1F* hrsp = grsp->GetHistogram();
  hrsp->SetXTitle("p_{T}^{REF}");
  hrsp->SetYTitle("p_{T}/p_{T}^{REF}");
  grsp->SetLineWidth(2);
  grsp->SetMarkerStyle(kFullCircle);

  string fitcor_as_str;
  fitcor_as_str = "[0]+[1]/(pow(log10(x),2)+[2])+[3]*exp((-[4]*(log10(x)-[5])*(log10(x)-[5]))+([6]*(log10(x)-[5])))";
  TF1* fitcor = new TF1("fitcor",fitcor_as_str.c_str(),
			gcor->GetX()[1],gcor->GetX()[gcor->GetN()-1]);
  //gcor->GetX()[0],gcor->GetX()[15]);

  // fitcor->SetParameter(0,1.04);
  // fitcor->SetParameter(1,0.011);
  // fitcor->SetParameter(2,-0.10);
  // fitcor->SetParameter(3,0.010);
  // fitcor->SetParameter(4,0.020);
  // fitcor->SetParameter(5,1.10);
  // fitcor->SetParameter(6,0.010);

  fitcor->SetParameter(0,1.04);
  fitcor->SetParameter(1,.033);
  fitcor->SetParameter(2,-0.7);
  fitcor->SetParameter(3,0.02);
  fitcor->SetParameter(4,1.02);
  fitcor->SetParameter(5,2.7);
  fitrsp->SetParameter(6,0.016);

  fitcor->SetLineWidth(2);
  gcor->Fit(fitcor,"QR");


  for(int i=0; i<np; i++){
    cout << "Ref pT : " << l3rsp[i][0] << "\t Rec pT : "<< l3cor[i][0] << "\t Corr fac : "<< fitcor->Eval(l3cor[i][0]) << endl;
  }

  TCanvas* ccor = new TCanvas("ccor","ccor",715,0,700,600);
  ccor->cd();
  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.05);
  gPad->SetBottomMargin(0.15);
  //if (logx) gPad->SetLogx();
  //if (logy) gPad->SetLogy();

  gcor->Draw("AP");
  gcormod->Draw("PSAME");
  TH1F* hcor = gcor->GetHistogram();
  hcor->SetXTitle("p_{T}");
  hcor->SetYTitle("L3 correction");
  gcor->SetLineWidth(2);
  gcor->SetMarkerStyle(kFullCircle);
  

}
