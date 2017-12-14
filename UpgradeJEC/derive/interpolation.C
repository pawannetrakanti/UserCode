// Declare Global Variables, 
// which can be referred from anywhere after their declaration: 
const Long64_t n = 4;
Double_t Point_x[n] = {0, 0.5, 1, 2};
Double_t Point_y[n] = {1, 0.368, 0.135, 0.018};

Double_t LagPol(Double_t *x, Double_t *par) ;
Double_t LogChi2(Double_t *x, Double_t *par) ;

// main function with the same name of this file
// it will be executed when you type ".x interpolation.C" in ROOT
void interpolation()
{
// Section 1. Draw the points on a canvas
   TCanvas *c1 = new TCanvas("c1","interpolation",0,0,1000,800);

   TGraph *g1 = new TGraph(n,Point_x,Point_y);

   g1->SetMarkerStyle(20);
   g1->SetMarkerSize(2);
   
   g1->GetXaxis()->SetLimits(-1, 3);        // set real range
   g1->GetXaxis()->SetRangeUser(-0.5, 2.5); // set visible range
   g1->GetXaxis()->SetTitle("X");
   g1->GetXaxis()->CenterTitle();

   g1->GetYaxis()->SetLimits(-1, 2.0);
   g1->GetYaxis()->SetRangeUser(-0.4, 1.6);
   g1->GetYaxis()->SetTitle("Y");
   g1->GetYaxis()->CenterTitle();

   g1->Draw("ap"); // options to draw a graph are described on 
                   // http://root.cern.ch/root/html/TGraph.html#TGraph:PaintGraph
   

// Section 2. Draw the Lagrange Polynomial to the same canvas
   TF1 *lag = new TF1("lag", LagPol, -0.5, 2.5, 0); 
   // refer to http://root.cern.ch/root/html/TF1.html#TF1:TF1 for the usage of TF1
   // function "LagPol" is defined outside the main function

   lag->SetLineColor(kBlue); // you can find several predefined color variables on
                             // http://root.cern.ch/root/html/Gtypes.h
   lag->Draw("lsame"); // options to draw a function are described on
                       // http://root.cern.ch/root/html/TF1.html#TF1:Draw


// Section 3. Draw the Cubic Spline to the same canvas
   TSpline3 *sp = new TSpline3("Cubic Spline", Point_x, Point_y, n, "b2e2", 0, 0);
   // refer to http://root.cern.ch/root/html/Tspline3.html for the usage of TSpline3
   // "b2e2" together with the last two "0" means that the second derivatives 
   // of the begin and end points equal to zero
   sp->SetLineColor(kRed);
   sp->Draw("lsame");


// Section 4. Use A*e^(-Bx) to fit the points
   TF1 *ex = new TF1("ex","[0]*exp(-[1]*x)",-3,3);
   ex->SetLineColor(kGreen);
   g1->Fit("ex", "W", "lsame", -0.5, 2.5); 
   // refer to http://root.cern.ch/root/html/TGraph.html#TGraph:Fit


// Section 5. Draw the legends for the curves and points
   TLegend* leg = new TLegend(0.5, 0.6, 0.85, 0.8);
   leg->AddEntry(g1,"Experiment Points","p");
   leg->AddEntry(lag,"Lagrange Polynomial","l");
   leg->AddEntry(sp,"Cubic Spline","l");
   leg->AddEntry(ex,"Fit Function","l");
   leg->Draw();


// Section 6. Draw lg(chi^2) on another canvas
   TCanvas *c2 = new TCanvas("c2","lg(Chi^2) Distribution",600,0,1000,800);
   TF2 *ch = new TF2("ch", LogChi2, 0.9, 1.1, 1.9, 2.1, 0);

   ch->GetXaxis()->SetTitle("par 1");
   ch->GetXaxis()->CenterTitle();
   
   ch->GetYaxis()->SetTitle("par 2");
   ch->GetYaxis()->CenterTitle();

   ch->GetZaxis()->SetTitle("lg(#chi^{2})");
   ch->GetZaxis()->CenterTitle();

   ch->Draw("surf1");
}

//____________________________________________________________________________
Double_t LagPol(Double_t *x, Double_t *par) 
{
   Double_t Lagrange, term;
   
   Lagrange = 0;
   for (Long64_t i=0; i<n; i++) {
      term = 1;
      for (Long64_t j=0; j<n; j++) {
         if (j != i) 
            term = term*(x[0]-Point_x[j])/(Point_x[i]-Point_x[j]);
      }
      Lagrange = Lagrange + Point_y[i]*term;
   }
   return Lagrange;
}

//____________________________________________________________________________
Double_t LogChi2(Double_t *x, Double_t *par) 
{
   Double_t ChiSqr = 0;
   for (Long64_t i=0; i<n; i++)
     ChiSqr = ChiSqr + pow(Point_y[i] - x[0]*exp(-x[1]*Point_x[i]),2);
   return log(ChiSqr);
}
