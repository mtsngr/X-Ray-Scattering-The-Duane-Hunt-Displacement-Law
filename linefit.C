// Line fitting macro written for ROOT
{

//       Data codes for the accelerating potentials

   TCanvas *c1 = new TCanvas("c1","A Simple Graph with error bars",200,10,700,500);


   c1->SetGrid();
   c1->GetFrame()->SetBorderSize(12);
    gStyle->SetOptFit(11111);

      const int ndata = 99;

  float x[ndata] = {3.0,	3.1,	3.2,	3.3,	3.4,	3.5,	3.6,	3.7,	3.8,	3.9,	4.0,	4.1,	4.2,	4.3,	4.4,	4.5,	4.6,	4.7,	4.8,	4.9,	5.0,	5.1,	5.2,	5.3,	5.4,	5.5,	5.6,	5.7,	5.8,	5.9,	6.0,	6.1,	6.2,	6.3,	6.4,	6.5,	6.6,	6.7,	6.8,	6.9,	7.0,	7.1,	7.2,	7.3,	7.4,	7.5,	7.6,	7.7,	7.8,	7.9,	8.0,	8.1,	8.2,	8.3,	8.4,	8.5,	8.6,	8.7,	8.8,	8.9,	9.0,	9.1,	9.2,	9.3,	9.4,	9.5,	9.6,	9.7,	9.8,	9.9,	10.0,	10.1,	10.2,	10.3,	10.4,	10.5,	10.6,	10.7,	10.8,	10.9,	11.0,	11.1,	11.2,	11.3,	11.4,	11.5,	11.6,	11.7,	11.8,	11.9,	12.0,	12.1,	12.2,	12.3,	12.4,	12.5,	12.6,	12.7,	12.8};


  float y[ndata] = {37.0,	25.0,	20.0,	18.0,	10.0,	7.0,	4.0,	4.0,	2.0,	1.0,	2.0,	1.0,	1.0,	1.0,	1.0,	3.0,	2.0,	2.0,	0.0,	1.0,	1.0,	2.0,	1.0,	2.0,	2.0,	0.0,	1.0,	1.0,	1.0,	1.0,	1.0,	0.0,	1.0,	1.0,	1.0,	0.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	0.0,	2.0,	0.0,	2.0,	0.0,	2.0,	1.0,	1.0,	1.0,	0.0,	2.0,	2.0,	4.0,	2.0,	6.0,	5.0,	13.0,	13.0,	17.0,	16.0,	20.0,	27.0,	25.0,	33.0,	25.0,	23.0,	28.0,	27.0,	31.0,	29.0,	29.0,	29.0,	33.0,	28.0,	34.0,	36.0,	30.0,	32.0,	33.0,	31.0,	34.0,	27.0,	32.0,	32.0,	31.0,	31.0,	33.0,	36.0,	30.0,	30.0,	28.0,	26.0,	30.0,	32.0,	31.0};

  float sx[ndata]= {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  float sy[ndata]= {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};



  float weight = 0;
  float totw = 0; // Total weight
  float xybar = 0, xbar = 0, ybar = 0, x2bar = 0; // weighted averages

  for (int i=0; i<ndata; ++i) {
    weight = 1./(sy[i]*sy[i]);
    totw += weight;
    xybar += (x[i]*y[i]*weight);
    xbar += (x[i]*weight);
    ybar += (y[i]*weight);
    x2bar += (x[i]*x[i]*weight);
  }

  xybar /= totw;
  xbar /= totw;
  ybar /= totw;
  x2bar /= totw;

  float sy2bar = ndata / totw; // weighted average error squared
  float slope = (xybar - xbar*ybar) / (x2bar - xbar*xbar);
  float itcpt = ybar - slope * xbar;
  float slopeerr = sqrt ( sy2bar / (ndata * (x2bar - xbar*xbar) ) );
  float itcpterr = sqrt ( x2bar ) * slopeerr;
  cout << "slope of fit line = " << slope << " +- " << slopeerr << endl;
  cout << "intercept of fit line = " << itcpt << " +- " << itcpterr << endl;


  TGraphErrors *mygraph = new TGraphErrors(ndata,x,y,sx,sy);
  mygraph->Draw("A*");

  /*
  TF1 *ffitline = new TF1("ffitline","[0]*x+[1]",0,6); // x ekseninin aralığı [0,6]
  ffitline->SetParameter(0,slope);
  ffitline->SetParameter(1,itcpt);
  ffitline->SetLineColor(kBlue); // draw in blue color
  ffitline->SetLineStyle(2); // draw dotted line
  //ffitline->Draw("same");
  */

  mygraph->SetTitle("; Theta (degree); Accumulation Rate (Imp/s) ");


  // let's also try the same with ROOT's own fitter
  TF1 *fnew = new TF1("fnew","[0]*x+[1]",0,6);
  fnew->SetParameters(3.1416,2.7182); // arbitrary starting parameters
  mygraph->Fit("pol2", "F");
  mygraph->Fit(fnew);


  // TF1 *fpol = mygraph->GetFunction("pol2", 8.5, 9.3);
  // fpol->SetLineWidth(1);


  //mygraph->Draw("A*");
  

  gStyle->SetOptStat(111110);
  gStyle->SetOptFit(1111);

  mygraph->SetTitle("The intensity of reflected X-Rays as a function of theta for V = 30 kV; Theta (degree); Accumulation Rate (Imp/s) ");   


   mygraph->SetMarkerColor(1);
   mygraph->SetMarkerStyle(3);
   mygraph->Draw("ALP");

   c1->Update();
  


}





