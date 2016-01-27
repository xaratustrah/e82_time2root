//////////////////////////////////////////////////////////////////////////
// 								        //
// This program reads IQT Tektronix files and after applying multitaper //
// spectrum estimation shows the resutls using root.		        //
// 								        //
// IQT interface and multitaper: F. Nolden @ GSI		        //
// root interface: N. Winckler @ GSI				        //
// 								        //
// 2010-2012							        //
// GSI Darmstadt						        //
// 								        //
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>
#include <stdlib.h>
#include <TH2D.h>
#include <TRootCanvas.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TColor.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TSystem.h>
#include <gsl/gsl_fft_complex.h>

#include "iqtdata.h"
#include "FritzDPSS.h"

#include <iostream>
#include <sstream> // for ostringstream
#include <string>

using namespace std;

int main(int argc, char *argv[])
{
     IQData__TekRSA3303B my_iq_data;
     string filename, file_basename, iq_filename, out_filename;

     bool draw=false;


     int ip;
     //cout << "Please enter name of *.iqt-file without extension: ";
     //cin >> filename;
     if (argc == 2) draw = false;
     else if (argc == 3){
	  draw = true;
	  out_filename = argv[2];
	  gROOT->SetBatch(kTRUE);
     }
     else{
	  cout << "Too many arguments." << endl;
	  cout << "Please provice either iqt filename for visual display, or iqt filename AND output filename for storing in a PNG file in batch mode." << endl;
	  cout << "Aborting..." << endl;
	  
	  exit(0);
     }
     
     filename=argv[1];
     file_basename = gSystem->TSystem::BaseName(argv[1]);
     iq_filename = filename; //+ ".iqt";

     cout << "Data " << iq_filename << " is going to be read" << endl;
     if (!my_iq_data.ReadFile(iq_filename))
     {
	  cout << "Program has been stopped." << endl;
	  exit(0);
     }
     double delta_t;

     if(!my_iq_data.GetDeltaT(delta_t)) 
     {
	  cout << "Error during the check of time interval" << endl;
	  cout << "Program has been stopped." << endl;
	  exit(0);
     }
     int valid_frames = my_iq_data.ValidFrames();
 
     int block_size = my_iq_data.BlockSize();

     if (valid_frames == block_size) 
     {
	  cout << "All possible " << block_size
	       << " Frames with each 1024 Points are valid" << endl;
     }
     else 
     {
	  cout << "From the " << block_size
	       << "  Frames with each 1024 Points, only "
	       << valid_frames << " are valid." << endl;
	  cout << "Program has been stopped." << endl;
	  exit(0);
     }
     cout << "Es gibt also insgesamt " << block_size*1024 << " gemessene Punkte" << endl;


     ///--------------FFT-PARAMETERS---->
     int ifreq;
     int Npoint = 1024; // Number of points for making one FFT frame
     double center_frequency = my_iq_data.CenterFrequency();
     int first_index = 0;
     complex <double> * iq_data_p = new complex <double> [Npoint];

     ///-----------------------FFT---->
     /// for Multitaper, delta_t cover on RSA always 1024 Points
     double dtm = (double) Npoint / 1024. * delta_t; // Multitaper time frame
     double dfreq_multi =  1. / dtm; // Nyquist!
     cout<<"1 frame in frequency="<<dfreq_multi<<" Hz"<<endl;
     double start_frequency = center_frequency - Npoint * dfreq_multi / 2.;
     int ifmin_multi = (Npoint * 3) / 16;
 
     int ifmax_multi = 1 + (Npoint * 13) / 16;
     FritzDPSS my_dpss(Npoint, 4);
     my_dpss.GetFromDisk();
 
     int n_of_multi = block_size;
     int i_multi;
     double ** multi_taper_spektrum = new double * [n_of_multi]; 
     for (i_multi=0; i_multi < n_of_multi; i_multi++) multi_taper_spektrum[i_multi] = new double [Npoint];

     i_multi=0;
     for (first_index=0; first_index<=1024*(n_of_multi-1); first_index+=1024)
     {
	  int index = first_index;
	  complex <double> iq_value;

	  for (ip=0; ip < Npoint; ip++) 
	  {
	       if (!my_iq_data.GetIQ(iq_value, index))
	       {
		    cout << "Error retrieving the value at index " << index << endl;
		    exit(0);
 
	       }
	       index++;
	       iq_data_p[ip] = iq_value; 
	  }

	  // Spectrum calculation
 
	  my_dpss.GetSpectrum(iq_data_p, multi_taper_spektrum[i_multi]);

	  // progress bar
	  if (!(i_multi % 10)) {
	       printf("processing... %5.2f%%\r", (double)i_multi/n_of_multi*100);
	       fflush(stdout);
	  }

	  i_multi++;
     }

     cout << endl;

     int binf=ifmax_multi-ifmin_multi;
     double fstart=start_frequency+(double)ifmin_multi*dfreq_multi;
     double fend= start_frequency+(double)(ifmax_multi-1)*dfreq_multi;
 
     int bint=n_of_multi;
     TString histoname="Schottky_Spectrum_"+file_basename;
     TH2D *RAW = new TH2D(histoname,histoname,binf,fstart,fend,bint,0.0,(double)(bint)*delta_t);

     // constant of multiplication for z axis
     double constant = 1.0;
     for (ifreq = ifmin_multi; ifreq < ifmax_multi; ifreq++)
     {
	  for (i_multi = 0; i_multi < n_of_multi; i_multi++)
	  {
	       RAW->SetBinContent(ifreq-ifmin_multi+1,i_multi+1, multi_taper_spektrum[i_multi][ifreq] * constant);
	  } 
     }

     /// DIPLAY SPECTRUM //////////////////////////////////

     // custom palette specially good for number of particles,
     // or where you have many and few particles in the same spectrum,
     // in which case log scale is also useful -- D. Shubina 2011

     Int_t MyPalette[100];
     Double_t r[]    = {1., .9, .0, 0.0, 0.0};
     Double_t g[]    = {1., .9, 1.0, .0, .0};
     Double_t b[]    = {1., .9, .8, 1.0, .0};
     Double_t stop[] = {0., .25, .50, .75, 1.0};
     Int_t FI = TColor::CreateGradientColorTable(5, stop, r, g, b, 100);
     for (int i=0;i<100;i++) MyPalette[i] = FI+i;

     //gStyle->SetPalette(100, MyPalette); // use the custom palette above or one of the standard ones below
     gStyle->SetPalette(1); // Standard color scheme
     //gStyle->SetPalette(20); // Jet Color scheme
     //gStyle->SetPalette(500); // Blue Color scheme

     gStyle->SetTitleFontSize(0.03);

     gStyle->SetCanvasBorderMode(0);
     gStyle->SetCanvasColor(0);
     gStyle->SetFrameFillColor(0);
     gStyle->SetFrameBorderMode(0);
     gStyle->SetPadBorderMode(0);
     gStyle->SetPadColor(0);
     gStyle->SetOptStat("");
  
     // Title of the histogramm
     std::ostringstream titlestring;  
     titlestring <<  file_basename << ";Frequency [Hz] ("<< dfreq_multi <<" [Hz/bin]);Time [s] (" << delta_t << " [s/bin]);Intensity [a.u.]";
     TString title=titlestring.str();

     TApplication app("App", &argc, argv);
     TCanvas *canva_2dh = new  TCanvas("Time Resolved spectrum","Time Resolved spectrum",1000,800);
     canva_2dh->ToggleEditor();
     canva_2dh->SetCrosshair();
     canva_2dh->ToggleEventStatus();
     canva_2dh->ToggleToolBar();
     RAW->SetTitle(title);
     RAW->SetStats(kFALSE);
     RAW->Draw("zcol");
     RAW->GetXaxis()->SetLabelSize(0.025);
     RAW->GetXaxis()->SetTitleSize(0.025);
     RAW->GetXaxis()->SetTitleOffset(1.5);
     RAW->GetYaxis()->SetLabelSize(0.025);
     RAW->GetYaxis()->SetTitleSize(0.025);
     RAW->GetYaxis()->SetTitleOffset(1.5);
     RAW->GetZaxis()->SetLabelSize(0.025);
     RAW->GetZaxis()->SetTitleSize(0.025);
     RAW->GetZaxis()->SetTitleOffset(1.3);
      
     canva_2dh->Update(); // this line updates the canvas automatically, should come after Draw()

     if (draw){
	  canva_2dh->SaveAs(Form("%s.png", out_filename.c_str()));
     }
     else {
	  // The following line to connect the close button of the window manager to the main frame, in order to close properly.
	  ((TRootCanvas *)canva_2dh->GetCanvasImp())->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");

	  app.Run();
     }
     
     for (i_multi=0; i_multi < n_of_multi; i_multi++) delete [] multi_taper_spektrum[i_multi];
     delete [] multi_taper_spektrum;
     delete [] iq_data_p;
     delete RAW;


     return 0;
}

