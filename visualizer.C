//
// (c) Copyright:
// X. Chen and M.S. Sanjari 2014
//
// this is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// barion is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
//
// along with this software. If not, see <http://www.gnu.org/licenses/>.


////////////////////////////////////
// IQT TIQ Visualizer Macro	  //
// Works also as compiled version //
// Use Makefile			  //
////////////////////////////////////


// includes for the compiled version

#include <TH1D.h>
#include <TH1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TRootCanvas.h>
#include <TApplication.h>
#include <TGFileDialog.h>
#include <TStyle.h>
#include <TColor.h>
#include <THStack.h>
#include <Riostream.h>


//______________________________________________________________________________
void draw_on_screen(TH1 * h)
{
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

     TCanvas * c = new TCanvas("c", "Measurements", 1600, 700);
     c->SetGrid();
     c->ToggleEditor();
     c->ToggleEventStatus();
     c->ToggleToolBar();

     h->Draw("colz");

     c->Modified();
     c->Update(); // this line updates the canvas automatically, should come after Draw()
     // The following line to connect the close button of the window manager to the main frame, in order to close properly.
     ((TRootCanvas *)c->GetCanvasImp())->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");

     return;
}

//______________________________________________________________________________
void compare(TH1* h1, TH1* h2) {
     TCanvas* c = new TCanvas("c", "comparison", 1280, 720);
     c->ToggleEditor();
     c->ToggleEventStatus();
     c->ToggleToolBar();
     c->SetGrid();

     c->Divide(2, 1);
     c->cd(1);
     h1->Draw("colz");
     c->cd(2);
     h2->Draw("colz");
}

//______________________________________________________________________________
TH1 * do_read_from_file(const char* filename, const char * histo_name) // TH1 is the mother of all histos
{
     TFile f (filename, "read");
     if(!f.GetListOfKeys()->Contains(histo_name))
     {
	  cout << "The binary file does not contain a histogram with the name of " << histo_name << "." << endl;
	  cout << "Try one of the following:\n\n";
	  f.GetListOfKeys()->ls();
	  return 0;
     }
     TH1 * h = (TH1*)f.Get(histo_name);
     h->SetDirectory(0);	// magic! histo should not belong to any directory in memory
     f.Close();
     return h;
}

//______________________________________________________________________________
//
char * file_name_gui (){

     TGFileInfo fFileInfo;
     const char *types[]={"ROOT files","*.root",0,0};
     fFileInfo.fFileTypes=types;
     TString parentDir(".");
     fFileInfo.fIniDir=StrDup(parentDir);
     new TGFileDialog (gClient->GetRoot(), 0, kFDOpen, &fFileInfo);

     return ((char*)fFileInfo.fFilename);
}


//______________________________________________________________________________
//
void draw_stack(const char * filename){
     THStack *hs = new THStack("h_kicker_stack","Kicker signals");
     TFile f (filename, "read");
     hs->SetTitle("Kicker signals;Voltage [v];Time [s]");
     TH1D * h;
     string histo_name;
     int color = 2;
     for (int i = 1; i <= 4; ++i)
     {
	  histo_name = Form("h_kicker_time_C%d", i);
	  if(f.GetListOfKeys()->Contains(histo_name.c_str()))
	  {
	       h = (TH1D*)do_read_from_file(filename, histo_name.c_str());
	       h->SetLineColor(color);
	       hs->Add(h);
	       color = color + 2; // avoid ugly green
	  }
     }

     TCanvas * c = new TCanvas("c", "Measurements", 1600, 700);
     c->SetGrid();
     c->ToggleEditor();
     c->ToggleEventStatus();
     c->ToggleToolBar();

     hs->Draw("nostack");
//     hs->Draw("pads");

     c->Modified();
     c->Update(); // this line updates the canvas automatically, should come after Draw()
     // The following line to connect the close button of the window manager to the main frame, in order to close properly.
     ((TRootCanvas *)c->GetCanvasImp())->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");

     return;

}

//______________________________________________________________________________
inline bool exists (const std::string& name) {
    ifstream f(name.c_str());
    if (f.good()) {
        f.close();
        return true;
    } else {
        f.close();
        return false;
    }   
}


//______________________________________________________________________________
//
int main(int argc, char **argv)
{
     if (argc != 3){
     	  cout << "Usage:\n\n";
      	  cout << "    visualizer histo_name file_name\n\n";
     	  cout << "or\n\n";
     	  cout << "    visualizer stack file_name\n\n";
     	  return 1;
     }

     char * which_histo;
     char * filename;

     which_histo = argv[1];
     filename = argv[2];

     if (!exists(argv[2])) {cout << "No such file " << argv[2] << "\nAborting ... \n" << endl; return 2;}

     TApplication theApp("theApp", &argc, argv);

     if(!strcmp(which_histo, "stack"))
	  draw_stack(filename);
     else{
     TH1 * h = do_read_from_file(filename, which_histo);
     if(h == 0) return 3;
     cout << "Showing plot " << which_histo << " from file " << filename << "." << endl;
     draw_on_screen(h);

     }


     theApp.Run();

     return 0;
}



