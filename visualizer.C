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


//////////////////////////////////////
// IQT TIQ Visualizer Macro         //
// Works also as compiled version   //
// Use Makefile                     //
//////////////////////////////////////


// includes for the compiled version

#include <cstdlib>
#include <cstring>
#include <TH1D.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TRootCanvas.h>
#include <TApplication.h>
#include <TGFileDialog.h>
#include <TStyle.h>
#include <TColor.h>
#include <THStack.h>
#include <Riostream.h>
#include "header.h"


//______________________________________________________________________________
void draw_on_screen(TH1 * h, double begin, double end)
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

    int nx, ny1, ny2, i, j;
    if (h->GetDimension() == 1) {
        nx = h->GetNbinsX();
        TH1F* hcopy = new TH1F(h->GetName(), h->GetTitle(),
                nx, h->GetXaxis()->GetBinLowEdge(1), h->GetXaxis()->GetBinUpEdge(nx));
        hcopy->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
        for (i = 1; i <= nx; i++)
            hcopy->SetBinContent(i, (Float_t) h->GetBinContent(i));
        hcopy->Draw();
    } else {
        nx = h->GetNbinsX();
        double w = h->GetYaxis()->GetBinWidth(1);
        if (begin < 0) {
            cout << "warning: the begin time is underflow" << endl;
            ny1 = 0;
        } else {
            ny1 = (int) (begin / w);
        }
        if (end > h->GetYaxis()->GetBinLowEdge(h->GetNbinsY())) {
            cout << "warning: the end time is overflow" << endl;
            ny2 = h->GetNbinsY();
        } else {
            ny2 = (int) (end / w);
        }
        TH2F* hcopy = new TH2F(h->GetName(), h->GetTitle(),
                nx, h->GetXaxis()->GetBinLowEdge(1), h->GetXaxis()->GetBinUpEdge(nx),
                ny2 - ny1, ny1 * w, ny2 * w);
        hcopy->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
        hcopy->GetYaxis()->SetTitle(h->GetYaxis()->GetTitle());
        for (i = 1; i <= nx; i++)
            for (j = 1; j <= ny2 - ny1 + 1; j++)
                hcopy->SetBinContent(i, j, (Float_t) h->GetBinContent(i, j+ny1));
        hcopy->Draw("colz");
    }

    //h->Draw("colz");

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
TH1* do_read_from_file(const char* filename, const char* histo_name) // TH1 is the mother of all histos
{
    TH1* h = NULL;
    TFile f(filename);
    string dir(histo_name, 5);

    if (!gFile->GetDirectory(dir.c_str())) {
        cout << "The root file does not contain a histogram with the name of `"
            << histo_name << "'" << endl;
        cout << "Try one of the following:\n" << endl;
        TListIter* iter = (TListIter*) gFile->GetListOfKeys()->MakeIterator();
        while (TDirectoryFile* dirf = (TDirectoryFile*) iter->Next()) {
            gDirectory = gFile;
            gDirectory->cd(dirf->GetName());
            gDirectory->ls();
        }
        f.Close();
        return h;
    }

    gDirectory->cd(dir.c_str());
    if (!gDirectory->GetKey(histo_name)) {
        cout << "The root file does not contain a histogram with the name of `"
            << histo_name << "'" << endl;
        cout << "Try one of the following:\n" << endl;
        gDirectory->ls();
        f.Close();
        return h;
    }

    if (!strcmp(dir.c_str(), "Oscil")) {
        h = (TH1*)gDirectory->Get(histo_name);
        h->SetDirectory(0);	// magic! histo should not belong to any directory in memory
    } else {
        Header* header = (Header*) gDirectory->Get(Form("%s_header", dir.c_str()));
        header->Show();
        cout << "Data from instrument " << dir << endl;
        h = (TH1*)gDirectory->Get(histo_name);
        h->SetDirectory(0);
    }

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
bool draw_stack(const char* filename) {
    THStack* hsi = new THStack("Oscil_stack_inj", "Kicker Signals at Injection;Time [s];Voltage [V]");
    THStack* hse = new THStack("Oscil_stack_ext", "Kicker Signals at Extraction;Time [s];Voltage [V]");
    TH1D* hi = NULL;
    TH1D* he = NULL;

    int color = 2;
    for (int i = 1; i <= 4; i++) {
        hi = (TH1D*)do_read_from_file(filename, Form("Oscil_C%d_inj", i));
        he = (TH1D*)do_read_from_file(filename, Form("Oscil_C%d_ext", i));
        if (!(hi && he))
            return false;
        hi->SetLineColor(color);
        he->SetLineColor(color);
        hsi->Add(hi);
        hse->Add(he);
        color = color + 2; // avoid ugly green
    }

    TCanvas* c = new TCanvas("c", "kicker", 1280, 720);
    c->ToggleEditor();
    c->ToggleEventStatus();
    c->ToggleToolBar();
    c->Divide(2, 1);

    c->cd(1);
    gPad->SetGrid();
    hsi->Draw("nostack");
    //hsi->Draw("pads");
    c->cd(2);
    gPad->SetGrid();
    hse->Draw("nostack");
    //hse->Draw("pads");

    // The following line to connect the close button of the window manager to the main frame, in order to close properly.
    ((TRootCanvas *)c->GetCanvasImp())->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
    return true;
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
int main(int argc, char **argv) {
    if (argc != 3 && argc != 5){
        cout << "Usage:\n\n";
        cout << "    visualizer histo_name file_name\n\n";
        cout << "or\n\n";
        cout << "    visualizer histo_name file_name begin_time end_time\n\n";
        cout << "or\n\n";
        cout << "    visualizer stack file_name\n\n";
        return 1;
    }

    char * which_histo;
    char * filename;
    double begin = 1;
    double end = 0;

    which_histo = argv[1];
    filename = argv[2];
    if (argc == 5) {
        begin = atof(argv[3]);
        end = atof(argv[4]);
    }

    if (!exists(argv[2])) {cout << "No such file " << argv[2] << "\nAborting ... \n" << endl; return 2;}

    TApplication theApp("theApp", &argc, argv);

    if(!strcmp(which_histo, "stack")) {
        if (!draw_stack(filename))
            return 3;
    } else {
        TH1* h = do_read_from_file(filename, which_histo);
        if(h == NULL)
            return 4;

        cout << "Showing plot `" << which_histo << "' from file `" << filename << "'" << endl;
        if (begin > end)
            draw_on_screen(h, 0, h->GetYaxis()->GetBinLowEdge(h->GetNbinsY()));
        else
            draw_on_screen(h, begin, end);
    }

    theApp.Run();

    return 0;
}
