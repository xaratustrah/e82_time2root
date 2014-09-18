//
// (c) Copyright:
// X. Chen, M. S. Sanjari and P. Buehler 2014
// N. Winckler and M. S. Sanjari 2010 - 2012
// F. Nolden 2009
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


////////////////////////////////////////////////////////////////////////////////
// This program time domain data based on Tektronix IQT and TIQ               //
// spectrum estimation and stores the results in a root file.                 //
////////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <iomanip>
#include <string>
#include <cstdio>
#include <cstdlib> // for exit
#include <cstdarg>
#include <gsl/gsl_fft_complex.h>

#include <TH2D.h>
#include <THStack.h>
#include <TFile.h>
#include <TComplex.h>
#include <TDatime.h>

#include "iqtdata.h"
#include "FritzDPSS.h"
#include "multitaper.h"
#include "header.h"
#include "tiq.h"

using namespace std;


//______________________________________________________________________________
void do_append_to_file(const char* filename, TNamed* obj)
{
    TFile f(filename, "update");
    string dir(obj->GetName(), 5);

    if (!gFile->GetDirectory(dir.c_str()))
        gFile->mkdir(dir.c_str());
    gDirectory->cd(dir.c_str());
    if (gDirectory->GetKey(obj->GetName())) {
        cout <<"An object with the same name `" << obj->GetName() <<
            "' already exists in the file `" << filename << "'. No overwrite!" << endl;
        return;
    }

    if (gDirectory->WriteTObject(obj)) {
        cout << "Object `" << obj->GetName() << "' successfully written to file `" << filename << "'" << endl;
        f.Close();
        return;
    }
    else {
        cerr << "Object `" << obj->GetName() << "' failed to be written to file `" << filename << "'" << endl;
        f.Close();
        exit(EXIT_FAILURE);
    }
}

//______________________________________________________________________________
void do_process_iqt(const char* outfile, const char* infile, const char* basename,
        Int_t ntap, bool* towrite)
{
    IQData__TekRSA3303B my_iq_data;

    cout << "Processing IQT file: " << infile << endl;

    // check the file
    if (!my_iq_data.ReadFile(infile))
    {
        cout << "Error reading file." << endl;
        exit(EXIT_FAILURE);
    }

    Double_t delta_t;

    if(!my_iq_data.GetDeltaT(delta_t)) 
    {
        cout << "Error during the check of time interval" << endl;
        cout << "Program stopped." << endl;
        exit(EXIT_FAILURE);
    }

    Header* header = new Header("RSA30_header", basename);
    header->SetValidFrames(my_iq_data.ValidFrames());
    header->SetFrameLength(my_iq_data.FrameLength());
    header->SetCenterFrequency(my_iq_data.CenterFrequency());
    header->SetSpan(my_iq_data.Span());
    header->SetGainOffset(my_iq_data.GainOffset());
    header->SetDateTime(my_iq_data.DateTime());

    // Fill out the constants

    Int_t valid_frames = my_iq_data.ValidFrames();
    Double_t center_frequency = my_iq_data.CenterFrequency(); 
    Int_t block_size = my_iq_data.BlockSize();

    if (valid_frames == block_size) 
    {
        cout << "All possible " << block_size
            << " Frames with each 1024 Points are valid." << endl;
    }
    else 
    {
        cout << "From the " << block_size
            << "  Frames with each 1024 Points, only "
            << valid_frames << " are valid." << endl;
        cout << "Program stopped." << endl;
        exit(EXIT_FAILURE);
    }

    cout << "There are " << block_size*1024 << " measured points in total." << endl;

    // FFT parameters

    Int_t Npoint = 1024; // Number of points for making one FFT frame
    Int_t first_index = 0;
    complex <Double_t> * iq_data_p = new complex <Double_t> [Npoint];
    fftw_complex* sgn = new fftw_complex[Npoint];
    fftw_plan p = fftw_plan_dft_1d(Npoint, sgn, sgn, FFTW_FORWARD, FFTW_ESTIMATE);


    /// for Multitaper, delta_t cover on RSA always 1024 Points

    Double_t dtm = (Double_t) Npoint / 1024. * delta_t; // Multitaper time frame
    Double_t dfreq_multi =  1. / dtm; // Nyquist!

    cout << "1 frame in frequency=" << dfreq_multi << " Hz" << endl;

    Double_t start_frequency = center_frequency - Npoint * dfreq_multi / 2.;
    Int_t ifmin_multi = (Npoint * 3) / 16;
    Int_t ifmax_multi = 1 + (Npoint * 13) / 16;
    Int_t nw = (ntap + 2) / 2;
    FritzDPSS my_dpss(Npoint, nw);
    my_dpss.GetFromDisk();

    Int_t n_of_multi = block_size;
    Int_t i_multi;
    Double_t ** multi_taper_spektrum = new Double_t * [n_of_multi]; 
    Double_t ** fft = new Double_t * [n_of_multi];
    for (i_multi=0; i_multi < n_of_multi; i_multi++) {
        multi_taper_spektrum[i_multi] = new Double_t [Npoint];
        fft[i_multi] = new Double_t [Npoint];
    }

    i_multi=0;

    // define time histograms

    TH1D * hmag = new TH1D("RSA30_mag", Form("%s;Time [s];Magnitude [a.u.]", basename), 1024 * n_of_multi, 0, n_of_multi * dtm);
    TH1D * hphs = new TH1D("RSA30_phs", Form("%s;Time [s];Phase [a.u.]", basename), 1024 * n_of_multi, 0, n_of_multi * dtm);

    for (first_index=0; first_index <= 1024 * (n_of_multi-1); first_index += 1024)
    {
        Int_t index = first_index;
        complex <Double_t> iq_value;

        for (Int_t ip=0; ip < Npoint; ip++) 
        {
            if (!my_iq_data.GetIQ(iq_value, index))
            {
                cout << "Error retrieving the value at index " << index << endl;
                exit(EXIT_FAILURE);

            }
            index++;
            iq_data_p[ip] = iq_value;
            sgn[ip][0] = real(iq_value);
            sgn[ip][1] = imag(iq_value);

            // Fill out time histograms

            hmag->SetBinContent(index + ip, abs(iq_value));
            hphs->SetBinContent(index + ip, arg(iq_value));

        }

        // Spectrum calculation

        my_dpss.GetSpectrum(iq_data_p, multi_taper_spektrum[i_multi]);
        fftw_execute(p);
        for (int index = 0; index < Npoint; index++)
            fft[i_multi][index] = SQ(sgn[index][0]) + SQ(sgn[index][1]);

        // progress bar
        if (!(i_multi % 10))
            cout << "processing... " << fixed << setw(5) << setprecision(2) << right << (double)i_multi/n_of_multi*100 << "%\r" << flush;

        i_multi++;
    }

    cout << endl;

    Int_t binf=ifmax_multi-ifmin_multi;
    Double_t fstart=start_frequency+(Double_t)ifmin_multi*dfreq_multi;
    Double_t fend= start_frequency+(Double_t)(ifmax_multi-1)*dfreq_multi;

    Int_t bint = n_of_multi;

    // Create the 2D histogram

    // add number_of_tapers to the histogram name

    TH2D* hmtpsd = new TH2D(Form("RSA30_mtpsd%d", ntap),Form("%s;Frequency[Hz] (%g [Hz/bin]);Time [s] (%g [s/bin]);Intensity [a.u.]", basename, dfreq_multi, delta_t),binf,fstart,fend,bint,0.0,(Double_t)(bint)*delta_t);
    TH2D* hfft = new TH2D("RSA30_fft",Form("%s;Frequency[Hz] (%g [Hz/bin]);Time [s] (%g [s/bin]);Intensity [a.u.]", basename, dfreq_multi, delta_t),binf,fstart,fend,bint,0.0,(Double_t)(bint)*delta_t);

    // constant of multiplication for z axis
    Double_t constant = 1.0;
    for (Int_t ifreq = ifmin_multi; ifreq < ifmax_multi; ifreq++)
    {
        for (i_multi = 0; i_multi < n_of_multi; i_multi++)
        {
            hmtpsd->SetBinContent(ifreq-ifmin_multi+1,i_multi+1, multi_taper_spektrum[i_multi][ifreq] * constant);
            hfft->SetBinContent(ifreq-ifmin_multi+1,i_multi+1, fft[i_multi][(Npoint-ifmax_multi/2-ifmin_multi/2+ifreq)%Npoint] * constant);
        } 
    }


    // Store the spectra

    do_append_to_file(outfile, header); // header is written in any case
    if (towrite[0]) do_append_to_file(outfile, hmtpsd);
    if (towrite[1]) do_append_to_file(outfile, hfft);
    if (towrite[2]) do_append_to_file(outfile, hmag);
    if (towrite[3]) do_append_to_file(outfile, hphs);

    // clean up

    delete header;
    delete hmtpsd;
    delete hfft;
    delete hmag;
    delete hphs;
    fftw_destroy_plan(p);
    for (i_multi=0; i_multi < n_of_multi; i_multi++) {
        delete [] multi_taper_spektrum[i_multi];
        delete [] fft[i_multi];
    }
    delete [] multi_taper_spektrum;
    delete [] fft;
    delete [] iq_data_p;
    delete [] sgn;

}

//______________________________________________________________________________
bool do_process_tiq(const char* outfile, FILE* fp, Info_t* pInfo, int ntap, bool* towrite) {
    fseek(fp, pInfo->Offset, SEEK_SET);
    int bins = pInfo->Span * pInfo->Intvl * cFrmPt;
    int blksz = pInfo->NumPt / cFrmPt;
    pInfo->NumPt = blksz * cFrmPt;
    cout << "number of points: " << pInfo->NumPt << endl;
    cout << "center frequency: " << pInfo->CenFreq << " Hz" << endl;
    cout << "acquisition bandwidth: " << pInfo->Span << " Hz" << endl;

    Header* header = new Header("", pInfo->File);
    header->SetValidFrames(blksz); 
    header->SetFrameLength(cFrmPt * pInfo->Intvl);
    header->SetCenterFrequency(pInfo->CenFreq);
    header->SetSpan(pInfo->Span);
    header->SetScaling(pInfo->Scaling);
    header->SetDateTime(pInfo->DaTm.AsSQLString());
    header->SetSerialNumber(pInfo->ID);

    string prefix;
    if (!strcmp(header->GetSerialNumber(), "B010426"))
        prefix = "RSA51";
    else
        prefix = "RSA52";
    header->SetName(Form("%s_header", prefix.c_str()));

    TH1D* hmag = new TH1D(Form("%s_mag", prefix.c_str()), pInfo->File, pInfo->NumPt,
            0, pInfo->NumPt * pInfo->Intvl);
    TH1D* hphs = new TH1D(Form("%s_phs", prefix.c_str()), pInfo->File, pInfo->NumPt,
            0, pInfo->NumPt * pInfo->Intvl);
    TH2D* hmtpsd = new TH2D(Form("%s_mtpsd%d", prefix.c_str(), ntap),
            Form("%s;frequency[Hz] (%g[Hz/bin]);time[s] (%g[s/bin]);"
                "intensity[a.u.]", pInfo->File, 1. / pInfo->Intvl / cFrmPt,
                pInfo->Intvl * cFrmPt), bins, pInfo->CenFreq - pInfo->Span / 2,
            pInfo->CenFreq + pInfo->Span / 2, blksz, 0,
            pInfo->NumPt * pInfo->Intvl);
    TH2D* hfft = new TH2D(Form("%s_fft", prefix.c_str()),
            Form("%s;frequency[Hz] (%g[Hz/bin]);time[s] (%g[s/bin]);"
                "intensity[a.u.]", pInfo->File, 1. / pInfo->Intvl / cFrmPt,
                pInfo->Intvl * cFrmPt), bins, pInfo->CenFreq - pInfo->Span / 2,
            pInfo->CenFreq + pInfo->Span / 2, blksz, 0,
            pInfo->NumPt * pInfo->Intvl);
    //h->GetYaxis()->SetTimeDisplay(true);
    //h->GetYaxis()->SetTimeFormat("#splitline{%H:%M:%S}{%d.%m.%y}");
    //h->GetYaxis()->SetTimeOffset(dt->Convert());

    int j, k, tmp;
    Multitaper multitaper(cFrmPt, ntap, ntap/2+1);
    double* mtpsd = (double*) malloc(sizeof(double) * cFrmPt);
    fftw_complex* sgn = (fftw_complex*) fftw_malloc(
            sizeof(fftw_complex) * cFrmPt);
    fftw_plan p = fftw_plan_dft_1d(cFrmPt, sgn, sgn,
            FFTW_FORWARD, FFTW_ESTIMATE);

    for (j = 0; j < blksz; j++) {
        for (k = 0; k < cFrmPt; k++) {
            fread(&tmp, 4, 1, fp);
            sgn[k][0] = pInfo->Scaling * tmp;
            fread(&tmp, 4, 1, fp);
            sgn[k][1] = pInfo->Scaling * tmp;
            hmag->SetBinContent(j*cFrmPt + k + 1,
                    sqrt(SQ(sgn[k][0]) + SQ(sgn[k][1])));
            hphs->SetBinContent(j*cFrmPt + k + 1,
                    atan2(sgn[k][1], sgn[k][0]));
        }
        multitaper.estimate(sgn, mtpsd);
        fftw_execute(p);
        for (k = 0; k < bins; k++) {
            tmp = (cFrmPt - bins/2 + k) % cFrmPt;
            hfft->SetBinContent(k+1, j+1, SQ(sgn[tmp][0])+SQ(sgn[tmp][1]));
            hmtpsd->SetBinContent(k+1, j+1, mtpsd[tmp]);
        }
        if (!(j % 10))
            cout << "processing... " << fixed << setw(5) << setprecision(2)
               << right << (double)j/blksz*100 << "%\r" << flush;
    }
    cout << endl;

    do_append_to_file(outfile, header); // header is written in any case
    if (towrite[0]) do_append_to_file(outfile, hmtpsd);
    if (towrite[1]) do_append_to_file(outfile, hfft);
    if (towrite[2]) do_append_to_file(outfile, hmag);
    if (towrite[3]) do_append_to_file(outfile, hphs);

    delete header;
    delete hmag;
    delete hphs;
    delete hfft;
    delete hmtpsd;
    fftw_destroy_plan(p);
    fftw_free(sgn);
    free(mtpsd);
    return true;
}


//______________________________________________________________________________
bool prepare_tiq(const char * outfile, const char* infile, const char* basename,
        Int_t ntap, bool* towrite)
{
    Info_t Info = {basename, {}, 0, 0, .0, .0, .0, .0, TDatime(2001, 1, 1, 0, 0, 0)};
    FILE *fp = fopen(infile, "rb");
    if (!fp) {
        fprintf(stderr, "Error: can't open file `%s'!\n", infile);
    }

    /* file header */
    if (!SetInfo(fp, &Info)) {
        return false;
    }

    /* visualization */
    if (!do_process_tiq(outfile, fp, &Info, ntap, towrite)) {
        return false;
    }

    fclose(fp);

    return true;
}


//______________________________________________________________________________
// Make a histogram
TH1D* make_histo_csv_cpp_style(const char* filename, const char* basename) {
    ifstream data;
    data.open(filename);
    string tmp;
    getline(data, tmp); // skip the line
    getline(data, tmp); // skip the line
    getline(data, tmp); // skip the line
    getline(data, tmp); // skip the line
    getline(data, tmp); // skip the line

    // containers
    std::vector<Double_t> xvals, yvals;  
    Double_t xval, yval;
    char comma;
    for (data >> xval >> comma >> yval; data.good(); data >> xval >> comma >> yval) {
        xvals.push_back(xval);
        yvals.push_back(yval);
    }
    data.close();

    tmp = string(basename, 6);
    TH1D* h = new TH1D(Form("Oscil_%s", tmp.c_str()), basename, xvals.size(), xvals.front(), xvals.back());
    for (int i = 0; i < xvals.size(); i++)
        h->SetBinContent(i+1, yvals.at(i));	  

    return h;
}


//______________________________________________________________________________
// Make a histogram
TH1F * make_histo_csv_c_style(const char* filename){

    FILE * fp;
    fp=fopen(filename, "r");

    // skipping the header
    char s;
    int cnt = 0;
    while(cnt < 5){
        s = fgetc(fp);
        if (s == '\n') ++cnt;
    }

    // containers
    std::vector<Double_t> xvals, yvals;  
    Double_t xval, yval;

    while((s = fgetc(fp)) != EOF)
    {
        fscanf(fp, "%lf,%lf", &xval, &yval);
        xvals.push_back(xval);
        yvals.push_back(yval);
    }

    TH1F * h = new TH1F (Form("h_csv_time_%c%c", filename[0], filename[1]), "Kicker Plot", xvals.size(), xvals.front(), xvals.back());
    for(Int_t i=0; i < xvals.size(); ++i)
    {
        h->SetBinContent(i+1, yvals.at(i));	  
    }
    return h;
}

//______________________________________________________________________________
inline bool exists (const char* file) {
    ifstream f(file);
    if (f.good()) {
        f.close();
        return true;
    } else {
        f.close();
        return false;
    }   
}

//______________________________________________________________________________
void usage () {

  cout << "Usage:\n\n";
  cout << "    time2root root_filename time_filename\n\n";
  cout << "or\n\n";
  cout << "    time2root root_filename time_filename towrite\n\n";
  cout << "or\n\n";
  cout << "    time2root root_filename time_filename towrite number_of_tapers\n\n";
  cout << "    towrite: string of 4 digits - 0 or 1 - specifies type of outputs\n";
  cout << "        [0]: multi-taper estimation\n";
  cout << "        [1]: plain FFT             \n";
  cout << "        [2]: magnitude versus time \n";
  cout << "        [3]: phase versus time   \n\n";
  cout << "default value\n\n";
  cout << "        towrite = 1000, number_of_tapers = 6\n\n";
  
  exit(EXIT_FAILURE);

}

//______________________________________________________________________________
int main(int argc, char** argv)
{

    // towrite:
    // [0]: multi-taper estimation
    // [1]: plain FFT
    // [2]: magnitude versus time
    // [3]: phase versus time
    bool towrite[4] = {true, false, false, false};

    // by default 6 tapers
    Int_t ntap = 6;

    if (argc < 3 || argc > 5) usage();

    if (!exists(argv[2])) {
        cerr << "No such file `" << argv[2] << "'\nAborting ... \n" << endl;
        exit(EXIT_FAILURE);
    }

    /* extract the file basename and extension */
    const char* basename = argv[2];
    for (const char* p = argv[2]; *p; p++)
        if (*p == '/')
            basename = p + 1;
    const char* extension = basename;
    for (const char* p = basename; *p; p++)
        if (*p == '.')
            extension = p + 1;

    // parse towrite, if exists
    if (argc >= 4)
        for (int i = 0; i < 4; i++)
            towrite[i] = (bool) (argv[3][i] - '0');

    // parse number_of_tapers, if exists
    if (argc == 5)
        ntap = atoi(argv[4]);

    // process file
    if (!strcmp(extension, "iqt") || !strcmp(extension, "IQT")) // iqt file
        do_process_iqt(argv[1], argv[2], basename, ntap, towrite);
    else if (!strcmp(extension, "tiq") || !strcmp(extension, "TIQ")) // tiq file
        prepare_tiq(argv[1], argv[2], basename, ntap, towrite);
    else if (!strcmp(extension, "csv") || !strcmp(extension, "CSV")) // csv file
        do_append_to_file(argv[1], make_histo_csv_cpp_style(argv[2], basename));
    else {
        cerr << "Supported file types are:\n\n" << "    iqt, tiq and csv\n\n";
        exit(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;

}
