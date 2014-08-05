//
// (c) Copyright:
// X. Chen 2013
//
// this is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// time2root is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
//
// along with this software. If not, see <http://www.gnu.org/licenses/>.


/*------------------------------------------------------------------*
 | This program is aimed to provide multitaper method for spectral  |
 | density estimation. This method is developed by David J. Thomson |
 | in 1982, while he was working at the Bell Laboratories. The      |
 | algorithm is taken from one great book about spectral analysis:  |
 | "Spectral Analysis for Physical Applications: Multitaper and     |
 | Conventional Univariate Techniques, by Donald B. Percival and    |
 | Andrew T. Walden", to which one can refer for more information.  |
 *------------------------------------------------------------------*/


#include <TMath.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <fftw3.h>


#define cNW     4
#define cKMax   (2 * cNW - 2)
#define cFrmPt  1024
#define SQ(x)   ((x)*(x))

/* energy concentration weighted taper sequences */
static double Taper[cKMax][cFrmPt] = {};
/* generate DPSSs and corresponding concentrations */
void Eigen(void);

/* Seq is a time sequence whose spectral density shall be estimated, and
   PSD stores the estimated power spectral density. */
bool Multitaper(fftw_complex Seq[cFrmPt], double PSD[cFrmPt]) {
    if (!Taper[0][0]) { /* tapers not loaded to memory */
        fprintf(stdout, "taper sequences not found in memory.\n");
        FILE* fp = NULL;
        char FileName[32] = {};
        sprintf(FileName, "km%d_n%d_nw%d.egn", cKMax, cFrmPt, cNW);

        /* generate DPSSs and concentrations, then write tapers to disk */
        if (!(fp = fopen(FileName, "rb"))) {
            fprintf(stdout, "Generating DPSS right now.\n");
            Eigen();

            fprintf(stdout, "Saving taper sequences to disk.\n");
            fp = fopen(FileName, "wb");
            fwrite(Taper, sizeof(double), cKMax * cFrmPt, fp);
            fclose(fp);
        } else { /* read DPSSs and concentrations from disk */
            fprintf(stdout, "Loading taper sequences from disk.\n");
            fread(Taper, sizeof(double), cKMax * cFrmPt, fp);
            fclose(fp);
        }
    }

    int IdxVec, IdxElem;
    fftw_complex* SeqClone = (fftw_complex*) fftw_malloc(
            sizeof(fftw_complex) * cFrmPt);
    fftw_plan p = fftw_plan_dft_1d(cFrmPt, SeqClone, SeqClone,
            FFTW_FORWARD, FFTW_ESTIMATE);
    for (IdxElem = 0; IdxElem < cFrmPt; IdxElem++)
        PSD[IdxElem] = .0;

    for (IdxVec = 0; IdxVec < cKMax; IdxVec++) {
        for (IdxElem = 0; IdxElem < cFrmPt; IdxElem++) { /* tapering */
            SeqClone[IdxElem][0] = Seq[IdxElem][0] * Taper[IdxVec][IdxElem];
            SeqClone[IdxElem][1] = Seq[IdxElem][1] * Taper[IdxVec][IdxElem];
        }
        fftw_execute(p); /* discrete fourier transform */
        for (IdxElem = 0; IdxElem < cFrmPt; IdxElem++)
            PSD[IdxElem] += SQ(SeqClone[IdxElem][0]) + SQ(SeqClone[IdxElem][1]);
    }

    fftw_destroy_plan(p);
    fftw_free(SeqClone);
    return true;
}

void Eigen(void) {
    int i, j;
    double w = (double) cNW / cFrmPt;
    TMatrixDSym A(cFrmPt);
    TMatrixD Vec(cFrmPt, cFrmPt);
    TVectorD Val(cFrmPt);

    for (i = 0; i < cFrmPt; i++) {
        A(i, i) = 2 * w;
        for (j = i+1; j < cFrmPt; j++)
            A(i, j) = A(j, i) = TMath::Sin(TMath::TwoPi() * w * (j-i)) /
                (TMath::Pi() * (j-i));
    }
    Vec = A.EigenVectors(Val);

    w = .0;
    for (j = 0; j < cKMax; j++)
        w += Val(j);
    for (j = 0; j < cKMax; j++)
        for (i = 0; i < cFrmPt; i++)
            Taper[j][i] = Vec(i, j) * TMath::Sqrt(Val(j) / w);
}
