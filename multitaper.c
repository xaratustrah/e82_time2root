//
// (c) Copyright:
// X. Chen 2013
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
#include <TDatime.h>

#include "tiq.h"

#define cNW     4
#define cKMax   (2 * cNW - 1)

/* Discrete Prolate Spheroidal Sequences */
static double EigenVec[cKMax][cFrmPt] = {};
/* energy concentrations for DPSSs */
static double EigenVal[cKMax] = {};
/* generate DPSSs and corresponding concentrations */
void Eigen(void);

/* Seq is a time sequence whose spectral density shall be estimated, and
   PSD stores the estimated power spectral density. */
bool Multitaper(complex<double> Seq[cFrmPt], double PSD[cFrmPt])
{
    FILE* fp = NULL;
    char FileName[32] = {};
    sprintf(FileName, "km%d_n%d_nw%d.egn", cKMax, cFrmPt, cNW);

    if (!EigenVal[0]) /* DPSSs not loaded to memory */
    {
        fprintf(stdout, "DPSS not found in memory.\n");

        /* generate DPSSs and concentrations, then write to disk */
        if (!(fp = fopen(FileName, "rb")))
        {
            fprintf(stdout, "Generating DPSS right now.\n");
            Eigen();

            fprintf(stdout, "Saving DPSS to disk.\n");
            fp = fopen(FileName, "wb");
            fwrite(EigenVec, sizeof(double), cKMax * cFrmPt, fp);
            fwrite(EigenVal, sizeof(double), cKMax, fp);
            fclose(fp);
        }
        else /* read DPSSs and concentrations from disk */
        {
            fprintf(stdout, "Loading DPSS from disk.\n");
            fread(EigenVec, sizeof(double), cKMax * cFrmPt, fp);
            fread(EigenVal, sizeof(double), cKMax, fp);
            fclose(fp);
        }
    }

    int IdxVec, IdxElem;
    double Norm = .0;
    complex<double> SeqClone[cFrmPt] = {};
    fftw_plan p = fftw_plan_dft_1d(cFrmPt,
            reinterpret_cast<fftw_complex*> (SeqClone),
            reinterpret_cast<fftw_complex*> (SeqClone),
            FFTW_FORWARD, FFTW_MEASURE);

    for (IdxVec = 0; IdxVec < cKMax; IdxVec++)
    {
        Norm += EigenVal[IdxVec];

        for (IdxElem = 0; IdxElem < cFrmPt; IdxElem++) /* window the sequence */
            SeqClone[IdxElem] = Seq[IdxElem] * EigenVec[IdxVec][IdxElem];

        fftw_execute(p); /* discrete fourier transform */

        if (IdxVec) /* eigenvalue weighting scheme */
            for (IdxElem = 0; IdxElem < cFrmPt; IdxElem++)
                PSD[IdxElem] += norm(SeqClone[IdxElem]) * EigenVal[IdxVec];
        else
            for (IdxElem = 0; IdxElem < cFrmPt; IdxElem++)
                PSD[IdxElem] = norm(SeqClone[IdxElem]) * EigenVal[IdxVec];
    }
    for (IdxElem = 0; IdxElem < cFrmPt; IdxElem++)
        PSD[IdxElem] /= Norm;

    fftw_destroy_plan(p);
    return true;
}

void Eigen(void)
{
    int i, j;
    double w = (double) cNW / cFrmPt;
    TMatrixDSym A(cFrmPt);
    TMatrixD Vec(cFrmPt, cFrmPt);
    TVectorD Val(cFrmPt);

    for (i = 0; i < cFrmPt; i++)
    {
        A(i, i) = 2 * w;
        for (j = i+1; j < cFrmPt; j++)
            A(i, j) = A(j, i) = TMath::Sin(TMath::TwoPi() * w * (j-i)) /
                                (TMath::Pi() * (j-i));
    }
    Vec = A.EigenVectors(Val);
    for (j = 0; j < cKMax; j++)
    {
        EigenVal[j] = Val(j);
        for (i = 0; i < cFrmPt; i++)
            EigenVec[j][i] = Vec(i, j);
    }
}
