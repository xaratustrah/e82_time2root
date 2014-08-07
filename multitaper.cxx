//
// (c) Copyright:
// X. Chen 2014
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

#include "multitaper.h"

void Multitaper::build() {
    taper = new double* [n_taper];
    for (int i = 0; i < n_taper; i++) taper[i] = new double [n_sample];
    seqcopy = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n_sample);
    plan = fftw_plan_dft_1d(n_sample, seqcopy, seqcopy,
            FFTW_FORWARD, FFTW_ESTIMATE);
}

void Multitaper::load() {
    stringstream ss;
    ss << "n" << n_sample << "_k" << n_taper << "_nw" << nwidth << ".egn";
    string fname = ss.str();
    ifstream fin;
    fin.open(fname.c_str(), ios::binary);
    if (fin.is_open()) {
        for (int i = 0; i < n_taper; i++)
            fin.read((char*) taper[i], sizeof(double) * n_sample);
        cout << "multitaper functions read from disk" << endl;
        fin.close();
    } else {
        cout << "generating DPSS right now..." << endl;
        generate();
        ofstream fout;
        fout.open(fname.c_str(), ios::binary);
        for (int i = 0; i < n_taper; i++)
            fout.write((char*) taper[i], sizeof(double) * n_sample);
        cout << "multitaper functions written to disk" << endl;
        fout.close();
    }
}

void Multitaper::generate() {
    int i, j;
    double w = (double) nwidth / n_sample;
    TMatrixDSym A(n_sample);
    TMatrixD vec(n_sample, n_sample);
    TVectorD val(n_sample);
    for (i = 0; i < n_sample; i++) {
        A(i, i) = 2 * w;
        for (j = i + 1; j < n_sample; j++) {
            A(i, j) = A(j, i) = TMath::Sin(TMath::TwoPi() * w * (j-i)) /
                (TMath::Pi() * (j-i));
        }
    }
    vec = A.EigenVectors(val);

    w = .0;
    for (i = 0; i < n_taper; i++) w += val(i);
    for (i = 0; i < n_taper; i++) {
        for (j = 0; j < n_sample; j++)
            taper[i][j] = vec(j, i) * TMath::Sqrt(val(i) / w);
    }
}

void Multitaper::estimate(fftw_complex* seq, double* psd) {
    int i, j;
    for (i = 0; i < n_sample; i++) psd[i] = .0;
    for (i = 0; i < n_taper; i++) {
        for (j = 0; j < n_sample; j++) {
            seqcopy[j][0] = seq[j][0] * taper[i][j];
            seqcopy[j][1] = seq[j][1] * taper[i][j];
        }
        fftw_execute(plan);
        for (j = 0; j < n_sample; j++)
            psd[j] += sq(seqcopy[j][0]) + sq(seqcopy[j][1]);
    }
}
