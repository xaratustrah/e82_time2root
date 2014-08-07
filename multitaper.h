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



#ifndef _MULTITAPER_H_
#define _MULTITAPER_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <TMath.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <fftw3.h>

using namespace std;

class Multitaper {
    public:
        Multitaper(int n) :
            n_sample(n), n_taper(6), nwidth(4) {build(); load();}
        Multitaper(int n, int k) :
            n_sample(n), n_taper(k), nwidth(4) {build(); load();}
        Multitaper(int n, int k, int nw) :
            n_sample(n), n_taper(k), nwidth(nw) {build(); load();}
        ~Multitaper() {
            for (int i = 0; i < n_taper; i++) delete[] taper[i];
            delete[] taper;
            fftw_destroy_plan(plan);
            fftw_free(seqcopy);
        }
        void estimate(fftw_complex*, double*);

    private:
        const int n_sample;
        const int n_taper; // less than 2*nwidth
        const int nwidth;
        double** taper;
        fftw_complex* seqcopy;
        fftw_plan plan;
        double sq(double x) {return x*x;}
        void build();
        void load();
        void generate();
};

#endif // _MULTITAPER_H_
