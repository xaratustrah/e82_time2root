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
 | This code provides two different implementations of multitaper   |
 | estimation, which is widely used in digital signal analysis, and |
 | benchmarks their calculation speeds.                             |
 *------------------------------------------------------------------*/

#include <TRandom.h>
#include "FritzDPSS.h"
#include "multitaper.h"

int main(void) {
    int i, j;
    clock_t t;
    const int n_frame = 10000;
    const int n_point = 1024;
    double* spectrum = new double [n_point];
    complex<double>* fn_signal = new complex<double> [n_point];
    fftw_complex* xc_signal = new fftw_complex [n_point];
    TRandom event;

    cout << "FN's multitaper..." << endl;
    t=clock();
    FritzDPSS fn_dpss(n_point, 4);
    fn_dpss.GetFromDisk();
    for (i = 0; i < n_frame; i++) {
        for (j = 0; j < n_point; j++)
            fn_signal[j] = (event.Gaus(), event.Gaus());
        fn_dpss.GetSpectrum(fn_signal, spectrum);
    }
    t = clock() - t;
    cout << "it took " << (double)t/CLOCKS_PER_SEC << " seconds." << endl;

    cout << "XC's multitaper..." << endl;
    t=clock();
    Multitaper xc_multitaper(n_point);
    for (i = 0; i < n_frame; i++) {
        for (j = 0; j < n_point; j++) {
            xc_signal[j][0] = event.Gaus();
            xc_signal[j][1] = event.Gaus();
        }
        xc_multitaper.estimate(xc_signal, spectrum);
    }
    t = clock() - t;
    cout << "it took " << (double)t/CLOCKS_PER_SEC << " seconds." << endl;

    delete[] spectrum;
    delete[] fn_signal;
    delete[] xc_signal;
    return 0;
}
