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


/*-----------------------------------------*
 | TIQ Reader for Tektronix RSA5000 Series |
 | Reader for Header   		           |
 *-----------------------------------------*/


#ifndef _TIQ_H_
#define _TIQ_H_

#define cFrmPt  1024

#include <complex>
#include <fftw3.h>

class TDatime;
class TComplex;

typedef struct
{
    const char*   File;
    int     Offset;
    int     NumPt;
    double  Scaling;
    double  Intvl; /* second */
    double  CenFreq; /* Hz */
    double  Span; /* Hz */
    TDatime DaTm;
} Info_t;

using namespace std;

/* read general information from a .tiq file header, and store in Info_t */
bool SetInfo(FILE*, Info_t*);
/* apply FFT/IFFT to a TComplex array */
bool FFT(TComplex*, short);
/* apply multitaper method for spectral density estimation */
bool Multitaper(complex<double>*, double*);

#endif
