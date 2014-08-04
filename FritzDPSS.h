//
// (c) Copyright:
// F. Nolden 2008 - 2009
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


#ifndef _FritzDPSS_H
#define _FritzDPSS_H
#include <fstream>
#include <string>
#include <complex>
using namespace std;
class FritzDPSS{
//
//  oeffentliche Methoden
//
 public:
    FritzDPSS(int npoint, int nw) 
	: Npoint(npoint), NWIsInteger(1), NW(nw), MaxK(2 * NW - 2), AllReady(0),
	Pi(3.141592653), TwoPi(6.283185306), 
	TwoPiW(TwoPi * (double) nw / (double) npoint), 
	Wpar((double) nw / (double) npoint), Tol(5e-7) {};
    FritzDPSS(int npoint, double nw, int maxk) 
	: Npoint(npoint), NWIsInteger(0), NW(0), MaxK(maxk), AllReady(0),
	Pi(3.141592653), TwoPi(6.283185306), 
	TwoPiW(TwoPi * nw / (double) npoint),
	Wpar((double) nw / (double) npoint), Tol(5e-7) {};
    ~FritzDPSS();
    int Calculate();
    const int GetNpoint() const {return Npoint;};
    const int GetNW() const {return NW;};
    const int GetMaxK() const {return MaxK;};
    int IsReady() const {return AllReady;};
    double * GetTaper(int);
    double GetEigenvalue(int);
    int GetSpectrum(complex <double> *, double *);
    void MakeFileName();
    int StoreOnDisk();
    int GetFromDisk();
//
//  privater Teil
//
 private:
    FritzDPSS (const FritzDPSS&); /* never defined */
    FritzDPSS& operator= (const FritzDPSS&);/* never defined */
    int ConvergenceTest(double *, double *);
    void Slepian(double *, int);
    const int Npoint;
    const int NWIsInteger;
    const int NW;
    const int MaxK;
    int AllReady;
    const double Pi;
    const double TwoPi;
    const double TwoPiW;
    const double Wpar;
    const double Tol;
    double ** Eigenfunctions;
    double * Eigenvalues;
    string FileName;
    fstream * FileP;
};
#endif // _FritzDPSS_H
