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


#include <fstream>
#include <iostream>
#include <cmath>
#include <sstream>
#include <gsl/gsl_fft_complex.h>
#include "FritzDPSS.h"
using namespace std;
FritzDPSS::~FritzDPSS() {
    int k;
    for (k=0; k<MaxK+1; k++) delete[] Eigenfunctions[k];
    delete[] Eigenfunctions;
    delete[] Eigenvalues;
    if (!FileP) delete FileP;
}
int FritzDPSS::Calculate(){
//
//  Berechnung der dpss-taper nach Percival/Walden p. 381-384
//
    cout << "Hauptparameter der dpss-Taper" 
	 << "   N = "    << Npoint
	 << "   NW  = "  << NW
	 << "   kmax = " << MaxK << endl;
//
//  Npoint < 4 gibt unsinnige Spektren 
//
    if (Npoint < 4) {
	cout << "Npoint = " << Npoint << " macht keinen Sinn." << endl
	     << " Sollte groesser als 3 sein!" << endl
	     << "Daher werden KEINE Multitaper berechnet." << endl;
	return 0;
    }
//
//  Fall 1: NW ist integer (für multitaper)
//
    if (NWIsInteger == 1) {
	
	//
	//  fuer NW werden nur die Werte 2,3,4 akzeptiert
	//
	if (NW < 2 || NW > 4) {
		cout << "Parameter NW =" << NW 
		<< " ist keiner der Werte 2,3,4" << endl
		<< "Daher werden KEINE Multitaper berechnet." << endl;
		return 0;
	}
    }
    cout << "Berechnung der Taperfunktionen..." << endl;  
    // Deklaration der LINPLUS-Funktion (s.u.)
    int dto_sl(int, double *, double *, double *, int);
    //
    int k, n, niter, converged;
    double sinc, abs_u_orth, gamma;
    Eigenfunctions = new double * [MaxK+1];
    for (k=0; k<MaxK+1; k++) Eigenfunctions[k] = new double[Npoint];
    Eigenvalues = new double [MaxK+1];
    double * b_vec = new double [2*Npoint-1]; // Toeplitz Matrix
    double * v_vec = new double [Npoint];
    double * v_old_vec = new double [Npoint];
    double * u_vec; // wird immer wieder new gemacht
    double anf = pow((double) Npoint,-0.5);
//
//  Toeplitz-Matrix 
//  sin[2 pi N W (i-j)] / pi (i-j) fuer i!=j
//  (Kern der definierenden Integralgleichung der psf's in
//  digitaler Form
//
    for (n=1; n<Npoint; n++) {
	sinc = sin (TwoPiW * n) / (Pi * (double) n);
	b_vec[n] = sinc;
	b_vec[Npoint - 1 + n] = sinc;
    }
//
//  k = 0
//
    b_vec[0] =  2. * Wpar - 1.; // Hauptdiagonalelement
    // cout << "b(0) = " << b_vec[0] << endl;
    for (n=0; n<Npoint; n++) v_vec[n] = anf;
//
//  Iteration fuer k=0
//
    niter = 0;
    int kk;
    int maxit = 3 * (int) sqrt((double) Npoint); 
    do {
	cout << "k = 0 Iteration " << ++niter << endl;
        u_vec = new double [Npoint];
	dto_sl(Npoint, b_vec, v_vec, u_vec, 0); 
	abs_u_orth = 0.;
	for (n=0; n<Npoint; n++) abs_u_orth += u_vec[n] * u_vec[n];
	abs_u_orth = sqrt(abs_u_orth);
	for (n=0; n<Npoint; n++) v_vec[n] = u_vec[n] / abs_u_orth;
	/*
	for (n=0; n<Npoint; n++) cout << "v(" << n << ") = "
				      << v_vec[n] << endl;
	*/
//      Konvergenztest
	if (niter != 1) converged = ConvergenceTest(v_vec, v_old_vec);
	else converged = 0;
	if (converged == 0)
	    for (n=0; n<Npoint; n++) v_old_vec[n] =  v_vec[n];
	delete [] u_vec;
    } while (converged == 0 && niter < maxit);
    if (converged == 0) {
	cout << "Fuer k = 0 keine Konvergenz nach " 
	     << maxit << " Schritten" << endl;
	return 0;
    }
    gamma = -1. / abs_u_orth;
    Eigenvalues[0] = 1 + gamma;
    for (n=0; n<Npoint; n++) Eigenfunctions[0][n] = v_vec[n];
/*
    for (n=0; n<Npoint; n++) cout << "v(" << n << ") = "
				  << v_vec[n] << endl;
*/ 
//
//  Iteration fuer k>0
//
    double * skp_vec = new double[MaxK];
    for (k=1; k<MaxK+1; k++) {
//
//  Neues Diagonalemement der Toeplitz-Matrix
//
	b_vec[0] = 2. * Wpar - Eigenvalues[k-1];	
	niter = 0;
	maxit = (k+3) * (int) sqrt((double) Npoint);
	gamma = 1.;
//      Anfangswert fuer v (Nach Programm von Bell, Percival & Walden)
	int isig = 1;
	int kp1 = k+1;
	int ilow, ihigh;
	for (kk=1; kk<=kp1; kk++) {
	    ilow = ((kk-1) * Npoint) / kp1;
	    ihigh = ((kk * Npoint) / kp1) -1;
	    // cout << "ilow = " << ilow << "ihigh = " << ihigh << endl;
	    for (n = ilow; n <= ihigh; n++) {
		v_vec[n] = isig * anf;
		// cout << "v anf (" << n << ") = " << v_vec[n] << endl;
	    }
	    isig *= -1;
	}
	if(k % 2 == 0 && Npoint % 2 > 0) v_vec[Npoint/2] = 0.;
//
//      Iterationsschleife
//
	do {
	    cout << "k = " << k << "  Iteration " << ++niter << endl;
            u_vec = new double [Npoint]; 
	    dto_sl(Npoint, b_vec, v_vec, u_vec, 0);
	    // cout << "Loesung der Matrixgleichung" << endl;
//          Skalarprodukte
	    for (kk=0; kk<k; kk++) {
		skp_vec[kk] = 0.;
		for (n=0; n<Npoint; n++) 
		    skp_vec[kk] += u_vec[n] * Eigenfunctions[kk][n];
		//cout << kk << ".tes Skalarprodukt = " << skp_vec[kk] << endl;
	    }
//          Orthogonalitaet erzwingen
	    for (n=0; n<Npoint; n++) {
		v_vec[n] = u_vec[n];
		for (kk=0; kk<k; kk++) 
		    v_vec[n] -= skp_vec[kk] * Eigenfunctions[kk][n];
	    }
//          Normierung 
	    abs_u_orth = 0.;
	    for (n=0; n<Npoint; n++) 
		abs_u_orth += v_vec[n] * v_vec[n];
	    abs_u_orth = sqrt(abs_u_orth);
	    for (n=0; n<Npoint; n++) v_vec[n] = v_vec[n] / abs_u_orth;
	    gamma = -1. / abs_u_orth;
//          Konvergenztest
	    if (niter != 1) converged = ConvergenceTest(v_vec, v_old_vec);
	    else converged = 0;
	    delete [] u_vec; // siehe k=0
	    if (converged == 0)
		for (n=0; n<Npoint; n++) v_old_vec[n] =  v_vec[n];
	} while (converged == 0 && niter < maxit);
	if (converged == 0)  {
	    cout << "Fuer k = " << k << " keine Konvergenz nach " 
		 << maxit << " Schritten" << endl;
	    return 0;
	}
	Eigenvalues[k] = Eigenvalues[k-1] + gamma;
	for (n=0; n<Npoint; n++) Eigenfunctions[k][n] = v_vec[n];
/*
	for (n=0; n<Npoint; n++) cout << "v(" << n << ") = "
				      << v_vec[n] << endl;
*/ 
    }
    for (k=0; k<MaxK; k++) Slepian(Eigenfunctions[k], k);
    delete[] skp_vec;
    delete[] v_old_vec;
    delete[] v_vec;
    delete[] b_vec;
    AllReady = 1;
    return 1;
}
double * FritzDPSS::GetTaper(int k) {
    if (AllReady){
	if (k < 0 || k > MaxK) {
	    cout << "Taperindex k = " << k 
		 << " liegt ausserhalb von 0..." << MaxK
		 << endl;
	    return 0;
	}
	return Eigenfunctions[k];
    }
    else {
	cout << "Es stehen keine Taper zur Verfuegung" << endl;
	return 0;
    }
}
double FritzDPSS::GetEigenvalue(int k) {
    if (AllReady){
	if (k < 0 || k > MaxK) {
	    cout << "Taperindex k = " << k 
		 << " liegt ausserhalb von 0..." << MaxK
		 << endl;
	    return 0;
	}
	return Eigenvalues[k];
    }
    else {
	cout << "Es stehen keine Taper zur Verfuegung" << endl;
	return 0;
    }
}
int FritzDPSS::GetSpectrum(complex <double> * timeseries, double * spectrum) {
	//
	// Zeitserie mit Npoint Punkten wird vom aufrufenden Programm zur Verfügung gestellt
	// Spektrum ebenfalls
	//
	gsl_fft_complex_wavetable * wavetable = gsl_fft_complex_wavetable_alloc(Npoint);
	gsl_fft_complex_workspace * workspace= gsl_fft_complex_workspace_alloc(Npoint);
	int ik;
	int ip;
	complex <double> * kth_spectrum = new complex <double> [Npoint];
	double * kth_spectrum_gsl = new double[2*Npoint];
	double * kth_spectrum_gsl_iter;
	double fft_real, fft_imag;
	for (ip=0; ip<Npoint; ip++) spectrum[ip] = 0.;
	for(ik=0; ik<MaxK; ik++){
		kth_spectrum_gsl_iter = kth_spectrum_gsl;
		for (ip=0; ip<Npoint; ip++) {
			// multiply timeseries with taper
			kth_spectrum[ip] = timeseries[ip] * Eigenfunctions[ik][ip];
			*kth_spectrum_gsl_iter++ = kth_spectrum[ip].real();
			*kth_spectrum_gsl_iter++ = kth_spectrum[ip].imag();
		}
		// get tapered spectrum
		gsl_fft_complex_forward(kth_spectrum_gsl, 1, Npoint, wavetable, workspace);
		kth_spectrum_gsl_iter = kth_spectrum_gsl;
		// add up absolute values
		for (ip=0; ip<Npoint; ip++) {
			fft_real = *kth_spectrum_gsl_iter++;
			fft_imag = *kth_spectrum_gsl_iter++;
			spectrum[ip] += fft_real * fft_real + fft_imag * fft_imag;
		}
	}
	// normalize and copy back for reordering
	double dk = (double) MaxK;
	for (ip=0; ip<Npoint; ip++) kth_spectrum_gsl[ip] = spectrum[ip] / dk;
	// reorder by frequency (see /usr/share/doc/gsl-ref-psdoc/gsl-ref.ps.gz pp.141f)
	// this is valid only if Npoint is even!
	kth_spectrum_gsl_iter = kth_spectrum_gsl + Npoint/2; 
	for (ip=0; ip<Npoint/2; ip++) spectrum[ip] = *kth_spectrum_gsl_iter++;
	kth_spectrum_gsl_iter = kth_spectrum_gsl; 
	for (ip=Npoint/2; ip<Npoint; ip++) spectrum[ip] = *kth_spectrum_gsl_iter++;
	delete [] kth_spectrum_gsl;
	delete [] kth_spectrum;
	gsl_fft_complex_wavetable_free(wavetable);
	gsl_fft_complex_workspace_free(workspace);
	return 1;
    }
int FritzDPSS::ConvergenceTest(double * v_vec, double * v_old_vec){         
    int n;    
    double sumtest = 0.;
    double difftest = 0.;
    for (n=0; n<Npoint; n++) {
	// cout << v_vec[n] << "  " << v_old_vec[n] << endl;
	sumtest  += pow(v_vec[n] + v_old_vec[n], 2);
	difftest += pow(v_vec[n] - v_old_vec[n], 2);
    }
    if (sumtest > difftest) sumtest = difftest;
    // cout << "Test = " << sumtest << endl;
    return (sumtest < Tol);    
}
void FritzDPSS::Slepian(double * v, int k){
//
//  Vorzeichen der dpss nach Slepian-Konvention, siehe
//  Bell, B., Percival, D.B. and Walden, A.T.
//  "Calculating Thomson's Spectral Multitapers by Inverse Iteration",
//  J. Comput. and Graph. Stat., 1993. (Section 1.2.)
//
    int n;
    double dsum=0., dwsum=0.;
    for (n=0; n<Npoint; n++) {
	dsum += v[n];
	dwsum += v[n] * (double) (Npoint - 1 - 2 * n);
    }
    if(((k % 2 == 0) && (dsum < 0.)) || ((k % 2 == 1) && (dwsum < 0.)))
	for (n=0; n<Npoint; n++) v[n] *= -1.;
}
//
//
//
int FritzDPSS::StoreOnDisk(){
    int k;
    if (!AllReady){
	cout << "Taper sind noch nicht fertig. Es wird nichts geschrieben."
	     << endl;
	return 0;
    }
//
//      File oeffnen
//
    MakeFileName();
    FileP = new fstream(FileName.c_str(), ios::out | ios::binary); 
    if ((FileP->rdstate() & ifstream::failbit ) != 0 ) {
	cout << "Datei " << FileName.c_str() << " laesst sich nicht oeffnen. "
	     << endl;
	return 0;
    }
    cout << "tapers werden auf Datei " << FileName << " abgelegt" << endl;
    FileP->write((char*) &Npoint, sizeof(int));
    FileP->write((char*) &NW, sizeof(int));
    FileP->write((char*) &MaxK, sizeof(int));
    for (k=0; k<MaxK+1; k++) {
	FileP->write((char*) Eigenfunctions[k], sizeof(double)*Npoint);
    }
    FileP->close();
    cout << "Schreibvorgang beendet" << endl;
    return 1;
}
//
//
//
int FritzDPSS::GetFromDisk(){
    int k;
//
//      File oeffnen
//
    MakeFileName();
    FileP = new fstream(FileName.c_str(), ios::in | ios::binary); 
    if ((FileP->rdstate() & ifstream::failbit ) != 0 ) {
	cout << "Datei " << FileName << " laesst sich nicht oeffnen. "
	     << endl;
	return 0;
    }
    cout << "tapers werden von Datei " << FileName << " gelesen" << endl;
    int npoint=0, nw=0, maxk=0;
    FileP->read((char*) &npoint, sizeof(int));
    if (npoint != Npoint) {
	cout << "Inkonsistenz beim Lesen von Npoint: "
	     << " Eingelesen: " << npoint
	     << " Sollwert: " << Npoint << endl;
	return 0;
    }
    FileP->read((char*) &nw, sizeof(int));
    if (nw != NW) {
	cout << "Inkonsistenz beim Lesen von NW: "
	     << " Eingelesen: " << nw
	     << " Sollwert: " << NW << endl;
	return 0;
    }
    FileP->read((char*) &maxk, sizeof(int));
    if (maxk != MaxK) {
	cout << "Inkonsistenz beim Lesen von NW: "
	     << " Eingelesen: " << nw
	     << " Sollwert: " << NW << endl;
	return 0;
    }
    Eigenvalues = new double [MaxK+1];
    Eigenfunctions = new double * [MaxK+1];
    for (k=0; k<MaxK+1; k++) {
	Eigenfunctions[k] = new double[Npoint];
	FileP->read((char*) Eigenfunctions[k], sizeof(double)*Npoint);
    }
    FileP->close();
    cout << "Lesevorgang beendet" << endl;
    AllReady = 1;
    return 1;
}
//
//
//
void FritzDPSS::MakeFileName(){
    string s = "TAP";
    stringstream ss;
    ss << s << Npoint << "NW" << NW << "MaxK" << MaxK << ".tap";
    FileName = ss.str();
}
//
// Routine der LINPLUS-Bibliothek fuer Toeplitz-Matrizen
// Da die Matrix negative Eigenwerte hat, kann nicht die symmetrische
// Version benutzt werden, wo A positiv definit sein muesste.
//
//**********************************************************************
int dto_sl ( int n, double a[], double b[], double x[], int job )

//**********************************************************************
//
//  Purpose:
//
//    DTO_SL solves a DTO system.
//
//  Discussion:
//
//    The DTO storage format is used for a Toeplitz matrix, which is constant
//    along diagonals.  Thus, in an N by N Toeplitz matrix, there are at most 
//    2*N-1 distinct entries.  The format stores the N elements of the first
//    row, followed by the N-1 elements of the first column (skipping the
//    entry in the first row).
//
//  Modified:
//
//    23 September 2003
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[2*N-1], the DTO matrix.
//
//    Input, double B[N] the right hand side vector.
//
//    Input, int JOB,
//    0 to solve A*X=B,
//    nonzero to solve A'*X=B.
//
//    Output, double DTO_SL[N], the solution vector.  X and B may share the
//    same storage.
//
{
  double *c1;
  double *c2;
  int i;
  int nsub;
  double r1;
  double r2;
  double r3;
  double r5;
  double r6;

  if ( n < 1 )
  {
    return 0;
  }

//
//  Solve the system with the principal minor of order 1.
//
  r1 = a[0];
  x[0] = b[0] / r1;

  if ( n == 1 )
  {
    return 1;
  }

  c1 = new double[n-1];
  c2 = new double[n-1];
//
//  Recurrent process for solving the system with the Toeplitz matrix.
//
  for ( nsub = 2; nsub <= n; nsub++ )
  {
//
//  Compute multiples of the first and last columns of the inverse of
//  the principal minor of order NSUB.
//
    if ( job == 0 )
    {
      r5 = a[n+nsub-2];
      r6 = a[nsub-1];
    }
    else
    {
      r5 = a[nsub-1];
      r6 = a[n+nsub-2];
    }

    if ( 2 < nsub )
    {
      c1[nsub-2] = r2;

      for ( i = 1; i <= nsub-2; i++ )
      {
        if ( job == 0 )
        {
          r5 = r5 + a[n+i-1] * c1[nsub-i-1];
          r6 = r6 + a[i] * c2[i-1];
        }
        else
        {
          r5 = r5 + a[i] * c1[nsub-i-1];
          r6 = r6 + a[n+i-1] * c2[i-1];
        }
      }
    }

    r2 = - r5 / r1;
    r3 = - r6 / r1;
    r1 = r1 + r5 * r3;

    if ( 2 < nsub )
    {
      r6 = c2[0];
      c2[nsub-2] = 0.0;

      for ( i = 2; i <= nsub-1; i++ )
      {
        r5 = c2[i-1];
        c2[i-1] = c1[i-1] * r3 + r6;
        c1[i-1] = c1[i-1] + r6 * r2;
        r6 = r5;
      }
    }

    c2[0] = r3;
//
//  Compute the solution of the system with the principal minor of order NSUB.
//
    if ( job == 0 )
    {
      r5 = 0.0;
      for ( i = nsub-1; 1 <= i; i-- )
      {
        r5 = r5 + a[n+nsub-i-1] * x[i-1];
      }
    }
    else
    {
      r5 = 0.0;
      for ( i = nsub-1; 1 <= i; i-- )
      {
        r5 = r5 + a[nsub-i] * x[i-1];
      }
    }

    r6 = ( b[nsub-1] - r5 ) / r1;

    for ( i = 0; i < nsub-1; i++ )
    {
      x[i] = x[i] + c2[i] * r6;
    }
    x[nsub-1] = r6;
  }

  delete [] c1;
  delete [] c2;

  return 1;
}
