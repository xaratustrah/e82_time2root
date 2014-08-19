//
// (c) Copyright:
// F. Nolden 2008 - 2009
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


#ifndef _IQTDATA_H
#define _IQTDATA_H
#include<string>
#include<fstream>
#include<complex>
#include <gsl/gsl_fft_complex.h>
/*
C++ - Paket zur externen Auswertung von Messdaten, die mit dem
Tektronix RSA3303B Spektrum - Analysator aufgenommen wurden.
*/
using namespace std;
enum DimensionT {Hertz, Sekunden, dBm, mWatt, Dimensionslos};
// notwendige Deklarationen, Definitionen weiter unten
class FrameHeader__TekRSA3303B;
class FrameData__TekRSA3303B;
class ExtCorrData__TekRSA3303B;
//
// Diese Klasse verwaltet die kompletten Daten aus dem *.IQ-File 
//
class  IQData__TekRSA3303B {
	public:
	IQData__TekRSA3303B();
	~IQData__TekRSA3303B();
	int ReadFile(string s);
        int GetIQ(complex <double>&, int);
        int GetTimeStamp(double&, int);
	int GetDeltaT(double&);
        int Bins() const {return mBins;};
        int BlockSize() const {return mBlockSize;};
        int ValidFrames() const {return mValidFrames;};
        double FrameLength() const {return mFrameLength;};
        double GainOffset() const {return mGainOffset;};
        double CenterFrequency() const {return mCenterFrequency;};
        double Span() const {return mSpan;};
        const char* DateTime() const {return mDateTime;};
        ifstream * IQFileP() const {return mIQFileP;};
	private:
        IQData__TekRSA3303B (const IQData__TekRSA3303B&); /* never defined */
        IQData__TekRSA3303B& operator= (const IQData__TekRSA3303B&); /* never defined */
        // private methods
	int Aufdroeseln(string);
        int DecodeInt(int&, string);
	int DecodeDouble(double&, string);
        void DecodeDateTime(char*, string);
        int GetFrameNumber(int&, int);
        // private members
        // aus IQ-File direkt eingelesen
        char mDateTime[20];
        int mBins;
        double mMaxInputLevel;
        double mLevelOffset;
        double mCenterFrequency; // Hz
        double mFrequencyOffset;
        double mSpan; // Hz
        int mBlockSize;
        int mValidFrames;
        double mFramePeriod;
        double mUnitPeriod;
        double mFrameLength; // s
        double mGainOffset; // dB
        int mMultiFrames;
        int mMultiAddr;
	double mIOffset;
	double mQOffset;
        // sonstige
        int mReadCalled;
	ifstream * mIQFileP;
	double mIQScale;
	double mFirstFrequency;
	double mDeltaFrequency;
	//
        int mValidScans;
        int mValidFrequencies;
	// fuer die Korrektur im Frequenzbereich
	int mActualCorrectedScanNumber;
	complex <double> * ActualCorrectedScanP;
	double * mDataFuerFFT;
	gsl_fft_complex_wavetable * mWavetable;
	gsl_fft_complex_workspace * mWorkspace;
        // Frames: Header und Daten
        FrameHeader__TekRSA3303B ** mFrameHeaderPP;
        FrameData__TekRSA3303B ** mFrameDataPP;
	// Korrektur
	complex <double> * mCorrectionP;
};
//
//      abstrakte Basisklasse fuer frame-header und -daten
//
/*abstract*/ class FramePart__TekRSA3303B{
	public:
	FramePart__TekRSA3303B(IQData__TekRSA3303B * iq_data_p) :
	  mIQDataP(iq_data_p) {mIQDataFile = mIQDataP->IQFileP();};
	virtual ~FramePart__TekRSA3303B() {};
	virtual void Read() = 0;
	protected:
	short GetNextShort();
        IQData__TekRSA3303B * mIQDataP;
        ifstream * mIQDataFile;
};
//
//      Frame-Header-Klasse
//
class FrameHeader__TekRSA3303B : public FramePart__TekRSA3303B {
	public:
	FrameHeader__TekRSA3303B(IQData__TekRSA3303B * iq_data_p) :
	  FramePart__TekRSA3303B(iq_data_p), mLastFrame(0) {};
	virtual void Read();
	void PrintAll();
	short ValidA() const  {return mValidA;}; 
	short ValidP() const  {return mValidP;}; 
	short ValidI() const  {return mValidI;}; 
	short ValidQ() const  {return mValidQ;}; 
	short Bins() const  {return mBins;}; 
	short Triggered() const  {return mTriggered;}; 
        short Overload() const {return mOverLoad;};
	short LastFrame() const  {return mLastFrame;};
	unsigned long LongTicks() const  {return mLongTicks;};
	double Skalenfaktor() const  {return mSkalenfaktor;};
	private:
	// virtual ~FrameHeader__RSA3303B(); default
        // never defined
	FrameHeader__TekRSA3303B (const FrameHeader__TekRSA3303B&);
	FrameHeader__TekRSA3303B& operator=(const FrameHeader__TekRSA3303B&);
        //
	short mValidA;
	short mValidP; 
	short mValidI; 
	short mValidQ; 
	short mBins;
	short mTriggered; 
	short mOverLoad; 
	short mLastFrame;
	unsigned long mLongTicks;
        double mSkalenfaktor;
};
//
//      Frame-Daten-Klasse
//
class FrameData__TekRSA3303B : public FramePart__TekRSA3303B {
	public:
	FrameData__TekRSA3303B(IQData__TekRSA3303B *);
	virtual ~FrameData__TekRSA3303B();
	virtual void Read();
        short GetI(int index) {return IVectorP[index];};
        short GetQ(int index) {return QVectorP[index];};
	private:
        // never defined
	FrameData__TekRSA3303B (const FrameData__TekRSA3303B&);
	FrameData__TekRSA3303B& operator=(const FrameData__TekRSA3303B&);
        //
	short * IVectorP;
	short * QVectorP;
};
//
//      Extended Correction Daten-Klasse
//
class ExtCorrData__TekRSA3303B : public FramePart__TekRSA3303B {
	public:
	ExtCorrData__TekRSA3303B(IQData__TekRSA3303B *);
	virtual ~ExtCorrData__TekRSA3303B();
	virtual void Read();
        char GetA(int index) {return AVectorP[index];};
        char GetP(int index) {return PVectorP[index];};
	private:
        // never defined
	ExtCorrData__TekRSA3303B (const ExtCorrData__TekRSA3303B&);
	ExtCorrData__TekRSA3303B& operator=(const ExtCorrData__TekRSA3303B&);
        //
	char * AVectorP;
	char * PVectorP;
};
#endif //  _IQTDATA_H
