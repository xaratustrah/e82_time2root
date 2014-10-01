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


#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include "iqtdata.h"
using namespace std;
//
//               IQData__TekRSA3303B::IQData__TekRSA3303B()
//
IQData__TekRSA3303B::IQData__TekRSA3303B() : 
        mBins(0),
        mMaxInputLevel(100.),
        mLevelOffset(100.),
        mCenterFrequency(0.),
        mFrequencyOffset(1.0e10),
        mSpan(1.0e10),
        mBlockSize(0),
        mValidFrames(0),
        mFramePeriod(0.),
        mUnitPeriod(0.),
        mFrameLength(0.),
        mGainOffset(0.),
        mMultiFrames(0),
        mMultiAddr(0),
	mIOffset(0.),
	mQOffset(0.),
        mReadCalled(0),
	mActualCorrectedScanNumber(-1) {
	mIQFileP = new ifstream();
	};
//
//               IQData__TekRSA3303B::~IQData__TekRSA3303B()
//
IQData__TekRSA3303B::~IQData__TekRSA3303B() {
	int i;
	delete []  ActualCorrectedScanP;
	delete [] mCorrectionP;
	gsl_fft_complex_wavetable_free(mWavetable);
	gsl_fft_complex_workspace_free(mWorkspace);
	delete [] mDataFuerFFT;
	delete mIQFileP;
	if (mFrameHeaderPP) {
		for (i=0; i < mValidFrames; i++) {
			delete mFrameHeaderPP[i];
			delete mFrameDataPP[i];
		} 
	}
        delete []  mFrameHeaderPP;
	delete []  mFrameDataPP;
}
//
//               IQData__TekRSA3303B::ReadFile()
//
int IQData__TekRSA3303B::ReadFile(string filename) {
	if (mReadCalled) {
	  cout << "Zweimal lesen geht nicht." << endl;
	  return 0;
	}
	mReadCalled = 1;
//
//      File oeffnen
//
	mIQFileP->open(filename.c_str(), ios::in | ios::binary);
	if (mIQFileP->is_open())
		cout << "File ist offen" << endl;
	else {
		cout << "File ist nicht offen!" << endl;
		return 0;
	}
	
	mIQFileP->seekg (0, ios::beg);
//
//      Erste 5 Bytes
//
	char some[3] = {};
	mIQFileP->read(some, 2);
	first = (some[0] - '0')* 10 + (some[1] - '0');
	string what(some, 2);
	stringstream wandler(some);
	int first=0;
	wandler >> first;
	if (first != 40) {
		cout << "keine 40 am Anfang" << endl;
		return 0;
	}
	// cout << "40 am Anfang: o.k." << endl;
	mIQFileP->read(some, 3);
	what = some;
	wandler.clear();
	wandler.str(what);
	wandler >> first;
	// cout << "Zahl der Bytes" << first << endl;
	//
	// Rest des Headers
	//
	int ibyte;
	what.clear();
	for (ibyte = 0; ibyte < first; ibyte++) {
		mIQFileP->read(some, 1);
		what += some[0];
	}
	wandler.clear();
	wandler.str(what);
	string zeile;
		wandler >> zeile;
	while (!wandler.eof()) {
		wandler >> zeile;
		if (!wandler.fail()) {
			Aufdroeseln(zeile);
		}
	}
	cout << "header gelesen" << endl;
	//
	// Vektoren fuer Frame-Header und -Daten erzeugen
	//
	int i;
	mFrameHeaderPP = new FrameHeader__TekRSA3303B * [mValidFrames];
	mFrameDataPP = new FrameData__TekRSA3303B * [mValidFrames];
	for (i=0; i < mValidFrames; i++) {
	  mFrameHeaderPP[i] = new FrameHeader__TekRSA3303B(this);
	  mFrameDataPP[i] = new FrameData__TekRSA3303B(this);
	}
	//
	// alle Blocks einlesen
	// Header auf Besonderheiten checken und Trigger identifizieren
	//
	int frames_read;
	FrameHeader__TekRSA3303B * actual_header_p;
	unsigned long ticks_at_trigger = 0;
	cout << "Bitte warten..." << endl;
	for (frames_read=0; frames_read < mValidFrames; frames_read++) {
	// Header auf Besonderheiten checken und Trigger identifizieren
	  actual_header_p = mFrameHeaderPP[frames_read];
	  actual_header_p->Read();
	  //actual_header_p->PrintAll();
	  if (actual_header_p->Overload())
	    cout << "In  Frame " << frames_read 
		 << " wurde der maximale Eingangspegel von "
		 << mMaxInputLevel << " dBm ueberschritten" << endl;
	  if (actual_header_p->LastFrame())
	    cout << "Blockende in Frame " << frames_read << endl;
	  if (!ticks_at_trigger) {
	    if (actual_header_p->Triggered())
	      ticks_at_trigger = actual_header_p->LongTicks();
	      }
        // Daten lesen
	  mFrameDataPP[frames_read]->Read();
	}
	cout << frames_read << " Frames gelesen" << endl;
	//
	// jetzt kommt der dümmliche correction data block
	//
	cout << "correction data" << endl;
	FrameHeader__TekRSA3303B * correction_header_p = new FrameHeader__TekRSA3303B(this);
	correction_header_p->Read();
	// correction_header_p->PrintAll();
	delete correction_header_p; // der Header ist fuer nichts gut
	FrameData__TekRSA3303B * correction_data_p = new FrameData__TekRSA3303B(this);
	correction_data_p->Read();
	cout << "Correction Data Block gelesen" << endl;
	//
	// jetzt kommt der noch dümmlichere "dummy header" 40000
	//
	char cd;
	what.clear();
	for (i=0; i<5; i++) {
	cd = mIQFileP->get();
	what += cd;
	}
	cout << "Dummy Header " << what << endl;
	// cout << "get pointer position: " << mIQFileP->tellg() << endl;
	//
	// und jetzt der noch dümmlichere "extended correction block"
	//
	ExtCorrData__TekRSA3303B * ext_correction_data_p = new ExtCorrData__TekRSA3303B(this);
	ext_correction_data_p->Read();
	cout << "Extended Correction Data Block gelesen" << endl;
	//
	// und was nicht kommt, ist Date and Time!
	//
	mIQFileP->close();
	//
	// Berechnung der linearen Korrekturfaktoren
	//
	mCorrectionP = new complex <double> [1024];
	double ac, pc;
	double pid180 = 0.01745329252;
	int gcd, egcd;
	for (i=0; i<1024; i++) {
		gcd = correction_data_p->GetI(i);
		gcd = gcd << 8;
		egcd = ext_correction_data_p->GetA(i);
		egcd = egcd & 0x000000ff;
		gcd = gcd | egcd;
		ac = -((double) gcd) / (128*256); // - Vorzeichen
		ac = pow(10., ac/20.); // Umrechnung auf linear
		gcd = correction_data_p->GetQ(i);
		gcd = gcd << 8;
		egcd = ext_correction_data_p->GetP(i);
		egcd = egcd & 0x000000ff;
		gcd = gcd | egcd;
		pc = ((double) gcd ) / (128*256); // + Vorzeichen
		/*
		cout << i << " "
		<< mAmplitudeCorrectionP[i]  << " dB"
		<< mPhaseCorrectionP[i]  << " Grad" << endl;
		*/
		pc *= pid180; // Korrektur im Bogenmass
		mCorrectionP[i] = complex <double> (ac * cos(pc), ac * sin(pc));
	}
	delete ext_correction_data_p;
	delete correction_data_p;
	mReadCalled = 1;
	// es gibt im RSA3303B-Handbuch (Seite 3-272) zwei Skalierungsvorschriften,
	// eine für die Amplitude in dBm und eine für den Spannungwert in V. Letztere wird in
	// mIQScale benutzt. Sie ist nur in einfacherer Form programmiert und laesst sich aus
	// P = U^2 / (2 ZL) herleiten. Der Term -1. zum Schluss ergibt sich mit -1=-3+2 einmal
	// aus der Umrechnung von mW in W und weiter aus der Multiplikation mit 2*50 Ohm.
		mIQScale = pow(10., ((mGainOffset+mMaxInputLevel+mLevelOffset)/10-1.)/2.);
	//	mIQScale = sqrt(pow(10., (mGainOffset + mMaxInputLevel+ mLevelOffset) / 10) / 20 * 2);
	//
	// Vektoren fuer FFT und Korrektur
	//
	mDataFuerFFT = new double [2048]; // doppelter Platz fuer Re und Im
	mWavetable = gsl_fft_complex_wavetable_alloc(1024);
	mWorkspace= gsl_fft_complex_workspace_alloc(1024);
	ActualCorrectedScanP = new complex <double> [1024];
	return 1;
}
//
//               IQData__TekRSA3303B::GetIQ()
//
int IQData__TekRSA3303B::GetIQ(complex <double>& amplitude, int npt) {
	if (!mReadCalled) {
		cout << "File ist noch nicht eingelesen" << endl;
		return 0;
	}
	int max_npt = mValidFrames*1024;
	if (npt < 0 || npt >= max_npt) {
		cout << "in GetIQ: " << npt << " ist ausserhalb von (0, "
			<< max_npt << ")" << endl;
		return 0;
	}
	//
	// jetzt muss die Korrektur auf den Frequenzgang beruecksichtigt werden.
	// Das geschieht Frame fuer Frame folgendermassen:
	// 1. FFT des gesamten Frames
	// 2. Anbringen der Korrektur im Frequenzbereich
	// 3. Ruecktransformation
	//
	int frame_nummer = npt / 1024;
	int frame_index = npt % 1024;
	int ip;
	double i_v, q_v;
	complex <double> amp;
	double * data_fuer_fft_iter;
	if (frame_nummer != mActualCorrectedScanNumber) {
		// cout << "frame-Nummer = " << frame_nummer << endl;
		data_fuer_fft_iter = mDataFuerFFT;
		for (ip=0; ip<1024; ip++) {
			*data_fuer_fft_iter++ = (double) mFrameDataPP[frame_nummer]->GetI(ip);
			*data_fuer_fft_iter++ = (double) mFrameDataPP[frame_nummer]->GetQ(ip);
		}
		gsl_fft_complex_forward(mDataFuerFFT, 1, 1024, mWavetable, mWorkspace);
		data_fuer_fft_iter = mDataFuerFFT;
		for (ip=0; ip<1024; ip++) {
			i_v = *data_fuer_fft_iter++;
			q_v = *data_fuer_fft_iter;
			amp = complex <double> (i_v,q_v);
			// Korrektur
			amp *= mCorrectionP[ip];
			data_fuer_fft_iter--;
			*data_fuer_fft_iter++ = amp.real();
			*data_fuer_fft_iter++ = amp.imag();
		}
		gsl_fft_complex_inverse(mDataFuerFFT, 1, 1024, mWavetable, mWorkspace);
		data_fuer_fft_iter = mDataFuerFFT;
		for (ip=0; ip<1024; ip++) {
			i_v = *data_fuer_fft_iter++;
			q_v = *data_fuer_fft_iter++;
			ActualCorrectedScanP[ip] = complex <double> (i_v,q_v);
		}
		mActualCorrectedScanNumber = frame_nummer;
	}
	amplitude =  mIQScale * ActualCorrectedScanP[frame_index];
	return 1;	
}
/*
//
//               IQData__TekRSA3303B::GetTimeStamp()
//
int IQData__TekRSA3303B::GetTimeStamp(double& time_stamp, int scan_nummer) {
        int frame_nummer;
	if (!GetFrameNumber(frame_nummer, scan_nummer)) return 0;
	if (mMultiFrames > 1) frame_nummer = scan_nummer * mMultiFrames;
	time_stamp = mFrameHeaderPP[frame_nummer]->LongTicks() * mUnitPeriod
	           - mTriggerTime;
	return 1;
}
*/
//
//               IQData__TekRSA3303B::Aufdroeseln()
//
int IQData__TekRSA3303B::Aufdroeseln(string zeile) {
	//
	// Vor dem "=" steht der Name,
	// dann kommt der Inhalt
	//
	size_t  equal_position = zeile.find("=");
	if (equal_position == string::npos) {
		cout << "kein Gleichheitszeichen in Header-Zeile:" << zeile << endl;
		return 0;
	}
	size_t length = zeile.length();
	string name = zeile.substr(0,equal_position);
	string inhalt = zeile.substr(equal_position+1, length-equal_position-1);
	// cout << "zeile: " << zeile << ", Name :" << name << ", Inhalt : " << inhalt << endl;
	//
	// Jetzt werden Namen und Inhalte interpretiert
	//
	if (name == "Type") {
	  if      (inhalt == "RSA3408AIQT")  {
		cout << "iqt-File von Analyzer RSA3408A" << endl;
	  }
	  else {
		cout << "Ungueltiger Datentyp" << endl;
		return 0;
	  }
	}
	else if (name == "FrameReverse") {
	  if  (inhalt != "Off") {
		cout << "FrameReverse steht NICHT auf Off" << endl;
		return 0;
	  }
	}
	else if (name == "FramePadding") {
	  if      (inhalt !="Before") {
		cout << "FramePadding steht NICHT auf Before" << endl;
		return 0;
	  }
	}
	else if (name == "Band") {} // Uninteressant
	else if (name == "MemoryMode") {} // Uninteressant
	else if (name == "FFTPoints") {
	  if  (inhalt != "1024") {
		cout << "FFTPoints steht NICHT auf 1024" << endl;
		return 0;
	  }
	}
	else if (name == "Bins") {
	  if (DecodeInt(mBins, inhalt)) {
	    cout << mBins << " gueltige Bin-Daten" << endl;
	  }
	  else return 0;
	}
	
	else if (name == "MaxInputLevel") {
	  if (DecodeDouble(mMaxInputLevel, inhalt)) {
	    cout  << "Maximaler Eingangspegel " << mMaxInputLevel 
	          <<  " dBm" << endl;
	  }
	  else return 0;
	}
	else if (name == "LevelOffset") {
	  if (DecodeDouble(mLevelOffset, inhalt)) {
	    cout  << "Pegel-Offset " << mLevelOffset
	          <<  " dBm" << endl;
	  }
	  else return 0;
	}
	else if (name == "CenterFrequency") {
	  if (DecodeDouble(mCenterFrequency, inhalt)) {
	    cout  << "Mittenfrequenz " << mCenterFrequency
	          <<  " Hz" << endl;
	  }
	  else return 0;
	}
	else if (name == "FrequencyOffset") {
	  if (DecodeDouble(mFrequencyOffset, inhalt)) {
	    cout  << "Frequenz-Offset " << mFrequencyOffset    
	          <<  " Hz" << endl;
	  }
	  else return 0;
	}
	else if (name == "Span") {
	  if (DecodeDouble(mSpan, inhalt)) {
	    cout  << "Frequenzbandbreite " << mSpan    
	          <<  " Hz" << endl;
	  }
	  else return 0;
	}
	else if (name == "BlockSize") {
	  if (DecodeInt(mBlockSize, inhalt)) {
	    cout << "Blocklaenge " << mBlockSize << endl;
	  }
	  else return 0;
	}
	else if (name == "ValidFrames") {
	  if (DecodeInt(mValidFrames, inhalt)) {
	    cout << "gueltige Frames " << mValidFrames << endl;
	  }
	  else return 0;
	}
	else if (name == "FramePeriod") {
	  if (DecodeDouble(mFramePeriod, inhalt)) {
	    cout  << "Frameabstand " << mFramePeriod 
		  << " s" << endl;
	  }
	  else return 0;
	}
	else if (name =="UnitPeriod") {
	  if (DecodeDouble(mUnitPeriod, inhalt)) {
	    cout  << "interne Zeit-Einheit " << mUnitPeriod  
		  << " s" << endl;
	  }
	  else return 0;
	}
	else if (name == "FrameLength") {
	  if (DecodeDouble(mFrameLength, inhalt)) {
	    cout  << "Laenge eines frames " << mFrameLength  
		  << " s" << endl;
	  }
	  else return 0;
	}
	else if (name == "DateTime") {
        DecodeDateTime(mDateTime, inhalt);
		cout << "Aufnahmedatum: " << inhalt << endl;
	} 
	else if (name == "GainOffset") {
	  if (DecodeDouble(mGainOffset, inhalt)) {
	    cout  << "Gain Offset " << mGainOffset   
	          <<  " dBm" << endl;
	  }
	  else return 0;
	}
	else if (name == "MultiFrames") {
	  if (DecodeInt(mMultiFrames, inhalt)) {
	    cout << mMultiFrames << " Frames pro Scan" << endl;
	  }
	  else return 0;
	}
	else if (name == "MultiAddr") {
	  if (DecodeInt(mMultiAddr, inhalt)) {
	    cout  << "Multi Addr " << mMultiAddr << endl;
	    int fehlende_frames = mMultiFrames - mMultiAddr - 1;
	    if (fehlende_frames)
	      cout << "Im letzten Scan fehlen " << fehlende_frames 
		   << " frames" << endl;
	  }
	  else return 0;
	}
	else if (name == "IOffset") {
	  if(DecodeDouble(mIOffset, inhalt)) {
	    cout << "Offset von I ist " << mIOffset << endl;
		}
	  else return 0;
	}
	else if (name == "QOffset") {
	  if(DecodeDouble(mQOffset, inhalt)) {
	    cout << "Offset von Q ist " << mIOffset << endl;
		}
	  else return 0;
	}
	return 1;
}

void IQData__TekRSA3303B::DecodeDateTime(char* ca, string s) {
    for (int i = 0; i < 10; i++) { // date
        if (s[i] >= '0' && s[i] <= '9')
            ca[i] = s[i];
        else
            ca[i] = '-';
    }
    ca[10] = ' ';
    for (int i = 11; i < 19; i++) { // time
        if ((s[i] >= '0' && s[i] <= '9') || s[i] == ':')
            ca[i] = s[i];
        else
            ca[i] = '\0';
    }
    ca[19] = '\0';
}
//
//               IQData__TekRSA3303B::Ganzzahl()
//
int IQData__TekRSA3303B::DecodeInt(int& result, string inhalt) {
	stringstream wandler(inhalt);
	wandler >> result;
	if (wandler.fail() | wandler.bad()) {
		cout << inhalt << " kann nicht in integer umgewandelt werden" << endl;
		return 0;
	}
	else return 1;
}
//
//               IQData__TekRSA3303B::DecodeDouble()
//
int IQData__TekRSA3303B::DecodeDouble(double& result, string inhalt) {
	double faktor = 1.;
	// letztes Zeichen herausholen
	int length = inhalt.length();
	string letztes_zeichen = inhalt.substr(length-1,1);
	//
	// ist das letzte Zeichen ein Buchstabe?
	// Wenn ja, daraus Faktor berechnen...
	int letztes_zeichen_ist_einheit = 1;
	if      (letztes_zeichen == "p") faktor = 1e-12;
	else if (letztes_zeichen == "n") faktor = 1e-9;
	else if (letztes_zeichen == "u") faktor = 1e-6;
	else if (letztes_zeichen == "m") faktor = 1e-3;
	else if (letztes_zeichen == "k") faktor = 1e3;
	else if (letztes_zeichen == "M") faktor = 1e6;
	else if (letztes_zeichen == "G") faktor = 1e9;
	else letztes_zeichen_ist_einheit = 0;
	string zahlenstring = inhalt.substr(0, length - letztes_zeichen_ist_einheit);
	stringstream wandler(zahlenstring);
	wandler >> result;
	if (wandler.fail() | wandler.bad()) {
		cout << zahlenstring << " kann nicht in double umgewandelt werden" << endl;
		cout << "letztes_zeichen_ist_einheit " << letztes_zeichen_ist_einheit << endl;
		cout << "length" << length << endl;
		return 0;
	}
	//
	// in jedem Fall den Rest nach double umrechnen
	//
	result *= faktor;
	return 1;
}
//
//               IQData__TekRSA3303B::GetFrameNumber()
//
int IQData__TekRSA3303B::GetFrameNumber(int& frame_nummer, int scan_nummer) {
        if (scan_nummer < 0 || scan_nummer > mValidScans) {
	  cout << "Scannummer " << scan_nummer 
	       << " ist ungueltig." << endl;
	  return 0;
	}
	frame_nummer = scan_nummer;
	return 1;
}
//
//               IQData__TekRSA3303B::GetDeltaT(double&);
//
int IQData__TekRSA3303B::GetDeltaT(double& delta_t) {
	if (!mReadCalled) {
		cout << "File ist noch nicht eingelesen" << endl;
		return 0;
	}
	unsigned long lt0 = mFrameHeaderPP[0]->LongTicks();
	unsigned long lt1 = mFrameHeaderPP[1]->LongTicks();
	delta_t = mUnitPeriod * (double) (lt1 - lt0);
	return 1;
}
//
//               Basisklasse von FrameHeader__TekRSA3303B und FrameData__TekRSA3303B
//
//
//               FramePart__TekRSA3303B::GetNextShort()
//
short FramePart__TekRSA3303B::GetNextShort() {
	union HiLo{
	char sc[2];
	short s;
	};
	HiLo hilo;
	mIQDataFile->get(hilo.sc[0]);
	mIQDataFile->get(hilo.sc[1]);
	return hilo.s;
}
//
//               FrameHeader__TekRSA3303B::Read()
//
void FrameHeader__TekRSA3303B::Read() {
	short dummy;
	dummy = GetNextShort();
	mValidA = GetNextShort();
	mValidP = GetNextShort();
	mValidI = GetNextShort();
	mValidQ = GetNextShort();
	mBins = GetNextShort();
	dummy = GetNextShort();
	mTriggered = GetNextShort();
	mOverLoad = GetNextShort();
	mLastFrame = GetNextShort();
	union UnsignedLong{
	     char ulc[sizeof(unsigned int)]; // changed for compatibility to 64 bit machines
		unsigned int ul;
	};
	UnsignedLong ul;
	unsigned int ib;
	for (ib=0; ib < 4; ib++) 
		mIQDataFile->get(ul.ulc[ib]);
	mLongTicks = ul.ul;
}
//
//               FrameHeader__TekRSA3303B::PrintAll()
//
void FrameHeader__TekRSA3303B::PrintAll() {
	cout << "ValidA : "     << mValidA     << endl;
	cout << "ValidP : "     << mValidP     << endl;
	cout << "ValidI : "     << mValidI     << endl;
	cout << "ValidQ : "     << mValidQ     << endl;
	cout << "Bins : "       << mBins       << endl;
	cout << "Triggered : "  << mTriggered  << endl;
	cout << "OverLoad : "   << mOverLoad   << endl;
	cout << "LastFrame : "  << mLastFrame  << endl;
	cout << "LongTicks : "  << mLongTicks  << endl;
}
//
//               FrameData__TekRSA3303B::FrameData__TekRSA3303B()
//
FrameData__TekRSA3303B::FrameData__TekRSA3303B(IQData__TekRSA3303B * iq_data_p)
	: FramePart__TekRSA3303B(iq_data_p) {
	IVectorP = new short[1024];
	QVectorP = new short[1024];
}
//
//               FrameData__TekRSA3303B::~FrameData__TekRSA3303B()
//
FrameData__TekRSA3303B::~FrameData__TekRSA3303B() {
	delete [] IVectorP;
	delete [] QVectorP;
}
//
//               FrameData__TekRSA3303B::Read()
//
void FrameData__TekRSA3303B::Read() {
	int i;
	short * i_p = IVectorP;
	short * q_p = QVectorP;
	for (i = 0; i < 1024; i++) {
	*q_p++ = GetNextShort();
        *i_p++ = GetNextShort(); 
	}
}
//
//               ExtCorrData__TekRSA3303B::FrameData__TekRSA3303B()
//
ExtCorrData__TekRSA3303B::ExtCorrData__TekRSA3303B(IQData__TekRSA3303B * iq_data_p)
	: FramePart__TekRSA3303B(iq_data_p) {
	AVectorP = new char[1024];
	PVectorP = new char[1024];
}
//
//               ExtCorrData__TekRSA3303B::~FrameData__TekRSA3303B()
//
ExtCorrData__TekRSA3303B::~ExtCorrData__TekRSA3303B() {
	delete [] AVectorP;
	delete [] PVectorP;
}
//
//               ExtCorrData__TekRSA3303B::Read()
//
void ExtCorrData__TekRSA3303B::Read() {
	int i;
	char * a_p = AVectorP;
	char * p_p = PVectorP;
	for (i = 0; i < 1024; i++) 
		*a_p++ = mIQDataFile->get();
	for (i = 0; i < 1024; i++) 
		*p_p++ = mIQDataFile->get();
}
