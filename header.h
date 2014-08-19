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



#ifndef _HEADER_H_
#define _HEADER_H_

#include <iostream>
#include <TNamed.h>

using namespace std;

class Header : public TNamed {
    public:
        Header() :
            fValidFrames(0),
            fFrameLength(.0),
            fCenterFrequency(.0),
            fSpan(.0),
            fGainOffset(.0),
            fScaling(1.) { SetDateTime("1970-01-01 00:00:00"); }
        Header(const char* name, const char* title) : TNamed(name, title),
            fValidFrames(0),
            fFrameLength(.0),
            fCenterFrequency(.0),
            fSpan(.0),
            fGainOffset(.0),
            fScaling(1.) { SetDateTime("1970-01-01 00:00:00"); }
        Header(const Header&);
        Header& operator=(const Header&);
        virtual ~Header() { }

        void SetValidFrames(int n) { fValidFrames = n; }
        void SetFrameLength(double x) { fFrameLength = x; }
        void SetCenterFrequency(double x) { fCenterFrequency = x; }
        void SetSpan(double x) { fSpan = x; }
        void SetGainOffset(double x) { fGainOffset = x; }
        void SetScaling(double x) { fScaling = x; }
        void SetDateTime(const char*);

        int GetValidFrames() const { return fValidFrames; }
        double GetFrameLength() const { return fFrameLength; }
        double GetCenterFrequency() const { return fCenterFrequency; }
        double GetSpan() const { return fSpan; }
        double GetGainOffset() const { return fGainOffset; }
        double GetScaling() const { return fScaling; }
        const char* GetDateTime() const { return fDateTime; }

        void Show() const {
            ShowValidFrames();
            ShowFrameLength();
            ShowCenterFrequency();
            ShowSpan();
            ShowGainOffset();
            ShowScaling();
            ShowDateTime();
        }

    protected:
        int    fValidFrames;
        double fFrameLength; // s
        double fCenterFrequency; // Hz
        double fSpan; // Hz
        double fGainOffset; // dB, for IQT
        double fScaling; // lin, for TIQ
        char   fDateTime[20]; // yyyy-mm-dd hh:mm:ss

        void ShowValidFrames() const;
        void ShowFrameLength() const;
        void ShowCenterFrequency() const;
        void ShowSpan() const;
        void ShowGainOffset() const;
        void ShowScaling() const;
        void ShowDateTime() const;

        ClassDef(Header, 1)
};

#endif // _HEADER_H_
