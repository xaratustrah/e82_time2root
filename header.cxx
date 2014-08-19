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



#include "header.h"

ClassImp(Header)

Header::Header(const Header& h) : TNamed(h) {
    fValidFrames = h.fValidFrames;
    fFrameLength = h.fFrameLength; 
    fCenterFrequency = h.fCenterFrequency;
    fSpan = h.fSpan;
    fGainOffset = h.fGainOffset;
    fScaling = h.fScaling;
    SetDateTime(h.GetDateTime());
}

Header& Header::operator=(const Header& rhs) {
    if (this != &rhs) {
        TNamed::operator=(rhs);
        fValidFrames = rhs.fValidFrames;
        fFrameLength = rhs.fFrameLength; 
        fCenterFrequency = rhs.fCenterFrequency;
        fSpan = rhs.fSpan;
        fGainOffset = rhs.fGainOffset;
        fScaling = rhs.fScaling;
        SetDateTime(rhs.GetDateTime());
    }
    return *this;
}

void Header::SetDateTime(const char* s) {
    for (int i = 0; i < 19; i++)
        fDateTime[i] = s[i];
    fDateTime[19] = '\0';
}

void Header::ShowValidFrames() const {
    cout << "Valid Frames: " << fValidFrames << endl;
}

void Header::ShowFrameLength() const {
    cout << "Frame Length: " << fFrameLength * 1e3 << " ms" << endl;
}

void Header::ShowCenterFrequency() const {
    cout << "Center Frequency: " << fCenterFrequency / 1e6 << " MHz" << endl;
}

void Header::ShowSpan() const {
    cout << "Span: " << fSpan / 1e3 << " kHz" << endl;
}

void Header::ShowGainOffset() const {
    if (fGainOffset == .0)
        return;
    cout << "Gain Offset: " << fGainOffset << " dB" << endl;
}

void Header::ShowScaling() const {
    if (fScaling == 1.)
        return;
    cout << "Scaling: " << fScaling << endl;
}

void Header::ShowDateTime() const {
    cout << "Date and Time: " << fDateTime << endl;
}
