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
 | Reader for Header	      	           |
 *-----------------------------------------*/


#include <TString.h>
#include <TDatime.h>
#include <cstdio>
#include "tiq.h"

/* get the specified attribute of an element from xml header */
void    GetAttrib(FILE*, TString*, int);
/* get the value of an element from xml header */
void    GetValue(FILE*, TString*);
/* get the TDatime from a TString */
TDatime GetDaTm(TString*);

bool SetInfo(FILE *fp, Info_t *pInfo) {
    char c;
    int countdown = 8;
    int nsmpl, nfrm;
    TString key;
    TString* info = new TString();

    if (fgetc(fp) != '<') {
        fprintf(stderr, "Error: wrong file type!\n");
        return false;
    }

    rewind(fp);
    while ((c = fgetc(fp)) != EOF && countdown) {
        if (c != '<') continue;
        else {
            key.Resize(0);
            while (isalpha(c = fgetc(fp)))
                key += c;

            if (!key.Length()) continue;
            else if (key.EqualTo("DataFile")) {
                GetAttrib(fp, info, 1);
                pInfo->Offset = info->Atoi();
                countdown--;
            } else if (key.EqualTo("DataSets")) {
                GetAttrib(fp, info, 2);
                nfrm = info->Atoi();
                countdown--;
            } else if (key.EqualTo("SamplingFrequency")) {
                GetValue(fp, info);
                pInfo->Intvl = 1. / info->Atof();
                countdown--;
            } else if (key.EqualTo("NumberSamples")) {
                GetValue(fp, info);
                nsmpl = info->Atoi();
                countdown--;
            } else if (key.EqualTo("Scaling")) {
                GetValue(fp, info);
                //pInfo->Scaling = sqrt(info->Atof());
                pInfo->Scaling = info->Atof();
                countdown--;
            } else if (key.EqualTo("DateTime")) {
                GetValue(fp, info);
                pInfo->DaTm = GetDaTm(info);
                countdown--;
            } else if (key.EqualTo("Frequency")) {
                GetValue(fp, info);
                pInfo->CenFreq = info->Atof();
                countdown--;
            } else if (key.EqualTo("AcquisitionBandwidth")) {
                GetValue(fp, info);
                pInfo->Span = info->Atof();
                countdown--;
            }
        }
    }
    delete info;

    if (countdown) {
        fprintf(stderr, "Error: no enough information in the header!\n");
        return false;
    } else {
        pInfo->NumPt = nsmpl * nfrm;
        return true;
    }
}

void GetAttrib(FILE* fp, TString* s, int natt) {
    char c;
    int nquo = 2 * natt;
    s->Resize(0);
    while (nquo) {
        c = fgetc(fp);
        if (c == '"') nquo--;
        if (nquo == 1 && c != '"') *s += c;
    }
}

void GetValue(FILE* fp, TString* s) {
    char c;
    s->Resize(0);
    fseek(fp, -1, SEEK_CUR);
    while (fgetc(fp) != '>') ;
    while ((c = fgetc(fp)) != '<') *s += c;
}

TDatime GetDaTm(TString* s) {
    s->Replace(10, 1, " ");
    TString t((*s)(0, 19));
    return TDatime(t.Data());
}
