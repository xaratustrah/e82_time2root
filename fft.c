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


/*------------------------------------------------------------------*
 | This program is aimed to provide discrete Fourier transform,     |
 | based on decimation-in-time fast Fourier transform algorithm.    |
 | The number of points of DFT must be an integer power of 2, and   |
 | no greater than 2^15.                                            |
 *------------------------------------------------------------------*/


#include <TComplex.h>
#include <TDatime.h>

#include "tiq.h"

#define cMaxFFTPt   2048 /* no greater than 2^15 */

/* index vector for bit-reversal rearrangement */
static unsigned short Idx[cMaxFFTPt] = {};
/* initialize an index vector for the purpose of bit reversing */
void SetIdx(unsigned short);

/* Seq is the sequence to be transformed, |Num| is the dimension of the
   sequence (negative Num indicates an inverse FFT is applied). */
bool FFT(TComplex Seq[], short Num)
{
    bool IsInv;
    unsigned char Exp, Step;
    unsigned short Block, Base, Offset, Elem;

    if (Num == 0)
    {
        fprintf(stderr, "Error: FFT point number can't be zero.\n");
        return false;
    }
    else if (Num > 0)
        IsInv = false;
    else /* inverse FFT will be applied */
    {
        IsInv = true;
        Num = -Num;
    }

    if (Num > cMaxFFTPt)
    {
        fprintf(stderr, "Error: FFT point number %d exceed limit %d.\n",
                Num, cMaxFFTPt);
        return false;
    }

    Exp = (unsigned char) TMath::Log2(Num);
    if (Num != TMath::Power(2, Exp))
    {
        fprintf(stderr, "Error: FFT point number "
                "%d must be an integer power of 2.\n", Num);
        return false;
    }

    /* build index vector */
    TComplex temp;
    if (Idx[1] != Num / 2)
        SetIdx(Num);

    /* rearrange input array to bit-reversed order */
    for (unsigned short i = 0; i < Num; i++)
    {
        if (i < Idx[i])
        {
            temp = Seq[Idx[i]];
            Seq[Idx[i]] = Seq[i];
            Seq[i] = temp;
        }
    }

    /* apply fast Fourier transform */
    TComplex W_N = TComplex(1, (IsInv ? 1 : -1) * TMath::TwoPi() / Num, true);
    for (Step = 0; Step < Exp; Step++)
    {
        Block = (unsigned short) TMath::Power(2, Step + 1);
        for (Base = 0; Base < Num; Base += Block)
        {
            Offset = Block / 2;
            for (Elem = 0; Elem < Offset; Elem++)
            {
                temp  = TComplex::Power(W_N, Elem * (unsigned short) 
                        TMath::Power(2, Exp - Step - 1));
                temp *= Seq[Base + Elem + Offset];
                Seq[Base + Elem + Offset] = Seq[Base + Elem] - temp;
                Seq[Base + Elem] += temp;
                if (IsInv && Step == Exp - 1)
                {
                    Seq[Base + Elem] /= Num;
                    Seq[Base + Elem + Offset] /= Num;
                }
            }
        }
    }
    return true;
}

/* initialize index vector, based on bit-reversal algorithm developed by
   M. Rubio, et al. [Int. J. Adapt. Control Signal Process, 16 (2002) 703] */
void SetIdx(unsigned short Num)
{
    unsigned short lower, upper;
    unsigned char Exp = (unsigned char) TMath::Log2(Num);
    Idx[1] = 1;
    Idx[1] <<= Exp - 1;

    for (unsigned char i = 2; i <= Exp; i++)
    {
        lower = (unsigned short) TMath::Power(2, i-1) - 1;
        upper = (unsigned short) TMath::Power(2, i) - 1;
        Idx[upper] = 1;
        Idx[upper] <<= Exp - i;
        Idx[upper] ^= Idx[lower];
        for (unsigned short j = 0; j < lower; j++)
            Idx[lower+j+1] = Idx[upper] ^ Idx[lower-j];
    }
}
