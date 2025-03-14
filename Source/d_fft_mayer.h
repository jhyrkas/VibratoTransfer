/*
 * d_fft_mayer.h
 */

/*
 This code is largely adapted from an older FFT implementation in Pure Data.
 The original FFT is implemented by Ron Mayer. This implementation may be protected
 by the BSD license, which copied below from the Pure Data official Githb repository:

 This software is copyrighted by Miller Puckette and others.  The following
 terms (the "Standard Improved BSD License") apply to all files associated with
 the software unless explicitly disclaimed in individual files:

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are
 met:

 1. Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.
 2. Redistributions in binary form must reproduce the above
    copyright notice, this list of conditions and the following
    disclaimer in the documentation and/or other materials provided
    with the distribution.
 3. The name of the author may not be used to endorse or promote
    products derived from this software without specific prior
    written permission.

 THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY
 EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR
 BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
 IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef D_FFT_MAYER_H_
#define D_FFT_MAYER_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef float t_sample;       /* a float type at most the same size */

#define REAL t_sample

void mayer_fht(REAL *fz, int n, REAL *coswrk, REAL *sinwrk);
void mayer_fft(int n, REAL *real, REAL *imag, REAL *coswrk, REAL *sinwrk);
void mayer_ifft(int n, REAL *real, REAL *imag, REAL *coswrk, REAL *sinwrk);
void mayer_realfft(int n, REAL *real, REAL *coswrk, REAL *sinwrk);
void mayer_realifft(int n, REAL *real, REAL *coswrk, REAL *sinwrk);


#ifdef __cplusplus
}
#endif

#endif /* D_FFT_MAYER_H_ */
