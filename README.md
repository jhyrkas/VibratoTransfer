This plug-in implements a real-time algorithm for vibrato transfer.

Right now, the plug-in takes stereo input with the assumption that the input sound (the audio on which to
apply vibrato) is in the left channel and the target sound (the audio that contains vibrato) is in the right
channel. In the future, the plug-in will handle stereo input and take the target audio in a sidechain.

The plug-in contains two sliders that control the amount of amplitude and frequency modulation that will be
applied from the target signal. When the slider for each parameter is set to 0, that modulation type will not
be applied at all. When the target signal is quiet or has an unstable fundamental frequency, no vibrato will
be applied to the input. 

The plug-in operates with a base latency of 512 samples between the input and output
signal even when no vibrato is being applied. Latency can increase or decreaes as vibrato is applied due to
the vibrato patterns in the target signal.

The algorithms in this plug-in are described in a previously published and upcoming research paper. The plug-in
implements a real-time vibrato transfer but is not optimized like a commercial plug-in. 
Original code in this repo uses the MIT License (see language below). 
**However**, two open-source implementations are used: Ron Mayer's FFT implementation, found in Pure Data and elsewhere, and
Edson Niu's implementation of a Butterworth bandpass filter which is based on the butter() function in MATLAB.
Sources for these implementations are cited in the code and these implementations carry their own licenses.

MIT License

Copyright (c) [2024] [authorname]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
