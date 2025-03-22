This plug-in implements a real-time algorithm for vibrato transfer.
The algorithm and impementation are detailed in an upcoming paper titled "Real-time implementation of vibrato transfer as an audio effect," which will be presented at ICMC 2025.

To build the plug-in, you will need to use download the JUCE framework and open VibratoTransfer.jucer file using the Projucer software.
This code is tested on macOS and has not yet been tested on Windows.

The plug-in takes an input sound in stereo (the audio on which to apply vibrato) and the target sound (the audio that contains vibrato) in the sidechain.
The plug-in contains three sliders that control the amount of amplitude and frequency modulation that will be
applied from the target signal and make-up gain. 
When the target signal is quiet or has an unstable fundamental frequency, no vibrato will be applied to the input. 

The plug-in operates with a base latency of 512 samples between the input and output
signal even when no vibrato is being applied. Latency can increase or decrease as vibrato is applied due to
the vibrato patterns in the target signal.

The plug-in implements a real-time vibrato transfer but is not optimized like a commercial plug-in. 
Two open-source implementations are used: Ron Mayer's FFT implementation, found in Pure Data and elsewhere, and
Ruoho Ruotsi's implementation of a Butterworth bandpass filter which is based on the butter() function in MATLAB.
Sources for these implementations are cited in the code and these implementations carry their own licenses.

Audio examples created from an early implementation can be viewed here: https://youtu.be/cW_WZEZrTmA
