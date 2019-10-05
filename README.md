#### A 4FSK Codec2 Digital Voice Modem Library in Java 8 SE
This is an experimental modem object for transmitting digital voice in 4FSK modulation.

This is a Java translation of the Codec2 Mode 700C Vocoder, and the 4FSK modem of the FreeDV project, and combined into a Library Object to be linked in to your GUI and audio interface threads.

Also included is a test program to exercise the modem. The data is 800 bit/s, with a symbol rate of 400 Hz, and the sample rate is 8 kHz (little-endian, mono, Signed 16-bit PCM).

The results is a spectrum with four FSK frequencies spaced by the symbol rate. These being: 800 Hz, 1200 Hz, 1600 Hz, and 2000 Hz.

A sample audio file is included that you can analyze (using Audacity for example).
