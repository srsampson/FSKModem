#### A 4FSK Codec2 Digital Voice Modem Library in Java 8 SE
This is an experimental modem object for transmitting digital voice in 4FSK modulation.

This is a Java translation of the Codec2 Mode 700C Vocoder, and the 4FSK modem of the FreeDV project, and combined into a Library Object to be linked in to your GUI and audio interface threads.

##### FSKTest Java Program
Also included is a test program to exercise the modem. The data is 800 bit/s, with a symbol rate of 400 baud, and the sample rate is 8 kHz (little-endian, mono, Signed 16-bit PCM).

The results is a spectrum with four FSK frequencies spaced by the symbol rate. These being: 800 Hz, 1200 Hz, 1600 Hz, and 2000 Hz.

A sample audio file is included that you can analyze (using Audacity for example).

I need to write a better test program, because you really can't just send individual frames and get a good output. You need to send a contiguous stream of frames, and then the output makes sense bit-wise.

In order to run the example, it uses the ```FSKModem.jar``` as a library. Thus, you need a ```lib``` directory and the jar file placed in it. For example:
```
mkdir ~/test/lib
cp FSKModem.jar ~/test/lib
cp FSKTest.jar ~/test
cd ~/test
java -jar FSKTest.jar
```

