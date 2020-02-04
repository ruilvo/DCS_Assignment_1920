close all;
clear all;

fs = 1e6;

fc = 1e3;

fD = 500;

npoints = 1000;

tvec = (1:npoints).' ./ fs;

xt = sin(2 * pi * fc * tvec);

pathDelays = [0, 10] * 1e-9;
avgPathGains = [0, -10];

rchan = comm.RayleighChannel(...
    'SampleRate', fs, ...
    'PathDelays', pathDelays, ...
    'AveragePathGains', avgPathGains, ...
    "FadingTechnique", "Sum of sinusoids", ...
    "MaximumDopplerShift", 0 ...
    );

yt = rchan(xt);

figure(1); hold on;
plot(tvec, xt);
plot(tvec, real(yt));
