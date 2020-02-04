h = spectrum.welch('Rectangular');
figure(1)
pwelch(x,[],[],[],fs);
figure(2)
pwelch(y,[],[],[],fs);
