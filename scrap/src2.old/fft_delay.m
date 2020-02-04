function outSig = fft_delay(inSig, td, fs)
    %FFT_DELAY delays inSig by td with information fs
    slen = length(inSig);
    nfft = 2^nextpow2(2 * slen); % To use max. computational efficiency of FF
    fax = fs * (-nfft / 2:nfft / 2 - 1)' / nfft; %  Create frequency shift vectors(bins)

    shft = exp(-1i * td * 2 * pi * fax); % w=2*pi*fax,, Frequency function for delay
    shft = ifftshift(shft); % Make axis compatable with numeric FFT
    fsd = fft(inSig, nfft); % Take FFT
    fsd = fsd .* shft; %  Apply delay
    dum = ifft(fsd); %  Return to time domain
    outSig = dum(1:slen); %  Trim time domain signal to required length
end
