% DCS_Assignment_1920
% Rui Oliveira
% Simulations for the AWGN channel

close all;
clear all;

%--- Parameters ---
% Number of tx bits
N = 2^18;

% SNRs to try
minsnr = -20;
maxsnr = 20;
snrstep = 5;
tsnrs = minsnr:snrstep:maxsnr;

% Modulation orders for PSK
dopsk = true;
pskMs = [2, 4];
% Modulation orders for QAM
doqam = true;
qamMs = [4, 64];
% Parameters for FSK
dofsk = true;
fskMs = [2];
freqsep = 8; % Frequency separation (Hz)
nsamp = 16; % Number of samples per symbol
Fs = freqsep * max(fskMs) * 4; % Sample rate (Hz)

% CRC codings to try
crcCodings = [
        8 4;
        12 7; ];

%--- Simulation ---
% Generate tx stream
bitStream = randi([0 1], N, 1);

% Open an array for the results
resultspskcrc = zeros(4, length(crcCodings), length(pskMs), length(tsnrs));
resultsqamcrc = zeros(4, length(crcCodings), length(qamMs), length(tsnrs));
resultsfskcrc = zeros(4, length(crcCodings), length(fskMs), length(tsnrs));

% CRC FEC
% Since generating the codes takes ages,
% it's better to do the iteration in "reverse"
for c = 1:length(crcCodings)
    % Acquire the code and message lengths
    clen = crcCodings(c, 1);
    mlen = crcCodings(c, 2);

    % Do the FEC encoding
    % First we need to adjust the bitstream so that it fits
    % We need to add mlen-mod(len,mlen) bits
    bitStreamFEC = [bitStream; randi([0 1], mlen - mod(length(bitStream), mlen), 1)];

    % Now the message can fit into the code
    % Create a generator polynomial for a cyclic code.
    gpol = cyclpoly(clen, mlen);
    % Create a parity-check matrix by using the generator polynomial.
    parmat = cyclgen(clen, gpol);
    % Create a syndrome decoding table by using the parity-check matrix.
    trt = syndtable(parmat);

    % Encode the data by using the generator polynomial.
    encData = encode(bitStreamFEC, clen, mlen, 'cyclic/binary', gpol);

    % Now we can transmit it in the various conditions

    for a = 1:length(tsnrs)
        tsnr = tsnrs(a);

        % PSK
        if dopsk

            for b = 1:length(pskMs)
                M = pskMs(b);
                nbits = log2(M);

                fomatSpec = "CRC = %d %d; SNR = %d; PSK, M = %d\n";
                fprintf(fomatSpec, clen, mlen, tsnr, M);

                % First, fit the array into the symbols
                % We need to add nbits-mod(len,nbits) bits
                bitStreamFixed = [encData; randi([0 1], nbits - mod(length(encData), nbits), 1)];

                % Now that the data fits, we transmit it

                % To do that, we reshape the stream to convert into symbols
                databin = reshape(bitStreamFixed, [length(bitStreamFixed) / nbits, nbits]);
                % Get symbols from the binary stream
                data = bi2de(databin);

                % Modulate data
                if M == 2
                    % For M = 2, initial phase is more convenient zero
                    txSig = pskmod(data, M, 0, 'gray');
                else
                    % Else, as per documentation
                    txSig = pskmod(data, M, pi / M, 'gray');
                end

                % Appy AWGN
                rxSig = awgn(txSig, tsnr, 'measured');

                % Demodulate data
                if M == 2
                    % For M = 2, initial phase is more convenient zero
                    dataOut = pskdemod(rxSig, M, 0, 'gray');
                else
                    % Else, as per documentation
                    dataOut = pskdemod(rxSig, M, pi / M, 'gray');
                end

                % Convert data back into binary
                dataRxBin = de2bi(dataOut, log2(M));
                % Convert back into line vector
                dataRxBin = reshape(dataRxBin, [], 1);

                % Now retreive the data
                % Adjust the size to the message length (aka, remove what we added)
                dataRxBin = dataRxBin(1:length(encData));
                % Decode the corrupted sequence
                dataRecBinary = decode(dataRxBin, clen, mlen, 'cyclic/binary', gpol, trt);
                % Cut the recovered binary into the original shape
                dataRecBinary = dataRecBinary(1:length(bitStream));

                % Now check what we got
                [~, ber] = biterr(bitStream, dataRecBinary);
                [~, ser] = symerr(bitStream, dataRecBinary);

                % And compute the theoretical value
                [berteo, serteo] = berawgn(tsnr - 10 * log10(log2(M)), "psk", M, 'nondiff');

                resultspskcrc(1, c, b, a) = ber;
                resultspskcrc(2, c, b, a) = ser;
                resultspskcrc(3, c, b, a) = berteo;
                resultspskcrc(4, c, b, a) = serteo;
            end % pskMs

        end % dopsk

        % QAM
        if doqam

            for b = 1:length(qamMs)
                M = qamMs(b);
                nbits = log2(M);

                fomatSpec = "CRC = %d %d; SNR = %d; QAM, M = %d\n";
                fprintf(fomatSpec, clen, mlen, tsnr, M);

                % First, fit the array into the symbols
                % We need to add nbits-mod(len,nbits) bits
                bitStreamFixed = [encData; randi([0 1], nbits - mod(length(encData), nbits), 1)];

                % Now that the data fits, we transmit it

                % To do that, we reshape the stream to convert into symbols
                databin = reshape(bitStreamFixed, [length(bitStreamFixed) / nbits, nbits]);
                % Get symbols from the binary stream
                data = bi2de(databin);

                % Modulate data
                txSig = qammod(data, M, 'gray');

                % Appy AWGN
                rxSig = awgn(txSig, tsnr, 'measured');

                % Demodulate data
                dataOut = qamdemod(rxSig, M, 'gray');

                % Convert data back into binary
                dataRxBin = de2bi(dataOut, log2(M));
                % Convert back into line vector
                dataRxBin = reshape(dataRxBin, [], 1);

                % Now retreive the data
                % Adjust the size to the message length (aka, remove what we added)
                dataRxBin = dataRxBin(1:length(encData));
                % Decode the corrupted sequence
                dataRecBinary = decode(dataRxBin, clen, mlen, 'cyclic/binary', gpol, trt);
                % Cut the recovered binary into the original shape
                dataRecBinary = dataRecBinary(1:length(bitStream));

                % Now check what we got
                [~, ber] = biterr(bitStream, dataRecBinary);
                [~, ser] = symerr(bitStream, dataRecBinary);

                % And compute the theoretical value
                [berteo, serteo] = berawgn(tsnr - 10 * log10(log2(M)), "qam", M);

                resultsqamcrc(1, c, b, a) = ber;
                resultsqamcrc(2, c, b, a) = ser;
                resultsqamcrc(3, c, b, a) = berteo;
                resultsqamcrc(4, c, b, a) = serteo;
            end % qamMs

        end % doqam

        % FSK
        if dofsk

            for b = 1:length(fskMs)
                M = fskMs(b);
                nbits = log2(M);

                fomatSpec = "CRC = %d %d; SNR = %d; FSK, M = %d\n";
                fprintf(fomatSpec, clen, mlen, tsnr, M);

                % First, fit the array into the symbols
                % We need to add nbits-mod(len,nbits) bits
                bitStreamFixed = [encData; randi([0 1], nbits - mod(length(encData), nbits), 1)];

                % Now that the data fits, we transmit it

                % To do that, we reshape the stream to convert into symbols
                databin = reshape(bitStreamFixed, [length(bitStreamFixed) / nbits, nbits]);
                % Get symbols from the binary stream
                data = bi2de(databin);

                % Modulate data
                txSig = fskmod(data, M, freqsep, nsamp, Fs, "cont", 'gray');

                % Appy AWGN
                rxSig = awgn(txSig, tsnr, 'measured');

                % Demodulate data
                dataOut = fskdemod(rxSig, M, freqsep, nsamp, Fs, "gray");

                % Convert data back into binary
                dataRxBin = de2bi(dataOut, log2(M));
                % Convert back into line vector
                dataRxBin = reshape(dataRxBin, [], 1);

                % Now retreive the data
                % Adjust the size to the message length (aka, remove what we added)
                dataRxBin = dataRxBin(1:length(encData));
                % Decode the corrupted sequence
                dataRecBinary = decode(dataRxBin, clen, mlen, 'cyclic/binary', gpol, trt);
                % Cut the recovered binary into the original shape
                dataRecBinary = dataRecBinary(1:length(bitStream));

                % Now check what we got
                [~, ber] = biterr(bitStream, dataRecBinary);
                [~, ser] = symerr(bitStream, dataRecBinary);

                % And compute the theoretical value
                % https://www.mathworks.com/help/comm/ug/awgn-channel.html#a1071501088
                [berteo, serteo] = berawgn(tsnr - 10 * log10(log2(M)) + 10 * log10(nsamp), "fsk", M, 'noncoherent');

                resultsfskcrc(1, c, b, a) = ber;
                resultsfskcrc(2, c, b, a) = ser;
                resultsfskcrc(3, c, b, a) = berteo;
                resultsfskcrc(4, c, b, a) = serteo;
            end % fskMs

        end % dofsk

    end % tsnrs

end % crcCodings

%% Plots
figure(1); hold on;
title("BER w/ CRC vs SNR");

% Experimental data
for c = 1:length(crcCodings)
    % Code and message lengths
    clen = crcCodings(c, 1);
    mlen = crcCodings(c, 2);

    if dopsk

        for b = 1:length(pskMs)
            plot(tsnrs, max(eps, reshape(resultspskcrc(1, c, b, :), [], 1)), '--*', ...
                'DisplayName', sprintf("PSK, M=%d; CRC, k=%d, n=%d", pskMs(b), clen, mlen));
        end

    end

    if doqam

        for b = 1:length(qamMs)
            plot(tsnrs, max(eps, reshape(resultsqamcrc(1, c, b, :), [], 1)), '--*', ...
                'DisplayName', sprintf("QAM, M=%d; CRC, k=%d, n=%d", qamMs(b), clen, mlen));
        end

    end

    if dofsk

        for b = 1:length(fskMs)
            plot(tsnrs, max(eps, reshape(resultsfskcrc(1, c, b, :), [], 1)), '--*', ...
                'DisplayName', sprintf("FSK, M=%d; CRC, k=%d, n=%d", fskMs(b), clen, mlen));
        end

    end

end

% Theoretical lines
if dopsk

    for b = 1:length(pskMs)
        plot(tsnrs, max(eps, reshape(resultspskcrc(3, 1, b, :), [], 1)), ...
            'DisplayName', sprintf("PSK, M=%d no FEC", pskMs(b)));
    end

end

if doqam

    for b = 1:length(qamMs)
        plot(tsnrs, max(eps, reshape(resultsqamcrc(3, 1, b, :), [], 1)), ...
            'DisplayName', sprintf("QAM, M=%d no FEC", qamMs(b)));
    end

end

if dofsk

    for b = 1:length(fskMs)
        plot(tsnrs, max(eps, reshape(resultsfskcrc(3, 1, b, :), [], 1)), ...
            'DisplayName', sprintf("FSK, M=%d no FEC", fskMs(b)));
    end

end

legend;
legend('Location', 'southwest', 'NumColumns', 1);
grid("minor");
hold off;
