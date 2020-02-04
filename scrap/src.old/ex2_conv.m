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
snrstep = 2;
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

% Conv encodings to try
% Must generate compatible arguments
% Two arguments = no feedback
% Three argument = feedback
trellisArgs = {
{[5 4], [23 35 0; 0 5 13]}, % 2/3 Feedforward Convolutional Encoder
{7, {'1 + x^3 + x^4 + x^5 + x^6', '1 + x + x^3 + x^4 + x^6'}}, % 1/2 Feedforward Convolutional Encoder
{5, [37 33], 37}}; % 1/2 Feedback Convolutional Encoders

%--- Simulation ---
% Generate tx stream
bitStream = randi([0 1], N, 1);

% Open an array for the results
resultspskconv = zeros(4, length(trellisArgs), length(pskMs), length(tsnrs));
resultsqamconv = zeros(4, length(trellisArgs), length(qamMs), length(tsnrs));
resultsfskconv = zeros(4, length(trellisArgs), length(fskMs), length(tsnrs));

for c = 1:length(trellisArgs)
    currArgs = trellisArgs{c};

    % Convert convolutional code polynomials to trellis description
    switch length(currArgs)
        case 2
            trellis = poly2trellis(currArgs{1}, currArgs{2});
        case 3
            trellis = poly2trellis(currArgs{1}, currArgs{2}, currArgs{3});
    end

    % Convolutionally encode the data
    encData = convenc(bitStream, trellis);

    % Now we can transmit it in the various conditions
    for a = 1:length(tsnrs)
        tsnr = tsnrs(a);

        % PSK
        if dopsk

            for b = 1:length(pskMs)
                M = pskMs(b);
                nbits = log2(M);

                fomatSpec = "Conv: %d / %d; SNR = %d; PSK, M = %d\n";
                fprintf(fomatSpec, c, length(trellisArgs), tsnr, M);

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
                % Decode the corrupted sequence using the Viterbi algorithm
                dataRecBinary = vitdec(dataRxBin, trellis, 34, 'trunc', 'hard');

                % Now check what we got
                [~, ber] = biterr(bitStream, dataRecBinary);
                [~, ser] = symerr(bitStream, dataRecBinary);

                % And compute the theoretical value
                [berteo, serteo] = berawgn(tsnr - 10 * log10(log2(M)), "psk", M, 'nondiff');

                resultspskconv(1, c, b, a) = ber;
                resultspskconv(2, c, b, a) = ser;
                resultspskconv(3, c, b, a) = berteo;
                resultspskconv(4, c, b, a) = serteo;
            end % pskMs

        end % dopsk

        % QAM
        if doqam

            for b = 1:length(qamMs)
                M = qamMs(b);
                nbits = log2(M);

                fomatSpec = "Conv: %d / %d; SNR = %d; QAM, M = %d\n";
                fprintf(fomatSpec, c, length(trellisArgs), tsnr, M);

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
                % Decode the corrupted sequence using the Viterbi algorithm
                dataRecBinary = vitdec(dataRxBin, trellis, 34, 'trunc', 'hard');

                % Now check what we got
                [~, ber] = biterr(bitStream, dataRecBinary);
                [~, ser] = symerr(bitStream, dataRecBinary);

                % And compute the theoretical value
                [berteo, serteo] = berawgn(tsnr - 10 * log10(log2(M)), "qam", M);

                resultsqamconv(1, c, b, a) = ber;
                resultsqamconv(2, c, b, a) = ser;
                resultsqamconv(3, c, b, a) = berteo;
                resultsqamconv(4, c, b, a) = serteo;
            end % qamMs

        end % doqam

        % FSK
        if dofsk

            for b = 1:length(fskMs)
                M = fskMs(b);
                nbits = log2(M);

                fomatSpec = "Conv: %d / %d; SNR = %d;FSK, M = %d\n";
                fprintf(fomatSpec, c, length(trellisArgs), tsnr, M);

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
                % Decode the corrupted sequence using the Viterbi algorithm
                dataRecBinary = vitdec(dataRxBin, trellis, 34, 'trunc', 'hard');

                % Now check what we got
                [~, ber] = biterr(bitStream, dataRecBinary);
                [~, ser] = symerr(bitStream, dataRecBinary);

                % And compute the theoretical value
                % https://www.mathworks.com/help/comm/ug/awgn-channel.html#a1071501088
                [berteo, serteo] = berawgn(tsnr - 10 * log10(log2(M)) + 10 * log10(nsamp), "fsk", M, 'noncoherent');

                resultsfskconv(1, c, b, a) = ber;
                resultsfskconv(2, c, b, a) = ser;
                resultsfskconv(3, c, b, a) = berteo;
                resultsfskconv(4, c, b, a) = serteo;
            end % fskMs

        end % dofsk

    end % tsnrs

end % trellisArgs

%% Plots
figure(1); hold on;
title("BER w/ Conv. Ecoding vs SNR");

% Expermental data
for c = 1:length(trellisArgs)

    if dopsk

        for b = 1:length(pskMs)
            plot(tsnrs, max(eps, reshape(resultspskconv(1, c, b, :), [], 1)), '--*', ...
                'DisplayName', sprintf("PSK, M=%d; Conv. enc. %d", pskMs(b), c));
        end

    end

    if doqam

        for b = 1:length(qamMs)
            plot(tsnrs, max(eps, reshape(resultsqamconv(1, c, b, :), [], 1)), '--*', ...
                'DisplayName', sprintf("QAM, M=%d; Conv. enc. %d", qamMs(b), c));
        end

    end

    if dofsk

        for b = 1:length(fskMs)
            plot(tsnrs, max(eps, reshape(resultsfskconv(1, c, b, :), [], 1)), '--*', ...
                'DisplayName', sprintf("FSK, M=%d; Conv. enc. %d", fskMs(b), c));
        end

    end

end

% Theoretical lines
if dopsk

    for b = 1:length(pskMs)
        plot(tsnrs, max(eps, reshape(resultspskconv(3, 1, b, :), [], 1)), ...
            'DisplayName', sprintf("PSK, M=%d no FEC", pskMs(b)));
    end

end

if doqam

    for b = 1:length(qamMs)
        plot(tsnrs, max(eps, reshape(resultsqamconv(3, 1, b, :), [], 1)), ...
            'DisplayName', sprintf("QAM, M=%d no FEC", qamMs(b)));
    end

end

if dofsk

    for b = 1:length(fskMs)
        plot(tsnrs, max(eps, reshape(resultsfskconv(3, 1, b, :), [], 1)), ...
            'DisplayName', sprintf("FSK, M=%d no FEC", fskMs(b)));
    end

end

set(gca, 'YScale', 'log');
legend('Location', 'southwest', 'NumColumns', 1);
grid("minor");
hold off;
