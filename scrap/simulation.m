% DCS_Assignment_1920
% Rui Oliveira

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

% Cyclic codes settigs
docrc = true;
% CRC codings to try
crcCodings = [
        8 4;
        12 7; ];

% Convolutional code settings
doconv = true;
% Convolutional encodings to try
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

% Create the AWGN channel
awgnchan = comm.AWGNChannel("NoiseMethod", "Signal to noise ratio (SNR)");

% Exercise 1: Do no corrections
% Open an array for the results
% (BER, SER, BERteo, SERteo)
resultspsk_ex1 = zeros(4, length(pskMs), length(tsnrs));
resultsqam_ex1 = zeros(4, length(qamMs), length(tsnrs));
resultsfsk_ex1 = zeros(4, length(fskMs), length(tsnrs));

% PSK
if dopsk

    for b = 1:length(pskMs)

        % Get properties
        M = pskMs(b);
        nbits = log2(M);

        % First, fit the array into the symbols
        bitStreamFixed = [bitStream; randi([0 1], nbits - mod(length(bitStream), nbits), 1)];

        % Then, reshape the stream to convert into symbols
        databin = reshape(bitStreamFixed, [length(bitStreamFixed) / nbits, nbits]);
        % Get symbols from the binary stream
        data = bi2de(databin);

        % Now, create the modem object, and modulate the data
        if M == 2
            % For M = 2, initial phase is more convenient zero
            modulator = comm.PSKModulator("ModulationOrder", M, ...
                "PhaseOffset", 0, ...
                "SymbolMapping", "Gray");
            txSig = modulator(data);

            demodulator = comm.PSKDemodulator("ModulationOrder", M, ...
                "PhaseOffset", 0, ...
                "SymbolMapping", "Gray");
        else
            % Else, as per documentation
            modulator = comm.PSKModulator("ModulationOrder", M, ...
                "PhaseOffset", pi / M, ...
                "SymbolMapping", "Gray");
            txSig = modulator(data);

            demodulator = comm.PSKDemodulator("ModulationOrder", M, ...
                "PhaseOffset", pi / M, ...
                "SymbolMapping", "Gray");
        end

        % Then iterate
        for a = 1:length(tsnrs)

            tsnr = tsnrs(a);

            % Display simulation info
            fomatSpec = "No FEC; SNR = %d; PSK, M = %d\n";
            fprintf(fomatSpec, tsnr, M);

            % Set the AWGN channel SNR and apply
            awgnchan.SNR = tsnr;
            rxSig = awgnchan(txSig);

            % Demodulate data
            dataOut = demodulator(rxSig);

            % Convert data back into binary
            dataRxBin = de2bi(dataOut, log2(M));
            % Convert back into line vector
            dataRxBin = reshape(dataRxBin, [], 1);
            % Cut the received binary into the original shape
            dataRxBin = dataRxBin(1:length(bitStream));

            % Now check what we got
            [~, ber] = biterr(bitStream, dataRxBin);
            [~, ser] = symerr(data, dataOut);

            % And compute the theoretical value
            [berteo, serteo] = berawgn(tsnr - 10 * log10(nbits), "psk", M, 'nondiff');

            % And save
            resultspsk_ex1(1, b, a) = ber;
            resultspsk_ex1(2, b, a) = ser;
            resultspsk_ex1(3, b, a) = berteo;
            resultspsk_ex1(4, b, a) = serteo;

        end % tsnrs

    end % pskMs

end % dopsk

% QAM
if doqam

    for b = 1:length(qamMs)

        % Get properties
        M = qamMs(b);
        nbits = log2(M);

        % First, fit the array into the symbols
        bitStreamFixed = [bitStream; randi([0 1], nbits - mod(length(bitStream), nbits), 1)];

        % Then, reshape the stream to convert into symbols
        databin = reshape(bitStreamFixed, [length(bitStreamFixed) / nbits, nbits]);
        % Get symbols from the binary stream
        data = bi2de(databin);

        % Now, create the modem object, and modulate the data
        if M == 2
            % For M = 2, initial phase is more convenient zero
            modulator = comm.PSKModulator("ModulationOrder", M, ...
                "PhaseOffset", 0, ...
                "SymbolMapping", "Gray");
            txSig = modulator(data);

            demodulator = comm.PSKDemodulator("ModulationOrder", M, ...
                "PhaseOffset", 0, ...
                "SymbolMapping", "Gray");
        else
            % Else, as per documentation
            modulator = comm.PSKModulator("ModulationOrder", M, ...
                "PhaseOffset", pi / M, ...
                "SymbolMapping", "Gray");
            txSig = modulator(data);

            demodulator = comm.PSKDemodulator("ModulationOrder", M, ...
                "PhaseOffset", pi / M, ...
                "SymbolMapping", "Gray");
        end

        % Then iterate
        for a = 1:length(tsnrs)

            tsnr = tsnrs(a);

            % Display simulation info
            fomatSpec = "No FEC; SNR = %d; PSK, M = %d\n";
            fprintf(fomatSpec, tsnr, M);

            % Set the AWGN channel SNR and apply
            awgnchan.SNR = tsnr;
            rxSig = awgnchan(txSig);

            % Demodulate data
            dataOut = demodulator(rxSig);

            % Convert data back into binary
            dataRxBin = de2bi(dataOut, log2(M));
            % Convert back into line vector
            dataRxBin = reshape(dataRxBin, [], 1);
            % Cut the received binary into the original shape
            dataRxBin = dataRxBin(1:length(bitStream));

            % Now check what we got
            [~, ber] = biterr(bitStream, dataRxBin);
            [~, ser] = symerr(data, dataOut);

            % And compute the theoretical value
            [berteo, serteo] = berawgn(tsnr - 10 * log10(nbits), "qam", M);

            % And save
            resultsqam_ex1(1, b, a) = ber;
            resultsqam_ex1(2, b, a) = ser;
            resultsqam_ex1(3, b, a) = berteo;
            resultsqam_ex1(4, b, a) = serteo;

        end % tsnrs

    end % qamMs

end % doqam

%% Plots
% Exercise 1
%% Plots
figure(11); hold on;
title("BER vs SNR");

if dopsk

    for b = 1:length(pskMs)
        plot(tsnrs, max(eps, reshape(resultspsk_ex1(1, b, :), [], 1)), '--*', ...
            'DisplayName', sprintf("PSK, M=%d", pskMs(b)));
    end

end

ax = gca;
ax.ColorOrderIndex = 1;

if dopsk

    for b = 1:length(pskMs)
        plot(tsnrs, max(eps, reshape(resultspsk_ex1(3, b, :), [], 1)), ...
            'DisplayName', sprintf("PSK, M=%d", pskMs(b)));
    end

end

set(gca, 'YScale', 'log');
legend('Location', 'southwest', 'NumColumns', 1);
grid("minor");
hold off;

figure(12); hold on;
title("SER vs SNR");

if dopsk

    for b = 1:length(pskMs)
        plot(tsnrs, max(eps, reshape(resultspsk_ex1(2, b, :), [], 1)), '--*', ...
            'DisplayName', sprintf("PSK, M=%d", pskMs(b)));
    end

end

ax = gca;
ax.ColorOrderIndex = 1;

if dopsk

    for b = 1:length(pskMs)
        plot(tsnrs, max(eps, reshape(resultspsk_ex1(4, b, :), [], 1)), ...
            'DisplayName', sprintf("PSK, M=%d", pskMs(b)));
    end

end

set(gca, 'YScale', 'log');
legend('Location', 'southwest', 'NumColumns', 1);
grid("minor");
hold off;
