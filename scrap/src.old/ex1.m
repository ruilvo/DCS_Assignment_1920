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

%--- Simulation ---
% Generate tx stream
bitStream = randi([0 1], N, 1);

% Open an array for the results
resultspsk = zeros(4, length(pskMs), length(tsnrs));
resultsqam = zeros(4, length(qamMs), length(tsnrs));
resultsfsk = zeros(4, length(fskMs), length(tsnrs));

for a = 1:length(tsnrs)
    tsnr = tsnrs(a);

    % PSK
    if dopsk

        for b = 1:length(pskMs)
            M = pskMs(b);
            nbits = log2(M);

            fomatSpec = "SNR = %d; PSK, M = %d\n";
            fprintf(fomatSpec, tsnr, M);

            % First, fit the array into the symbols
            % We need to add nbits-mod(len,nbits) bits
            bitStreamFixed = [bitStream; randi([0 1], nbits - mod(length(bitStream), nbits), 1)];

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

            % Now we only want the bits of the original stream
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
            [berteo, serteo] = berawgn(tsnr - 10 * log10(log2(M)), "psk", M, 'nondiff');

            resultspsk(1, b, a) = ber;
            resultspsk(2, b, a) = ser;
            resultspsk(3, b, a) = berteo;
            resultspsk(4, b, a) = serteo;
        end % pskMs

    end % dopsk

    % QAM
    if doqam

        for b = 1:length(qamMs)
            M = qamMs(b);
            nbits = log2(M);

            fomatSpec = "SNR = %d; QAM, M = %d\n";
            fprintf(fomatSpec, tsnr, M);

            % First, fit the array into the symbols
            % We need to add nbits-mod(len,nbits) symbols
            bitStreamFixed = [bitStream; randi([0 1], nbits - mod(length(bitStream), nbits), 1)];

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

            % Now we only want the bits of the original stream
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
            [berteo, serteo] = berawgn(tsnr - 10 * log10(log2(M)), "qam", M);

            resultsqam(1, b, a) = ber;
            resultsqam(2, b, a) = ser;
            resultsqam(3, b, a) = berteo;
            resultsqam(4, b, a) = serteo;
        end % qamMs

    end % doqam

    % FSK
    if dofsk

        for b = 1:length(fskMs)
            M = fskMs(b);
            nbits = log2(M);

            fomatSpec = "SNR = %d; FSK, M = %d\n";
            fprintf(fomatSpec, tsnr, M);

            % First, fit the array into the symbols
            % We need to add nbits-mod(len,nbits) symbols
            bitStreamFixed = [bitStream; randi([0 1], nbits - mod(length(bitStream), nbits), 1)];

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

            % Now we only want the bits of the original stream
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
            % https://www.mathworks.com/help/comm/ug/awgn-channel.html#a1071501088
            [berteo, serteo] = berawgn(tsnr - 10 * log10(log2(M)) + 10 * log10(nsamp), "fsk", M, 'noncoherent');

            resultsfsk(1, b, a) = ber;
            resultsfsk(2, b, a) = ser;
            resultsfsk(3, b, a) = berteo;
            resultsfsk(4, b, a) = serteo;
        end % fskMs

    end % dofsk

end % tsnrs

%% Plots
figure(1); hold on;
title("BER vs SNR");

if dopsk

    for b = 1:length(pskMs)
        plot(tsnrs, max(eps, reshape(resultspsk(1, b, :), [], 1)), '--*', ...
            'DisplayName', sprintf("PSK, M=%d", pskMs(b)));
    end

end

if doqam

    for b = 1:length(qamMs)
        plot(tsnrs, max(eps, reshape(resultsqam(1, b, :), [], 1)), '--*', ...
            'DisplayName', sprintf("QAM, M=%d", qamMs(b)));
    end

end

if dofsk

    for b = 1:length(fskMs)
        plot(tsnrs, max(eps, reshape(resultsfsk(1, b, :), [], 1)), '--*', ...
            'DisplayName', sprintf("FSK, M=%d", fskMs(b)));
    end

end

ax = gca;
ax.ColorOrderIndex = 1;

if dopsk

    for b = 1:length(pskMs)
        plot(tsnrs, max(eps, reshape(resultspsk(3, b, :), [], 1)), ...
            'DisplayName', sprintf("PSK, M=%d", pskMs(b)));
    end

end

if doqam

    for b = 1:length(qamMs)
        plot(tsnrs, max(eps, reshape(resultsqam(3, b, :), [], 1)), ...
            'DisplayName', sprintf("QAM, M=%d", qamMs(b)));
    end

end

if dofsk

    for b = 1:length(fskMs)
        plot(tsnrs, max(eps, reshape(resultsfsk(3, b, :), [], 1)), ...
            'DisplayName', sprintf("FSK, M=%d", fskMs(b)));
    end

end

set(gca, 'YScale', 'log');
legend('Location', 'southwest', 'NumColumns', 1);
grid("minor");
hold off;

figure(2); hold on;
title("SER vs SNR");

if dopsk

    for b = 1:length(pskMs)
        plot(tsnrs, max(eps, reshape(resultspsk(2, b, :), [], 1)), '--*', ...
            'DisplayName', sprintf("PSK, M=%d", pskMs(b)));
    end

end

if doqam

    for b = 1:length(qamMs)
        plot(tsnrs, max(eps, reshape(resultsqam(2, b, :), [], 1)), '--*', ...
            'DisplayName', sprintf("QAM, M=%d", qamMs(b)));
    end

end

if dofsk

    for b = 1:length(fskMs)
        plot(tsnrs, max(eps, reshape(resultsfsk(2, b, :), [], 1)), '--*', ...
            'DisplayName', sprintf("FSK, M=%d", fskMs(b)));
    end

end

ax = gca;
ax.ColorOrderIndex = 1;

if dopsk

    for b = 1:length(pskMs)
        plot(tsnrs, max(eps, reshape(resultspsk(4, b, :), [], 1)), ...
            'DisplayName', sprintf("PSK, M=%d", pskMs(b)));
    end

end

if doqam

    for b = 1:length(qamMs)
        plot(tsnrs, max(eps, reshape(resultsqam(4, b, :), [], 1)), ...
            'DisplayName', sprintf("QAM, M=%d", qamMs(b)));
    end

end

if dofsk

    for b = 1:length(fskMs)
        plot(tsnrs, max(eps, reshape(resultsfsk(4, b, :), [], 1)), ...
            'DisplayName', sprintf("FSK, M=%d", fskMs(b)));
    end

end

set(gca, 'YScale', 'log');
legend('Location', 'southwest', 'NumColumns', 1);
grid("minor");
hold off;
