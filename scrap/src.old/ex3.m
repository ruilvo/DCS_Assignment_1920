% DCS_Assignment_1920
% Rui Oliveira
% Simulations for the Multipath channel

close all;
clear all;

%--- Parameters ---
% Number of tx bits
N = 2^18;

% Modulation orders for PSK
dopsk = true;
pskMs = [2];
% Modulation orders for QAM
doqam = true;
qamMs = [4, 64];
% Parameters for FSK
dofsk = true;
fskMs = [2];
freqsep = 8; % Frequency separation (Hz)
nsamp = 16; % Number of samples per symbol
Fs = freqsep * max(fskMs) * 4; % Sample rate (Hz)

% Settings for the multipath channel
fs = 3.84e6; % Hz

for k = 1:1000
    pathDelays{k} = [0, 0.000005*k] * 1e-9;
    avgPathGains{k} = [0, -10];
end

% pathDelays = {...
%     [0, 1] * 1e-9, ...
%     [0, 0.05] * 1e-9, ...
%     [0, 0.005] * 1e-9, ...
%     [0, 0.005, 0.007, 0.010] * 1e-9, ...
%     }; % sec
% avgPathGains = {...
%     [0, -10], ...
%     [0, -3], ...
%     [0, -3], ...
%     [0, -6, -9, -12], ...
%     }; % dB
fD = 50; % Hz, maximum Doppler shift

% Generate tx stream
bitStream = randi([0 1], N, 1);

%--- Simulation ---
% Basic sanity check
if length(pathDelays) ~= length(avgPathGains)
    error("Wrong sizes");
end

% Open an array for the results
resultspsk = zeros(4, length(pskMs), length(avgPathGains));
resultsqam = zeros(4, length(qamMs), length(avgPathGains));
resultsfsk = zeros(4, length(fskMs), length(avgPathGains));

for a = 1:length(avgPathGains)

    % Create the channel
    rchan = comm.RayleighChannel(...
        'SampleRate', fs, ...
        'PathDelays', pathDelays{a}, ...
        'AveragePathGains', avgPathGains{a}, ...
        'MaximumDopplerShift', fD ...
        );

    % PSK
    if dopsk

        for b = 1:length(pskMs)

            M = pskMs(b);
            nbits = log2(M);

            fomatSpec = "Case = %d; PSK, M = %d\n";
            fprintf(fomatSpec, a, M);

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

            % Apply channel dispersion
            rxSig = rchan(txSig);

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

            resultspsk(1, b, a) = ber;
            resultspsk(2, b, a) = ser;

        end % pskMs

    end % dopsk

end % tdels

%% Plots
figure(1); hold on;
title("BER vs stuff");

if dopsk

    for b = 1:length(pskMs)
        plot(max(eps, reshape(resultspsk(1, b, :), [], 1)), '--*', ...
            'DisplayName', sprintf("PSK, M=%d", pskMs(b)));
    end

end
