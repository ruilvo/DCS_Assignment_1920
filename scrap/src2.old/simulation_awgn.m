% DCS_Assignment_1920
% Rui Oliveira

close all;
clear all;

%--- Parameters ---
% Number of tx bits
N = 2^13;
% Baud rate
br = 2e3; % Samples per second
% Sample frequency
fs = 1e4; % Hz
% Samples per symbol
sps = fs / br;

if ~(sps == fix(sps))
    error("Fix sampling/baud rate");
end % baud rate check

% RC filter parameters
rcbeta = 0.2;
rcspan = 10;

% Simulations to run
% Run basic simulation
donofec = true;

% Cyclic FEC settings
docyclic = true;
% CRC codings to try
crcCodings = {[10 2], ...
    [8 4]};

% Convolutional FEC
doconvolutional = true;
% Conv encodings to try
% Must generate compatible arguments
% Two arguments = no feedback
% Three argument = feedback
trellisArgsToTry = {
{[5 4], [23 35 0; 0 5 13]}, % 2/3 Feedforward Convolutional Encoder
{7, {'1 + x^3 + x^4 + x^5 + x^6', '1 + x + x^3 + x^4 + x^6'}}, % 1/2 Feedforward Convolutional Encoder
{5, [37 33], 37}}; % 1/2 Feedback Convolutional Encoders

% SNRs to try (dB)
minsnr = -30;
maxsnr = 20;
snrstep = 2;
tsnrs = minsnr:snrstep:maxsnr;

% PSK settings
dopsk = true;
% Modulation orders
pskMs = [2, 4, 64];

% QAM settings
doqam = true;
% Modulation orders
qamMs = [4, 64];

% FSK settings
dofsk = true;
fsep = 50; % Frequency separation (Hz)
% Modulation orders (number of frequencies)
fskMs = [2 32];

%--- Simulation ---
% Generate tx stream
bitStream = randi([0 1], N, 1);

% Results containers
resultspsk = zeros(2, length(pskMs), length(tsnrs));
resultsqam = zeros(2, length(qamMs), length(tsnrs));
resultsfsk = zeros(2, length(fskMs), length(tsnrs));

resultspsk_crc = zeros(length(crcCodings), length(pskMs), length(tsnrs));
resultsqam_crc = zeros(length(crcCodings), length(qamMs), length(tsnrs));
resultsfsk_crc = zeros(length(crcCodings), length(fskMs), length(tsnrs));

resultspsk_conv = zeros(length(trellisArgsToTry), length(pskMs), length(tsnrs));
resultsqam_conv = zeros(length(trellisArgsToTry), length(qamMs), length(tsnrs));
resultsfsk_conv = zeros(length(trellisArgsToTry), length(fskMs), length(tsnrs));

% Always get theoretical values
for b = 1:length(pskMs)
    M = pskMs(b);

    berteo_psk = zeros(length(tsnrs), 1);

    for a = 1:length(tsnrs)
        tsnr = tsnrs(a);
        [berteo_pski, serteo_pski] = berawgn(tsnr - 10 * log10(log2(M)) + 10 * log10(sps), "psk", M, 'nondiff');
        berteo_psk(a) = berteo_pski;
    end % tsnrs

    resultspsk(2, b, :) = berteo_psk;

end % pskMs

for b = 1:length(qamMs)
    M = qamMs(b);

    berteo_qam = zeros(length(tsnrs), 1);

    for a = 1:length(tsnrs)
        tsnr = tsnrs(a);
        [berteo_qami, serteo_qami] = berawgn(tsnr - 10 * log10(log2(M)) + 10 * log10(sps), "qam", M);
        berteo_qam(a) = berteo_qami;
    end % tsnrs

    resultsqam(2, b, :) = berteo_qam;

end % qamMs

for b = 1:length(fskMs)
    M = fskMs(b);

    berteo_fsk = zeros(length(tsnrs), 1);

    for a = 1:length(tsnrs)
        tsnr = tsnrs(a);
        [berteo_fski, serteo_fski] = berawgn(tsnr - 10 * log10(log2(M)) + 10 * log10(sps), "fsk", M, 'noncoherent');
        berteo_fsk(a) = berteo_fski;
    end % tsnrs

    resultsfsk(2, b, :) = berteo_fsk;

end % qamMs

% No FEC simulations
if donofec

    fprintf("No FEC simulations\n");

    % PSK
    if dopsk

        parfor b = 1:length(pskMs)
            M = pskMs(b);

            % Do the experiment
            berout_psk = modem_awgn('psk', pskMs(b), bitStream, tsnrs, ...
                'sps', sps, 'rcbeta', rcbeta, 'rcspan', rcspan);

            resultspsk(1, b, :) = berout_psk;

        end % pskMs

    end % dopsk

    % QAM
    if doqam

        parfor b = 1:length(qamMs)
            M = qamMs(b);

            % Do the experiment
            berout_qam = modem_awgn('qam', qamMs(b), bitStream, tsnrs, ...
                'sps', sps, 'rcbeta', rcbeta, 'rcspan', rcspan);

            resultsqam(1, b, :) = berout_qam;

        end % qamMs

    end % doqam

    % FSK
    if dofsk

        parfor b = 1:length(fskMs)
            M = fskMs(b);

            % Do the experiment
            berout_fsk = modem_awgn('fsk', fskMs(b), bitStream, tsnrs, ...
                'sps', sps, 'fs', fs, 'fsep', fsep);

            resultsfsk(1, b, :) = berout_fsk;

        end % qamMs

    end % dofsk

end % donofec

% Cyclic FEC simulations
if docyclic

    for c = 1:length(crcCodings)

        fprintf("Cyclic codes %d / %d\n", c, length(crcCodings));

        crcCoding = crcCodings{c};

        % PSK
        if dopsk

            for b = 1:length(pskMs)
                M = pskMs(b);

                % Do the experiment
                berout_psk = modem_awgn('psk', pskMs(b), bitStream, tsnrs, ...
                    'sps', sps, 'rcbeta', rcbeta, 'rcspan', rcspan, ...
                    "FECtype", "cyclic", "crcCoding", crcCoding);

                resultspsk_crc(c, b, :) = berout_psk;

            end % pskMs

        end % dopsk

        % QAM
        if doqam

            for b = 1:length(qamMs)
                M = qamMs(b);

                % Do the experiment
                berout_qam = modem_awgn('qam', qamMs(b), bitStream, tsnrs, ...
                    'sps', sps, 'rcbeta', rcbeta, 'rcspan', rcspan, ...
                    "FECtype", "cyclic", "crcCoding", crcCoding);

                resultsqam_crc(c, b, :) = berout_qam;

            end % qamMs

        end % doqam

        % FSK
        if dofsk

            for b = 1:length(fskMs)
                M = fskMs(b);

                % Do the experiment
                berout_fsk = modem_awgn('fsk', fskMs(b), bitStream, tsnrs, ...
                    'sps', sps, 'fs', fs, 'fsep', fsep, ...
                    "FECtype", "cyclic", "crcCoding", crcCoding);

                resultsfsk_crc(c, b, :) = berout_fsk;

            end % qamMs

        end % dofsk

    end % crcCodings

end % docyclic

% Convolutional FEC simulations
if doconvolutional

    for c = 1:length(trellisArgsToTry)

        fprintf("Convolutional codes %d / %d\n", c, length(trellisArgsToTry));

        trellisArgs = trellisArgsToTry{c};

        % PSK
        if dopsk

            for b = 1:length(pskMs)
                M = pskMs(b);

                % Do the experiment
                berout_psk = modem_awgn('psk', pskMs(b), bitStream, tsnrs, ...
                    'sps', sps, 'rcbeta', rcbeta, 'rcspan', rcspan, ...
                    "FECtype", "convolutional", "trellisArgs", trellisArgs);

                resultspsk_conv(c, b, :) = berout_psk;

            end % pskMs

        end % dopsk

        % QAM
        if doqam

            for b = 1:length(qamMs)
                M = qamMs(b);

                % Do the experiment
                berout_qam = modem_awgn('qam', qamMs(b), bitStream, tsnrs, ...
                    'sps', sps, 'rcbeta', rcbeta, 'rcspan', rcspan, ...
                    "FECtype", "convolutional", "trellisArgs", trellisArgs);

                resultsqam_conv(c, b, :) = berout_qam;

            end % qamMs

        end % doqam

        % FSK
        if dofsk

            for b = 1:length(fskMs)
                M = fskMs(b);

                % Do the experiment
                berout_fsk = modem_awgn('fsk', fskMs(b), bitStream, tsnrs, ...
                    'sps', sps, 'fs', fs, 'fsep', fsep, ...
                    "FECtype", "convolutional", "trellisArgs", trellisArgs);

                resultsfsk_conv(c, b, :) = berout_fsk;

            end % qamMs

        end % dofsk

    end % trellisArgs

end % doconv

%% Plots ----------------------------------------------------

% No FEC
if donofec
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

    % Theoretical lines
    ax = gca;
    ax.ColorOrderIndex = 1;

    if dopsk

        for b = 1:length(pskMs)
            plot(tsnrs, max(eps, reshape(resultspsk(2, b, :), [], 1)), ...
                'DisplayName', sprintf("PSK, M=%d", pskMs(b)));
        end

    end

    if doqam

        for b = 1:length(qamMs)
            plot(tsnrs, max(eps, reshape(resultsqam(2, b, :), [], 1)), ...
                'DisplayName', sprintf("QAM, M=%d", qamMs(b)));
        end

    end

    if dofsk

        for b = 1:length(fskMs)
            plot(tsnrs, max(eps, reshape(resultsfsk(2, b, :), [], 1)), ...
                'DisplayName', sprintf("FSK, M=%d", fskMs(b)));
        end

    end

    set(gca, 'YScale', 'log');
    legend('Location', 'southwest', 'NumColumns', 1);
    grid("minor");
    hold off;
end % donofec

% Cyclic codes
if docyclic
    figure(2); hold on;

    for c = 1:length(crcCodings)

        if dopsk

            for b = 1:length(pskMs)
                plot(tsnrs, max(eps, reshape(resultspsk_crc(c, b, :), [], 1)), '--*', ...
                    'DisplayName', sprintf("PSK, M=%d; CRC %d", pskMs(b), c));

            end % pskMs

        end % dopsk

        if doqam

            for b = 1:length(qamMs)
                plot(tsnrs, max(eps, reshape(resultsqam_crc(c, b, :), [], 1)), '--*', ...
                    'DisplayName', sprintf("QAM, M=%d; CRC %d", qamMs(b), c));

            end % qamMs

        end % doqam

        if dofsk

            for b = 1:length(fskMs)
                plot(tsnrs, max(eps, reshape(resultsfsk_crc(c, b, :), [], 1)), '--*', ...
                    'DisplayName', sprintf("FSK, M=%d; CRC %d", fskMs(b), c));

            end % qamMs

        end % doqam

    end % crcCodings

    % Theoretical lines
    if dopsk

        for b = 1:length(pskMs)
            plot(tsnrs, max(eps, reshape(resultspsk(2, b, :), [], 1)), ...
                'DisplayName', sprintf("PSK, M=%d", pskMs(b)));
        end

    end

    if doqam

        for b = 1:length(qamMs)
            plot(tsnrs, max(eps, reshape(resultsqam(2, b, :), [], 1)), ...
                'DisplayName', sprintf("QAM, M=%d", qamMs(b)));
        end

    end

    if dofsk

        for b = 1:length(fskMs)
            plot(tsnrs, max(eps, reshape(resultsfsk(2, b, :), [], 1)), ...
                'DisplayName', sprintf("FSK, M=%d", fskMs(b)));
        end

    end

    set(gca, 'YScale', 'log');
    legend('Location', 'southwest', 'NumColumns', 1);
    grid("minor");
    hold off;

end % docyclic

% Convolutional codes
if doconvolutional
    figure(3); hold on;

    for c = 1:length(trellisArgsToTry)

        if dopsk

            for b = 1:length(pskMs)
                plot(tsnrs, max(eps, reshape(resultspsk_conv(c, b, :), [], 1)), '--*', ...
                    'DisplayName', sprintf("PSK, M=%d; Conv %d", pskMs(b), c));

            end % pskMs

        end % dopsk

        if doqam

            for b = 1:length(qamMs)
                plot(tsnrs, max(eps, reshape(resultsqam_conv(c, b, :), [], 1)), '--*', ...
                    'DisplayName', sprintf("QAM, M=%d; Conv %d", qamMs(b), c));

            end % qamMs

        end % doqam

        if dofsk

            for b = 1:length(fskMs)
                plot(tsnrs, max(eps, reshape(resultsfsk_conv(c, b, :), [], 1)), '--*', ...
                    'DisplayName', sprintf("FSK, M=%d; Conv %d", fskMs(b), c));

            end % qamMs

        end % doqam

    end % trellisArgsToTry

    % Theoretical lines
    if dopsk

        for b = 1:length(pskMs)
            plot(tsnrs, max(eps, reshape(resultspsk(2, b, :), [], 1)), ...
                'DisplayName', sprintf("PSK, M=%d", pskMs(b)));
        end

    end

    if doqam

        for b = 1:length(qamMs)
            plot(tsnrs, max(eps, reshape(resultsqam(2, b, :), [], 1)), ...
                'DisplayName', sprintf("QAM, M=%d", qamMs(b)));
        end

    end

    if dofsk

        for b = 1:length(fskMs)
            plot(tsnrs, max(eps, reshape(resultsfsk(2, b, :), [], 1)), ...
                'DisplayName', sprintf("FSK, M=%d", fskMs(b)));
        end

    end

    set(gca, 'YScale', 'log');
    legend('Location', 'southwest', 'NumColumns', 1);
    grid("minor");
    hold off;

end % doconvolutional
