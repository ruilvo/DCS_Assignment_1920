% DCS_Assignment_1920
% Rui Oliveira

close all;
clear all;

%--- Parameters ---
% Number of tx bits
N = 2^10;
% Baud rate
br = 1/(0.4e-6); % Samples per second
% Sample frequency
fs = 20e6; % Hz
% Samples per symbol
sps = fs / br;

if ~(sps == fix(sps))
    error("Fix sampling/baud rate");
end % baud rate check

% RC filter parameters
rcbeta = 0.2;
rcspan = 10;

% PSK settings
dopsk = true;
% Modulation orders
pskMs = [4, 8, 16, 32, 64];

% QAM settings
doqam = true;
% Modulation orders
qamMs = [4, 16, 64];

% Generate cases to test

% Excess time delay
timedelays = linspace(0, 4e-6, 400);
scndchanpwr = -3; % dB

testsDelays = {};

for a = 1:length(timedelays)
    testsDelays{a} = {[0 timedelays(a)], [0 scndchanpwr]};
end

% Power of the reflected channel
relpowers = 5:-1:-20; % dB
scndchandelay = 0.5e-6; % s

testsPowers = {};

for b = 1:length(relpowers)
    testsPowers{b} = {[0 scndchandelay], [0 relpowers(b)]};
end

% Number of concurrent paths
npaths = 2:1:20;
delaystep = 0.1e-6;
otherchanpwrs = -9;

testsNPaths = {};

testsNPaths{1} = {[0], [0]};

for b = 1:length(npaths)
    np = npaths(b);
    delays = linspace(0, delaystep * np, np);
    powersvect = [0, otherchanpwrs * ones(1, np - 1)];
    testsNPaths{b + 1} = {delays, powersvect};
end

npaths = [1 npaths];

%% Simulation --------------
% Generate tx stream
bitStream = randi([0 1], N, 1);

% PSK
if dopsk
    results_psk_delays = zeros(length(pskMs), length(testsDelays));
    results_psk_powers = zeros(length(pskMs), length(testsPowers));
    results_psk_npaths = zeros(length(pskMs), length(testsNPaths));

    for d = 1:length(pskMs)

        M = pskMs(d);

        % Excess time delay
        switch M
            case 2
                scndchanpwr = -3;
            case 4
                scndchanpwr = -3;
            case 8
                scndchanpwr = -6;
            case 16
                scndchanpwr = -9;
            case 32
                scndchanpwr = -12;
            case 64
                scndchanpwr = -15;
            otherwise
                scndchanpwr = -9;
        end % M

        for a = 1:length(timedelays)
            testsDelays{a}{2} = [0 scndchanpwr];
        end

        for c = 1:length(testsDelays)

            fomatSpec = "PSK, M = %d; Delays, %d / %d\n";
            fprintf(fomatSpec, pskMs(d), c, length(testsDelays));
            drawnow;

            results_psk_delays(d, c) = modem_multipath("psk", pskMs(d), bitStream, ...
                "fs", fs, "sps", sps, "rcbeta", rcbeta, "rcspan", rcspan, ...
                "pathdelays", testsDelays{c}{1}, "pathpwr", testsDelays{c}{2});
        end % results

        for c = 1:length(testsPowers)

            fomatSpec = "PSK, M = %d; Powers, %d / %d\n";
            fprintf(fomatSpec, pskMs(d), c, length(testsPowers));
            drawnow;

            results_psk_powers(d, c) = modem_multipath("psk", pskMs(d), bitStream, ...
                "fs", fs, "sps", sps, "rcbeta", rcbeta, "rcspan", rcspan, ...
                "pathdelays", testsPowers{c}{1}, "pathpwr", testsPowers{c}{2});
        end % results

        for c = 1:length(testsNPaths)

            fomatSpec = "PSK, M = %d; N Paths, %d / %d\n";
            fprintf(fomatSpec, pskMs(d), c, length(testsNPaths));
            drawnow;

            results_psk_npaths(d, c) = modem_multipath("psk", pskMs(d), bitStream, ...
                "fs", fs, "sps", sps, "rcbeta", rcbeta, "rcspan", rcspan, ...
                "pathdelays", testsNPaths{c}{1}, "pathpwr", testsNPaths{c}{2});
        end % results

    end % pskMs

end % dopsk

% QAM
if doqam
    results_qam_delays = zeros(length(qamMs), length(testsDelays));
    results_qam_powers = zeros(length(qamMs), length(testsPowers));
    results_qam_npaths = zeros(length(qamMs), length(testsNPaths));

    for d = 1:length(qamMs)

        M = qamMs(d);

        % Excess time delay
        switch M
            case 2
                scndchanpwr = -3;
            case 4
                scndchanpwr = -3;
            case 16
                scndchanpwr = -8;
            case 64
                scndchanpwr = -10;
            otherwise
                scndchanpwr = -9;
        end % M

        for a = 1:length(timedelays)
            testsDelays{a}{2} = [0 scndchanpwr];
        end

        for c = 1:length(testsDelays)

            fomatSpec = "QAM, M = %d; Delays, %d / %d\n";
            fprintf(fomatSpec, qamMs(d), c, length(testsDelays));
            drawnow;

            results_qam_delays(d, c) = modem_multipath("qam", qamMs(d), bitStream, ...
                "fs", fs, "sps", sps, "rcbeta", rcbeta, "rcspan", rcspan, ...
                "pathdelays", testsDelays{c}{1}, "pathpwr", testsDelays{c}{2});
        end % results

        for c = 1:length(testsPowers)

            fomatSpec = "QAM, M = %d; Powers, %d / %d\n";
            fprintf(fomatSpec, qamMs(d), c, length(testsPowers));
            drawnow;

            results_qam_powers(d, c) = modem_multipath("qam", qamMs(d), bitStream, ...
                "fs", fs, "sps", sps, "rcbeta", rcbeta, "rcspan", rcspan, ...
                "pathdelays", testsPowers{c}{1}, "pathpwr", testsPowers{c}{2});
        end % results

        for c = 1:length(testsNPaths)

            fomatSpec = "QAM, M = %d; N Paths, %d / %d\n";
            fprintf(fomatSpec, qamMs(d), c, length(testsNPaths));
            drawnow;

            results_qam_npaths(d, c) = modem_multipath("qam", qamMs(d), bitStream, ...
                "fs", fs, "sps", sps, "rcbeta", rcbeta, "rcspan", rcspan, ...
                "pathdelays", testsNPaths{c}{1}, "pathpwr", testsNPaths{c}{2});
        end % results

    end % qamMs

end % doqam

% Save data
save("../results/resultsmultipath.mat", "timedelays", "results*", "pskMs", "qamMs", ...
    "relpowers", "npaths");

%% Plots -------------------------------

figure(1); hold on;

title("BER vs relative chan. delay (2 chan)");

if dopsk

    for d = 1:length(pskMs)
        plot(timedelays, max(1e-5, results_psk_delays(d, :)), ...
            'DisplayName', sprintf("PSK, M=%d", pskMs(d)));
    end

end % dopsk

if doqam

    for d = 1:length(qamMs)
        plot(timedelays, max(1e-5, results_qam_delays(d, :)), ...
            'DisplayName', sprintf("QAM, M=%d", qamMs(d)));
    end

end % doqam

set(gca, 'YScale', 'log');
legend('Location', 'best', 'NumColumns', 1);
grid("minor");
hold off;

figure(2); hold on;

title("BER vs relative chan. powers (2 chan)");

if dopsk

    for d = 1:length(pskMs)
        plot(relpowers, max(1e-5, results_psk_powers(d, :)), "--*", ...
            'DisplayName', sprintf("PSK, M=%d", pskMs(d)));
    end

end % dopsk

if doqam

    for d = 1:length(qamMs)
        plot(relpowers, max(1e-5, results_qam_powers(d, :)), "--*", ...
            'DisplayName', sprintf("QAM, M=%d", qamMs(d)));
    end

end % doqam

set(gca, 'YScale', 'log');
legend('Location', 'best', 'NumColumns', 1);
grid("minor");
hold off;

figure(3); hold on;

title("BER vs number of chans.");

if dopsk

    for d = 1:length(pskMs)
        plot(npaths, max(1e-5, results_psk_npaths(d, :)), "--*", ...
            'DisplayName', sprintf("PSK, M=%d", pskMs(d)));
    end

end % dopsk

if doqam

    for d = 1:length(qamMs)
        plot(npaths, max(1e-5, results_qam_npaths(d, :)), "--*", ...
            'DisplayName', sprintf("QAM, M=%d", qamMs(d)));
    end

end % doqam

set(gca, 'YScale', 'log');
legend('Location', 'best', 'NumColumns', 1);
grid("minor");
hold off;
