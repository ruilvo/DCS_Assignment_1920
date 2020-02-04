% DCS_Assignment_1920
% Rui Oliveira

close all;
clear all;

%--- Parameters ---
% Number of tx bits
N = 2^13;
% Sample frequency
fs = 20e6; % Hz

% OFDM parameters
nsubs = 64;
guardbands = [6; 6];
prefixlen = 16;
% Note on the assignment
% N = 64 - 6 - 6 = 52 useful subcariers, as intended
% 0.8us/(4us-0.8us) = 0.25
% 0.25 * 64 = 16, therefore 16 for the prefix length
% We don't have a notion of time here.
% It'll probably be important to find out how to select
% fs such as the symbol timings yield correct values.

% PSK settings
dopsk = true;
% Modulation orders
pskMs = [2, 4, 16];

% QAM settings
doqam = true;
% Modulation orders
qamMs = [4, 16];

% Generate cases to test
timedelays = linspace(0, 4e-6, 400);

testsDelays = {};

for a = 1:length(timedelays)
    testsDelays{a} = {[0 timedelays(a)], [0 -3]};
end

powers = 5:-1:-20; % dB

testsPowers = {};

for b = 1:length(powers)
    testsPowers{b} = {[0 0.4e-6], [0 powers(b)]};
end

npaths = 2:1:20;

testsNPaths = {};
delaystep = 0.1e-6;

for b = 1:length(npaths)
    np = npaths(b);
    delays = linspace(0, delaystep * np, np);
    powerscoise = [0, -9 * ones(1, np - 1)];
    testsNPaths{b} = {delays, powerscoise};
end

% Simulation
% Generate tx stream
bitStream = randi([0 1], N, 1);

% PSK
if dopsk
    results_psk_delays = zeros(length(pskMs), length(testsDelays));
    results_psk_powers = zeros(length(pskMs), length(testsPowers));
    results_psk_npaths = zeros(length(pskMs), length(testsNPaths));

    for d = 1:length(pskMs)

        for c = 1:length(testsDelays)

            fomatSpec = "PSK, M = %d; Delays, %d / %d\n";
            fprintf(fomatSpec, pskMs(d), c, length(testsDelays));
            drawnow;

            results_psk_delays(d, c) = modem_multipath_ofdm("psk", pskMs(d), bitStream, ...
                "fs", fs, "nsubs", nsubs, "guardbands", guardbands, "prefixlen", prefixlen, ...
                "pathdelays", testsDelays{c}{1}, "pathpwr", testsDelays{c}{2});
        end % results

        for c = 1:length(testsPowers)

            fomatSpec = "PSK, M = %d; Powers, %d / %d\n";
            fprintf(fomatSpec, pskMs(d), c, length(testsPowers));
            drawnow;

            results_psk_powers(d, c) = modem_multipath_ofdm("psk", pskMs(d), bitStream, ...
                "fs", fs, "nsubs", nsubs, "guardbands", guardbands, "prefixlen", prefixlen, ...
                "pathdelays", testsPowers{c}{1}, "pathpwr", testsPowers{c}{2});
        end % results

        for c = 1:length(testsNPaths)

            fomatSpec = "PSK, M = %d; N Paths, %d / %d\n";
            fprintf(fomatSpec, pskMs(d), c, length(testsNPaths));
            drawnow;

            results_psk_npaths(d, c) = modem_multipath_ofdm("psk", pskMs(d), bitStream, ...
                "fs", fs, "nsubs", nsubs, "guardbands", guardbands, "prefixlen", prefixlen, ...
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

        for c = 1:length(testsDelays)

            fomatSpec = "QAM, M = %d; Delays, %d / %d\n";
            fprintf(fomatSpec, qamMs(d), c, length(testsDelays));
            drawnow;

            results_qam_delays(d, c) = modem_multipath_ofdm("qam", qamMs(d), bitStream, ...
                "fs", fs, "nsubs", nsubs, "guardbands", guardbands, "prefixlen", prefixlen, ...
                "pathdelays", testsDelays{c}{1}, "pathpwr", testsDelays{c}{2});
        end % results

        for c = 1:length(testsPowers)

            fomatSpec = "QAM, M = %d; Powers, %d / %d\n";
            fprintf(fomatSpec, qamMs(d), c, length(testsPowers));
            drawnow;

            results_qam_powers(d, c) = modem_multipath_ofdm("qam", qamMs(d), bitStream, ...
                "fs", fs, "nsubs", nsubs, "guardbands", guardbands, "prefixlen", prefixlen, ...
                "pathdelays", testsPowers{c}{1}, "pathpwr", testsPowers{c}{2});
        end % results

        for c = 1:length(testsNPaths)

            fomatSpec = "QAM, M = %d; N Paths, %d / %d\n";
            fprintf(fomatSpec, qamMs(d), c, length(testsNPaths));
            drawnow;

            results_qam_npaths(d, c) = modem_multipath_ofdm("qam", qamMs(d), bitStream, ...
                "fs", fs, "nsubs", nsubs, "guardbands", guardbands, "prefixlen", prefixlen, ...
                "pathdelays", testsNPaths{c}{1}, "pathpwr", testsNPaths{c}{2});
        end % results

    end % qamMs

end % doqam

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
        plot(powers, max(1e-5, results_psk_powers(d, :)), "--*", ...
            'DisplayName', sprintf("PSK, M=%d", pskMs(d)));
    end

end % dopsk

if doqam

    for d = 1:length(qamMs)
        plot(powers, max(1e-5, results_qam_powers(d, :)), "--*", ...
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
