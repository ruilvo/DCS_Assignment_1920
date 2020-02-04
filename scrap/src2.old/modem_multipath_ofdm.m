function berout = modem_multipath_ofdm(modemtype, M, bitsin, varargin)
    % MODEM_AWGN Do the modem stuff
    % varargins fs, sps, pathdelays, pathpwr, nsubs, guardbands, prefixlen

    N = length(bitsin);
    nbits = log2(M);

    % Parse input
    p = inputParser;
    p.KeepUnmatched = true;
    addParameter(p, 'fs', 100);
    addParameter(p, 'pathdelays', [0, 2e-9]);
    addParameter(p, 'pathpwr', [0, -10]);
    addParameter(p, 'nsubs', 64);
    addParameter(p, 'guardbands', [6; 6]);
    addParameter(p, 'prefixlen', 16);

    parse(p, varargin{:});

    fs = p.Results.fs;
    pathdelays = p.Results.pathdelays;
    pathpwr = p.Results.pathpwr;
    nsubs = p.Results.nsubs;
    guardbands = p.Results.guardbands;
    prefixlen = p.Results.prefixlen;

    usedsubs = nsubs - sum(guardbands);
    frameSize = nbits * usedsubs;

    % First, fit the array into the symbols
    % The bitsin must fit within K frames
    symbstuffn = frameSize - mod(N, frameSize);
    bitStreamFixed = [bitsin; randi([0 1], symbstuffn, 1)];

    % Create the ODFM modem
    ofdmMod = comm.OFDMModulator("FFTLength", nsubs, "NumGuardBandCarriers", guardbands, ...
        "CyclicPrefixLength", prefixlen);
    ofdmDemod = comm.OFDMDemodulator("FFTLength", nsubs, "NumGuardBandCarriers", guardbands, ...
        "CyclicPrefixLength", prefixlen);
    ofdmDims = info(ofdmMod);

    % Modulate the data
    switch modemtype
        case "psk"
            % Create the PSK modem
            pskMod = comm.PSKModulator("ModulationOrder", M, "PhaseOffset", 0, ...
                "BitInput", true, "SymbolMapping", "Gray");
            pskDemod = comm.PSKDemodulator("ModulationOrder", M, "PhaseOffset", 0, ...
                "BitOutput", true, "SymbolMapping", "Gray", "DecisionMethod", "Hard decision");

            txSigUf = pskMod(bitStreamFixed);

        case "qam"
            txSigUf = qammod(bitStreamFixed, M, 'gray', "InputType", "bit");
    end % modemtype

    % Reshape symbols into frames
    nframes = length(txSigUf) / usedsubs;
    dataTxFramed = reshape(txSigUf, [usedsubs, nframes]);

    % Apply OFDM modulation;
    % We must do this frame by frame
    txSig = [];

    for k = 1:nframes
        tempSlice = ofdmMod(dataTxFramed(:, k));
        txSig = [txSig; tempSlice];
    end

    % Appy multipath
    rxSig = apply_multipath(txSig, pathdelays, pathpwr, fs);

    % Now to demodulate,
    % first, reshape rxSig to an approriate shape
    rxSigFramed = reshape(rxSig, ofdmDims.OutputSize(1), nframes);

    % Apply OFDM demod
    % We must do this frame by frame
    dataRx = [];

    for k = 1:nframes
        tempSlice = ofdmDemod(rxSigFramed(:, k));
        dataRx = [dataRx; tempSlice]; % Note the ;
    end

    % Now demodulate
    switch modemtype
        case "psk"
            dataOut = pskDemod(dataRx);
        case "qam"
            dataOut = qamdemod(dataRx, M, 'gray', "OutputType", "bit");
    end % modemtype

    [~, berout] = biterr(bitStreamFixed, dataOut);

end % function
