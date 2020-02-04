function rxSig = apply_multipath(txSig, pathdelays, pathpwr, fs)
    % Apply multipath dispersion
    %      !------------------Fixed gain-------!
    %      !                                   !
    %s(t)--!signal in                         Sum--------r(t)signal out
    %      !                                   !
    %      !--(xN)---Delay----Variable gain----!
    %                         and attenuation

    if length(pathdelays) ~= length(pathpwr)
        error("Different lengths");
    end

    pathpwrlin = db2pow(pathpwr);
    rxSig = zeros(length(txSig), 1);

    for a = 1:length(pathdelays)
        rxSig = rxSig + pathpwrlin(a) * fft_delay(txSig, pathdelays(a), fs);
    end % pathdelays

end
