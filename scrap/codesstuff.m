% CRC codings to try
cyclicCodes = {[10 2], ...
    [8 4]};

for i = 1:length(cyclicCodes)
    clen = cyclicCodes{i}(1);
    mlen = cyclicCodes{i}(2);
    
    % Now the message can fit into the code
    % Create a generator polynomial for a cyclic code.
    gpol = cyclpoly(clen, mlen);
    % Create a parity-check matrix by using the generator polynomial.
    parmat = cyclgen(clen, gpol);
    % Create a syndrome decoding table by using the parity-check matrix.
    trt = syndtable(parmat);
    
    wt = gfweight(parmat,'par')
    
    % This is the maximum correcteable errors
    t = (wt-1)/2
end