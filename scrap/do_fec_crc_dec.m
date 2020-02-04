function dataRecBinary = do_fec_crc_dec(encDataRxInt,nbits,n,k,gpol,trt)
%Do the recovery of the previous FEC

% First convert data back into binary
dataRxBin = de2bi(encDataRxInt,nbits);
% Then get back the proper size again
% Now I do it by removing data just becausue
dataRxBinFixed = reshape(dataRxBin,1,[]);
while mod(length(dataRxBinFixed),n)
    dataRxBinFixed(end) = []; % Fill in a bit
end
% Decode the corrupted sequence
dataRecBinary = decode(dataRxBinFixed,n,k,'cyclic/binary',gpol,trt);
end

