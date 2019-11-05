rng(130912);
userConfig = struct('xy', gr137mat, 'dmat', gr137mat, 'nSalesmen', 1, 'numIter', 5000, 'showProg',false,'showResult',false);
resultStruct = mtspf_ga(userConfig);
disp(resultStruct.minDist);
