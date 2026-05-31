function [learnR, confR] = learning_per_neuron(signal, np_fs, trigDelay, cW, normed)
    
    binT = 100;
    Ntemp = size(signal, 3);
    cT = round([trigDelay-cW trigDelay]*np_fs);
    frateW = 0.01;
    learnR = zeros(size(signal, 3), 1);
    confR = zeros(size(signal, 3), 2);
    marginW = 0.025;
    for i = 1:Ntemp


        meanR = sum(signal(cT(1):cT(2), :, i), 1)'/cW;

        if normed
            frate = smoothdata(signal(:,:,i), 'gaussian', frateW*np_fs)*np_fs;
            binSig = smoothdata(frate, 2, 'gaussian', binT);
            relR = meanR./median(binSig([1:round((trigDelay-cW-marginW)*np_fs) round((trigDelay+marginW)*np_fs):end], :), 1)';
            relR = relR*median(binSig([1:round((trigDelay-cW-marginW)*np_fs) round((trigDelay+marginW)*np_fs):end], 1), 1);
        else
            relR = meanR;
        end
            
        if ~anynan(relR)
            fitobject = fit((1:size(relR, 1))', relR, 'poly1');
            learnR(i) = fitobject.p1;
            c = confint(fitobject);
            confR(i, :) = c(:, 1);
        else
            learnR(i) = NaN;
            confR(i, :) = [NaN, NaN];
        end

    end

end