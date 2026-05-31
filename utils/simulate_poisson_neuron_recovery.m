
function simSignal = simulate_poisson_neuron_recovery(unitSignal, binT, np_fs, w_t, Nsim)

    frateW = 0.001;
    simSignal = zeros(size(unitSignal, 1), size(unitSignal, 2), Nsim, 'logical');

    signal_smooth = single(smoothdata(unitSignal, 1, 'movmean', frateW*np_fs))*np_fs;
    r = smoothdata(signal_smooth, 2, "gaussian", 5*binT);
    
    rVals = 0:1:1000;
    qVals = zeros(size(rVals));
    qVals(1) = 0;
    qVals(2:end) = min(1./integrate(w_t, 1./rVals(2:end), 0), 1000*rVals(2:end));

    q = zeros(size(r));
    Ni = size(r, 1);
    Nj = size(r, 2);
    parfor i = 1:Ni
        for j = 1:Nj
            q(i, j) = interp1(rVals, qVals, r(i, j));
        end
    end

    Nj = size(unitSignal, 2);
    Nt = size(unitSignal, 1);

    for i = 1:Nsim
        parfor j = 1:Nj
            thisSim = zeros(size(unitSignal, 1), 1, 'logical');
            tSpike = 1;
            ttest = 2;
            while ttest<=Nt
                alpha = -log(rand);
                found = false;
                ttest = tSpike+1;
                intOut = 0;
                while ~found
                    intOut = intOut+q(ttest, j)*w_t((ttest-tSpike)/np_fs)/np_fs;
                    if intOut>alpha
                        thisSim(ttest) = true;
                        tSpike = ttest;
                        found = true;
                    else
                        if ttest==Nt
                            found = true;
                        end
                    end
                    ttest=ttest+1;
                end
            end
            simSignal(:, j, i) = thisSim;
        end
    end

end
