
function [simSignal, burstTrials] = simulate_burst_ramp(fRate, Nsim, np_fs, cW, w_t, burstFR)

    deltaF = fRate(end)-fRate(1);
    f_poisson = 0;
    
    f_burst = 1/(burstFR - fRate(1));
   
    simSignal = zeros(round(cW*np_fs), numel(fRate), Nsim, 'logical');

    rVals = 0:1:1000;
    qVals = zeros(size(rVals));
    qVals(1) = 0;
    qVals(2:end) = min(1./integrate(w_t, 1./rVals(2:end), 0), 1000*rVals(2:end));


    q = zeros(size(fRate));
    Ni = numel(fRate);
    parfor i = 1:Ni
        q(i) = interp1(rVals, qVals, fRate(1)+(fRate(i)-fRate(1))*f_poisson);
    end
    qBurst = interp1(rVals, qVals, burstFR);

    pBurst = (fRate-fRate(1)) * f_burst;

    Nj = numel(fRate);
    Nt = round(cW*np_fs);

    burstTrials = zeros(numel(fRate), Nsim, 'logical');
    parfor k = 1:Nsim
        thisSim = zeros(round(cW*np_fs), Nj, 'logical');
        for j = 1:Nj
            tSpike = -0.005*rand*np_fs;
            ttest = 1;
            isBursting = rand<pBurst(j);
            if isBursting
                q_j = qBurst;
                burstTrials(j, k) = true;
            else
                q_j = q(j);
            end
            while ttest<=Nt
                alpha = -log(rand);
                found = false;
                intOut = 0;
                while ~found
                    intOut = intOut+q_j*w_t((ttest-tSpike)/np_fs)/np_fs;
                    if intOut>alpha
                        thisSim(ttest, j) = true;
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
        simSignal(:, :, k) = thisSim;
    end

end
