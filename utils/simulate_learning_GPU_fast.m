function learnRpredicted = simulate_learning_GPU_fast(eta, sigma, phi, unitCOM)


    distVals = gpuArray(zeros(size(unitCOM, 1), 'single'));
    for i = 1:3
        distVals = distVals + (repmat(unitCOM(:, i), 1, size(unitCOM, 1)) - repmat(unitCOM(:, i)', size(unitCOM, 1), 1)).^2;
    end
    distVals = sqrt(distVals);
    
    learnRpredicted = eta*sum(repmat(phi, 1, numel(phi)).*exp(-(distVals.^2)/(2*sigma^2)), 1);

end