function [communities, Q] = simple_louvain(corr_matrix, threshold)
    % Build adjacency from correlation
    %A = abs(corr_matrix) > threshold;
    %A = A - diag(diag(A));
    
    A = corr_matrix;

    % Initialize each node in its own community
    N = size(A, 1);
    communities = 1:N;
    
    % Calculate modularity
    m = sum(A(:)) / 2;
    k = sum(A, 2);
    
    improved = true;
    while improved
        improved = false;
        for i = 1:N
            best_comm = communities(i);
            best_dQ = 0;
            
            % Try moving node i to each neighbor's community
            neighbors = find(A(i, :));
            neighbor_comms = unique(communities(neighbors));
            
            for c = neighbor_comms
                dQ = modularity_gain(A, communities, i, c, k, m);
                if dQ > best_dQ
                    best_dQ = dQ;
                    best_comm = c;
                    improved = true;
                end
            end
            
            communities(i) = best_comm;
        end
    end
    
    % Relabel communities sequentially
    [~, ~, communities] = unique(communities);
    
    % Calculate final modularity
    %Q = compute_modularity(A, communities, k, m);
    Q = 0;
end

function dQ = modularity_gain(A, communities, node, new_comm, k, m)
    old_comm = communities(node);
    if old_comm == new_comm
        dQ = 0;
        return;
    end
    
    % Edges to new community
    ki_in_new = sum(A(node, communities == new_comm));
    sigma_new = sum(k(communities == new_comm));
    
    % Edges to old community
    ki_in_old = sum(A(node, communities == old_comm)) - A(node, node);
    sigma_old = sum(k(communities == old_comm));
    
    dQ = (ki_in_new - ki_in_old) / (2*m) - ...
         k(node) * (sigma_new - sigma_old + k(node)) / (2*m)^2;
end