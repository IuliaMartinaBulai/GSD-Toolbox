%--------------------------------------------------------------------------
% File: get_dataset.m
%
% Goal:
%   Return the sparse adjacency matrix A for one of the datasets used in the
%   numerical experiments of the paper (synthetic graphs and benchmark data).
%   For some datasets, coordinates are also exported to the base workspace
%   for visualization in the main script.
%
% Use:
%   A = get_dataset(dataset_name)
%   where dataset_name is one of:
%     'ring', 'path', 'swiss_roll', 'minnesota', 'bunny'
%
% Inputs:
%   dataset_name : (char/string) dataset identifier (case-insensitive).
%
% Outputs:
%   A : sparse adjacency matrix 
%
% Recalls:
%   - pdist2 (Statistics and Machine Learning Toolbox) for 'swiss_roll'
%   - load('minnesota.mat') for 'minnesota'
%   - load('data_bunny.mat') for 'bunny'
%
% Required data files:
%   - minnesota.mat : must contain variable A or W (sparse or dense)
%   - data_bunny.mat: must contain variable bunny (point cloud)
%
% Side effects (for plotting convenience):
%   - Assigns swiss_roll_coords (N-by-3) into base workspace for 'swiss_roll'
%   - Assigns bunny_coords (N-by-2) into base workspace for 'bunny'
%
% Authors:
%   Iulia Martina Bulai,
%   Isidoros Iakovidis,
%   Tim Steger
%
% Date last modified:
%   13-12-2025
%%-------------------------------------------------------------------------

function A = get_dataset(dataset_name)

    dataset_name = lower(string(dataset_name));

    % Dataset sizes (match the paper/experiments)
    % ring: 3000, path: 4000, swiss_roll: 3000
   
    switch dataset_name

        case "ring"
            % Synthetic ring (cycle) graph on N nodes (wrap-around).
            N = 3000;
            B = ones(N, 2);
            A = spdiags(B, [1 -1], N, N);
            A(1, end) = 1;
            A(end, 1) = 1;
            A = sparse(A);

        case "path"
            % Synthetic path (line) graph on N nodes (no wrap-around).
            N = 4000;
            B = ones(N, 2);
            A = spdiags(B, [1 -1], N, N);
            A(1, end) = 0;
            A(end, 1) = 0;
            A = sparse(A);

        case "swiss_roll"
            % Swiss roll manifold: 2D surface embedded in 3D, k-NN graph.
            N = 3000;
            fprintf('Generating Swiss roll with N=%d points...\n', N);

            % Parameters: t ∈ [1.5π, 4.5π], h ∈ [0, 1]
            t = 1.5*pi + 3*pi * rand(N, 1);
            h = rand(N, 1);

            % 3D embedding: (x, y, z) = (t cos(t), 21 h, t sin(t))
            X = zeros(N, 3);
            X(:, 1) = t .* cos(t);
            X(:, 2) = 21 * h;
            X(:, 3) = t .* sin(t);

            % k-nearest-neighbor graph (mutual k-NN)
            k = 8;
            fprintf('  Building %d-NN graph...\n', k);

            D = pdist2(X, X);

            A = sparse(N, N);
            for i = 1:N
                [~, idx] = sort(D(i, :));
                neighbors = idx(2:(k+1)); % exclude self
                A(i, neighbors) = 1;
            end

            A = double((A + A') > 0);
            A = sparse(A);

            % Export coordinates for visualization in main scripts
            assignin('base', 'swiss_roll_coords', X);

            fprintf('  Swiss roll generated: %d nodes, %d edges\n', ...
                    N, nnz(A)/2);

        case "bunny"
            % Stanford bunny graph (GBFlearn-style preprocessing + r-ball graph).
            fprintf('Loading GBF bunny from data_bunny.mat...\n');

            S = load('data_bunny.mat');
            if ~isfield(S,'bunny')
                error('data_bunny.mat must contain variable ''bunny''.');
            end

            bunny = S.bunny;                        % (Nbunny x 3)
            nodes = [bunny(:,1), bunny(:,2)];       % 2D projection

            % 1) Remove very close points (same logic as GBF_gengraph)
            thresh = 0.0005;
            N0 = size(nodes,1);
            idx = [];        % indices to remove
            stp = 0;
            for i = 1:N0
                for j = i+1:N0
                    if norm(nodes(i,:) - nodes(j,:), 2) <= thresh
                        stp = stp + 1;
                        idx(stp) = j;
                    end
                end
            end

            idx = unique(idx);
            if ~isempty(idx)
                nodes(idx,:) = [];
            end

            % 2) Radius graph: connect all pairs with distance <= r
            r = 0.01;
            N = size(nodes,1);
            fprintf('  Building r-ball graph on bunny: N = %d, r = %.4f...\n', N, r);

            A = sparse(N,N);
            for i = 1:N
                for j = i+1:N
                    if norm(nodes(i,:) - nodes(j,:), 2) <= r
                        A(i,j) = 1;
                        A(j,i) = 1;
                    end
                end
            end

            A = sparse(A);

            % Store coordinates for your plotting code
            assignin('base','bunny_coords',nodes);

            fprintf('  Bunny graph generated: %d nodes, %d edges\n', ...
                    N, nnz(A)/2);

        case "minnesota"
            % Load Minnesota graph
            S = load('minnesota.mat');
            if isfield(S, 'A')
                A = S.A;
            elseif isfield(S, 'W')
                A = S.W;
            else
                error('minnesota.mat must contain variable A or W.');
            end
            A = sparse(A);

        otherwise
            error(['Unknown dataset name: %s. ' ...
                   'Valid options: ''ring'', ''path'', ''swiss_roll'', ''minnesota'', ''bunny''.'], ...
                  dataset_name);
    end
end
