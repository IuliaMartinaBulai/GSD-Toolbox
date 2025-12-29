%--------------------------------------------------------------------------
% File: brain_demo.m
%
% Goal:
%   Estimate the spectral dimension of the brain graph using a stochastic
%   approximation of the eigenvalue counting function (no full eigendecomposition),
%   then fit a Weyl-type law on an automatic fraction-of-spectrum window.
%   The script also produces and saves two core figures:
%     (1) N_G(\ell) vs \ell
%     (2) log(N_G(\ell)) vs log(\ell) with fitted line
%
% Use:
%   Place 'Normalized_Laplacian.mat' (containing variable L) in the current folder
%   and run this script.
%
% Inputs:
%   - Normalized_Laplacian.mat : must contain variable L (sparse normalized Laplacian)
%
% Outputs:
%   In folder ./brain_results:
%     - brain_N_vs_L.pdf
%     - brain_loglog_fit.pdf
%     - brain_results_YYYYMMDD_HHMMSS.mat  (struct "results")
%     - brain_summary.txt
%
%
% Notes:
%   This script uses SGWT scripts. 
%   - sgwt_rough_lmax
%   - sgwt_cheby_coeff

% Authors:
%   Iulia Martina Bulai,
%   Isidoros Iakovidis,
%   Tim Steger
%
% Date last modified:
%   13-12-2025
%%-------------------------------------------------------------------------

%  BRAIN GRAPH SPECTRAL ANALYSIS (Stochastic-Only, No Eigenvalues)
%  Uses percentage-based automatic fitting

clc;
clear all;
close all;

addpath("Data");
addpath("Utils");

fprintf('BRAIN GRAPH SPECTRAL DIMENSION ANALYSIS \n');
fprintf('Start time: %s\n\n', datestr(now));

%  1 Load Graph
fprintf(' STEP 1: Loading Graph \n');
file_name = 'Normalized_Laplacian.mat';

tic;
load(file_name, 'L');
load_time = toc;

N = size(L, 1);
nnz_L = nnz(L);
density = nnz_L / (N^2);

fprintf('Loaded: %s (%.2f seconds)\n', file_name, load_time);
fprintf('Graph: N = %d nodes, %d edges\n', N, nnz_L/2);
fprintf('Density: %.2e, Memory: %.2f GB\n', density, (nnz_L * 16) / 1e9);

% Quick validation
fprintf('\nValidating Laplacian...\n');
d_sample = full(diag(L(1:min(1000,N), 1:min(1000,N))));
fprintf('  Diagonal: min=%.4f, median=%.4f, max=%.4f\n', ...
        min(d_sample), median(d_sample), max(d_sample));

[~, ~, vals] = find(L);
off_diag = vals(1:min(10000, length(vals)));  % Sample
frac_neg = sum(off_diag < 0) / length(off_diag);
fprintf('  Off-diagonal: %.1f%% negative (expect ~100%%)\n', frac_neg * 100);

clear vals off_diag d_sample;

%   2 : max eigenvalue to apply the theorems
fprintf('\n--- STEP 2: Estimating λ_max ---\n');
tic;
lmax = sgwt_rough_lmax(L);
fprintf('λ_max ≈ %.6f (%.2f seconds)\n', lmax, toc);

%  STEP 3: Stochastic and parameter setting
fprintf('\n STEP 3: Stochastic Spectral Scan \n');

% Parameters (optimized for brain-scale graph)
m_trace = 30;           % Chebyshev order
K_probes = 100;         % Hutchinson probes (balance speed/accuracy)
num_L_points = 80;      % Scan resolution
sigma_rel = 0.04;       % Filter width

% Define scan range
L_vals = linspace(0.01, lmax + 0.1, num_L_points);

% Normalize Laplacian for Chebyshev
L_tilde = (2*L - lmax*speye(N)) / lmax;

% Storage
NL_approx = zeros(size(L_vals));
scan_times = zeros(size(L_vals));

fprintf('Running scan (this will take a while)...\n');
total_tic = tic;

for i = 1:length(L_vals)
    if mod(i, 10) == 0 || i == 1
        elapsed = toc(total_tic);
        eta = (elapsed / i) * (length(L_vals) - i);
        fprintf('  %d/%d (%.1f%%) - ETA: %.1f min\n', ...
                i, length(L_vals), 100*i/length(L_vals), eta/60);
    end

    tic;
    Lthr = L_vals(i);
    sigma = max(sigma_rel * Lthr, 1e-8);

    % Soft-step filter
    g = @(x) 0.5 * (1 - tanh((x - Lthr) ./ sigma));

    % Chebyshev coefficients
    c = sgwt_cheby_coeff(g, m_trace, 2*(m_trace+1), [0, lmax]);
    c(1) = c(1) / 2;

    % Hutchinson trace estimation
    trace_est = 0;
    for k = 1:K_probes
        z = 2*(rand(N,1) > 0.5) - 1;
        T0 = z;
        T1 = L_tilde * z;

        yv = c(1)*T0;
        if m_trace >= 1, yv = yv + c(2)*T1; end

        for j = 2:m_trace
            Tn = 2*(L_tilde*T1) - T0;
            yv = yv + c(j+1)*Tn;
            T0 = T1; T1 = Tn;
        end

        trace_est = trace_est + (z' * yv);
    end

    NL_approx(i) = trace_est / K_probes;
    scan_times(i) = toc;
end

total_time = toc(total_tic);
fprintf('Scan complete: %.1f minutes (%.2f sec/point avg)\n', ...
        total_time/60, mean(scan_times));

%  STEP 4: Automatic Weyl Law Fitting
fprintf('\n--- STEP 4: Weyl Law Fitting ---\n');

% Fraction-of-spectrum window
alpha = 0.04;   % lower fraction of spectrum 
beta  = 0.70;   % upper fraction of spectrum 
% For reporting in percent
pct_low  = alpha * 100;
pct_high = beta  * 100;

% Estimated fraction N(λ)/N
frac_approx = NL_approx / N;   % in [0,1]

% Select fitting range based on N(λ)/N
mask_fit = (frac_approx >= alpha) & ...
           (frac_approx <= beta)  & ...
           (NL_approx > 0)        & ...
           (L_vals > 0)           & ...
           isfinite(L_vals)       & ...
           isfinite(NL_approx);

fprintf('Fraction range for fit: [%.2f, %.2f] of N\n', alpha, beta);
fprintf('Points in mask: %d\n', sum(mask_fit));

if sum(mask_fit) < 5
    error('Weyl fit: too few points in fraction-based mask. Increase num_L_points or relax alpha/beta.');
end

% Log–log fit on masked region
x_fit = log(L_vals(mask_fit));
y_fit = log(NL_approx(mask_fit));
valid = isfinite(x_fit) & isfinite(y_fit);

x_fit = x_fit(valid);
y_fit = y_fit(valid);

% Linear fit in log–log
p = polyfit(x_fit, y_fit, 1);

% Extract dimension and constant
d_est = 2 * p(1);
C_est = exp(p(2));

fprintf('\n=== WEYL LAW RESULTS ===\n');
fprintf('Fitting range in λ: [%.4f, %.4f]\n', ...
        min(L_vals(mask_fit)), max(L_vals(mask_fit)));
fprintf('Fitting points used: %d\n', numel(x_fit));
fprintf('Estimated dimension: d = %.4f\n', d_est);
fprintf('Weyl constant:       C = %.1f\n', C_est);
fprintf('Formula: N(λ) ≈ %.1f * λ^{%.3f}\n', C_est, d_est/2);

fprintf('\nINTERPRETATION (very rough):\n');
if d_est < 1.5
    fprintf('  → Low-dimensional structure (chain/tree-like)\n');
elseif d_est < 2.5
    fprintf('  → 2D-like structure (surface/cortical sheet)\n');
elseif d_est < 3.5
    fprintf('  → 3D-like structure (volumetric organization)\n');
else
    fprintf('  → High-dimensional (complex connectivity)\n');
end

% L-range and Weyl curve only on fit band
L_fit      = L_vals(mask_fit);
N_weyl_fit = C_est * (L_fit .^ (d_est/2));  % N(λ) ≈ C λ^{d/2}

%   5: Visualization
fprintf('\nSTEP 5: Generating Core Figures \n');

output_dir = fullfile(pwd, 'brain_results');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

%  Figure 1: N_G(ℓ) vs ℓ
fig1 = figure('Position', [100, 100, 900, 600]);
hold on; grid on;

plot(L_vals, NL_approx, 'b-', 'LineWidth', 2, ...
     'DisplayName', '$N_G(\ell)$ estimate');

xlabel('$\ell$','Interpreter','latex');
ylabel('$N_G(\ell)$','Interpreter','latex');
title(sprintf('Brain Graph Spectral Counting (N=%d)', N));
legend('Location', 'northwest', 'FontSize', 11, 'Interpreter','latex');
xlim([0, lmax]);
ylim([0, N]);
set(gca,'TickLabelInterpreter','latex');

exportgraphics(fig1, fullfile(output_dir, 'brain_N_vs_L.pdf'), ...
               'ContentType', 'vector');
fprintf('  Saved: brain_N_vs_L.pdf\n');

%  Figure 2: log(N_G(ℓ)) vs log(ℓ) with fitted line
fig2 = figure('Position', [150, 150, 900, 600]);
hold on; grid on;

x_all = log(L_vals);
y_all = log(NL_approx + eps);
valid_all = isfinite(x_all) & isfinite(y_all) & (NL_approx > 0);

x_line = x_all(valid_all);
y_line = y_all(valid_all);

plot(x_line, y_line, 'b-', 'LineWidth', 2, ...
     'DisplayName', '$\log(N_G(\ell))$');

xx = linspace(min(x_fit), max(x_fit), 200);
plot(xx, polyval(p, xx), 'k--', 'LineWidth', 3, ...
     'DisplayName', sprintf('Fit from approx. $N_G(\\ell)$: $d=%.2f$', d_est));

xlabel('$\log(\ell)$','Interpreter','latex');
ylabel('$\log(N_G(\ell))$','Interpreter','latex');
title('Weyl Law: log--log Scale');
legend('Location', 'southeast', 'FontSize', 11, 'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex');

exportgraphics(fig2, fullfile(output_dir, 'brain_loglog_fit.pdf'), ...
               'ContentType', 'vector');
fprintf('  Saved: brain_loglog_fit.pdf\n');

%  STEP 6: Save Results
fprintf('\n STEP 6: Saving Results \n');

results = struct();
results.graph.N       = N;
results.graph.nnz     = nnz_L;
results.graph.density = density;
results.graph.lmax    = lmax;

results.scan.L_vals     = L_vals;
results.scan.NL_approx  = NL_approx;
results.scan.times      = scan_times;
results.scan.total_time = total_time;
results.scan.params     = struct('m', m_trace, 'K', K_probes, ...
                                 'sigma_rel', sigma_rel, ...
                                 'num_L_points', num_L_points);

results.weyl.d            = d_est;
results.weyl.C            = C_est;
results.weyl.alpha_beta   = [alpha, beta];
results.weyl.fit_L_range  = [min(L_fit), max(L_fit)];
results.weyl.mask_fit     = mask_fit;

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
save(fullfile(output_dir, sprintf('brain_results_%s.mat', timestamp)), ...
     'results', '-v7.3');

% Text summary
fid = fopen(fullfile(output_dir, 'brain_summary.txt'), 'w');
fprintf(fid, '=== BRAIN GRAPH SPECTRAL DIMENSION ANALYSIS ===\n\n');
fprintf(fid, 'Date: %s\n', datestr(now));
fprintf(fid, 'Graph: N=%d nodes, %d edges\n\n', N, nnz_L/2);
fprintf(fid, 'WEYL LAW RESULTS:\n');
fprintf(fid, '  Dimension: d = %.4f\n', d_est);
fprintf(fid, '  Constant:  C = %.1f\n', C_est);
fprintf(fid, '  Fit fraction [alpha, beta]: [%.2f, %.2f]\n', alpha, beta);
fprintf(fid, '  Fit λ-range: [%.4f, %.4f]\n', min(L_fit), max(L_fit));
fprintf(fid, '\nCOMPUTATION:\n');
fprintf(fid, '  Total time: %.1f minutes\n', total_time/60);
fprintf(fid, '  Parameters: m=%d, K=%d, num_L_points=%d\n', ...
        m_trace, K_probes, num_L_points);
fclose(fid);

fprintf('  Saved: brain_results_%s.mat\n', timestamp);
fprintf('  Saved: brain_summary.txt\n');

fprintf('\n=== ANALYSIS COMPLETE ===\n');
fprintf('Spectral dimension of brain network: d = %.3f\n', d_est);
fprintf('All results saved in: %s\n', output_dir);
