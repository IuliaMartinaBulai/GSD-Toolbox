%--------------------------------------------------------------------------
% File: brain_data_manipulation.m
%
% Goal:
%   Load the latest brain_results_*.mat file and regenerate the two figures:
%     (1) N_G(\ell) vs \ell
%     (2) log(N_G(\ell)) vs log(\ell) with Weyl fit
%
% Use:
%   - Set results_dir to the folder containing brain_results_*.mat
%   - Run the script to export EPS/PDF figures to outdir
%
% Inputs:
%   - brain_results_*.mat (variable 'results')
%
% Outputs:
%   - NL_vs_L.pdf / NL_vs_L.eps
%   - Weyl_loglog_fit.pdf / Weyl_loglog_fit.eps
%
% Recalls:
%   None.
%
% Authors:
%   Iulia Martina Bulai,
%   Isidoros Iakovidis,
%   Tim Steger
%
% Date last modified:
%   13-12-2025
%%-------------------------------------------------------------------------

clc;
clear;
close all;

% Paths
results_dir = 'G:\.shortcut-targets-by-id\1fhidSzHlhR0M4Z2b4M-Q-llfyvf71Agj\Isidoros Iakovidis\Paper Graph Spectrum Dimension\Codes\brain_results';

if ~exist(results_dir, 'dir')
    error('Results directory not found: %s', results_dir);
end

files = dir(fullfile(results_dir, 'brain_results_*.mat'));
if isempty(files)
    error('No brain_results_*.mat files found in %s', results_dir);
end

[~, idx] = max([files.datenum]);
results_file = files(idx).name;

fprintf('Loading results from: %s\n', fullfile(results_dir, results_file));
load(fullfile(results_dir, results_file), 'results');

% Unpack
N        = results.graph.N;
nnz_L    = results.graph.nnz;
lmax     = results.graph.lmax;

L_vals    = results.scan.L_vals;
NL_approx = results.scan.NL_approx;

d_est     = results.weyl.d;
C_est     = results.weyl.C;
mask_fit  = results.weyl.mask_fit;

fprintf('Graph: N = %d nodes, %d edges\n', N, nnz_L/2);
fprintf('Estimated dimension d = %.4f, C = %.1f\n', d_est, C_est);

% Output folder
fig_w = 900;
fig_h = 600;
default_pos = [100 100 fig_w fig_h];

base_draft   = 'G:\.shortcut-targets-by-id\1fhidSzHlhR0M4Z2b4M-Q-llfyvf71Agj\Isidoros Iakovidis\Paper Graph Spectrum Dimension\Draft';
results_root = fullfile(base_draft, 'Results');
outdir       = fullfile(results_root, 'brain_results');

if ~exist(results_root, 'dir'), mkdir(results_root); end
if ~exist(outdir,       'dir'), mkdir(outdir);       end

fprintf('Output directory: %s\n', outdir);

% Figure 1: N_G(ℓ) vs ℓ
fig1 = figure;
set(fig1, 'Position', default_pos);
hold on; grid on;

plot(L_vals, NL_approx, 'r-', 'LineWidth', 1.5, ...
     'DisplayName', 'Stochastic approx.');

xlabel('$\ell$', 'FontSize', 16, 'Interpreter','latex');
ylabel('$N_G(\ell)$', 'FontSize', 16, 'Interpreter','latex');
title(sprintf('Brain graph spectral counting (N=%d)', N), 'FontSize', 16);
legend('Location', 'northwest', 'FontSize', 14);

xlim([0, lmax]);
ylim([0, N]);

set(gca, 'FontSize', 16, 'TickLabelInterpreter','latex');
set(fig1, 'Renderer', 'painters');

pdf_path1 = fullfile(outdir, 'NL_vs_L.pdf');
eps_path1 = fullfile(outdir, 'NL_vs_L.eps');

exportgraphics(fig1, pdf_path1, 'ContentType', 'vector');
exportgraphics(fig1, eps_path1, 'ContentType', 'vector');

fprintf('Saved: %s\n', pdf_path1);
fprintf('Saved: %s\n', eps_path1);

% Figure 2: log N_G(ℓ) vs log ℓ with Weyl fit
fig2 = figure;
set(fig2, 'Position', default_pos);
hold on; grid on;

x_all = log(L_vals);
y_all = log(NL_approx + eps);
valid_all = isfinite(x_all) & isfinite(y_all) & (NL_approx > 0);

x_line = x_all(valid_all);
y_line = y_all(valid_all);

plot(x_line, y_line, 'k-', 'LineWidth', 2, ...
     'DisplayName', '$\log N_G(\ell)$');

L_fit      = L_vals(mask_fit & valid_all);
x_fit      = log(L_fit);
x_fit_min  = min(x_fit);
x_fit_max  = max(x_fit);

xx = linspace(x_fit_min, x_fit_max, 200);
yy = log(C_est) + (d_est/2)*xx;

plot(xx, yy, 'b-', 'LineWidth', 2.5, ...
     'DisplayName', sprintf('Fit from approx. $N_G(\\ell)$: $d=%.2f$, $C=%.0f$', d_est, C_est));

xlabel('$\log(\ell)$', 'FontSize', 16, 'Interpreter','latex');
ylabel('$\log(N_G(\ell))$', 'FontSize', 16, 'Interpreter','latex');
title('Brain graph: Weyl law log--log fit', 'FontSize', 16);
legend('Location', 'northwest', 'FontSize', 14, 'Interpreter','latex');

set(gca, 'FontSize', 16, 'TickLabelInterpreter','latex');
set(fig2, 'Renderer', 'painters');

pdf_path2 = fullfile(outdir, 'Weyl_loglog_fit.pdf');
eps_path2 = fullfile(outdir, 'Weyl_loglog_fit.eps');

exportgraphics(fig2, pdf_path2, 'ContentType', 'vector');
exportgraphics(fig2, eps_path2, 'ContentType', 'vector');

fprintf('Saved: %s\n', pdf_path2);
fprintf('Saved: %s\n', eps_path2);
