%--------------------------------------------------------------------------
% File: graph_spectral_dimension_main.m
%
% Goal:
%   Diagnostic script for estimating the spectral (Weyl) dimension of a graph.
%   The script computes exact eigenvalues (when feasible), performs stochastic
%   eigenvalue counting, fits a Weyl law, and produces visualization figures.
%
% Use:
%   Select a dataset name below and run the script.
%
% Inputs:
%   - Dataset specified by variable 'dataset'
%
% Outputs:
%   - Eigenvalue histogram
%   - Graph visualization
%   - N_G(\ell) vs \ell
%   - log(N_G(\ell)) vs log(\ell) with Weyl fit
%   - EPS and PDF figures saved in ./results/<dataset>/
%
%  Recalls: sgwt_rough_lmax, sgwt_cheby_coeff
%  Authors:
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

addpath("Data");
addpath("Utils");

%% Parameters
m_trace   = 60;
q_probes  = 200;
sigma_rel = 0.04;
nL        = 500;

% Filter family: 'tanh' 
filter_type = 'tanh';

% Dataset: 'ring' | 'path' | 'minnesota' | 'swiss_roll' | 'bunny'
dataset = 'bunny';

% Avoid underscores in figure titles
ds_pretty = regexprep(dataset, '_', ' ');
ds_pretty = regexprep(ds_pretty, '(^.)', '${upper($1)}');

%% Load adjacency
A = get_dataset(dataset);
N = size(A,1);

%% Normalized Laplacian
L = sgwt_laplacian(A, 'opt', 'normalized');

%% Exact eigenvalues and histogram
fprintf('\nDiagnostic 1: Computing exact eigenvalues (ground truth)\n');
tic;
Lam_exact = sort(eig(full(L)));

lmx  = max(Lam_exact);
Lhat = (2*L - lmx*speye(N)) / lmx;

t_eig = toc;
fprintf('Eigenvalue computation: %.2f seconds\n', t_eig);
fprintf('λ_min = %.6f\n', min(Lam_exact(Lam_exact > 1e-17)));
fprintf('λ_max = %.6f\n', max(Lam_exact));

% Closed-form theory (ring/path only)
Lam_theory = [];
C_theory   = NaN;
k_vals = 0:(N-1);
switch lower(dataset)
    case 'ring'
        Lam_theory = sort(2 * sin(pi * k_vals / N).^2).';
        C_theory   = (sqrt(2) * N) / pi;
    case 'path'
        Lam_theory = sort(1 - cos(pi * k_vals / (N-1))).';
        C_theory   = (sqrt(2) * (N-1)) / pi;
end

if ~isempty(Lam_theory)
    fprintf('Max difference from theory: %.2e\n', max(abs(Lam_exact - Lam_theory)));
end

figure(1); clf;
%set(gcf,'Position',default_pos);
histogram(Lam_exact, 60, 'Normalization','pdf');
grid on;
set(gca,'FontSize',16);
xlabel('\lambda','FontSize',16);
ylabel('density','FontSize',16);
title(sprintf('Eigenvalue Histogram — %s (N=%d)', ds_pretty, N), 'FontSize',16);
legend('Eigenvalue density','Location','northwest','FontSize',14);
fh_hist = gcf;

%% Graph visualization
A_vis = [];
coords_vis = [];
N_vis = [];
graph_title = '';
is_swiss3d = false;
coords3d = [];

switch lower(dataset)
    case 'ring'
        N_vis = 30;
        B_vis = ones(N_vis,2);
        A_vis = spdiags(B_vis,[1 -1],N_vis,N_vis);
        A_vis(1,end) = 1;
        A_vis(end,1) = 1;
        A_vis = sparse(A_vis);
        t_vis = (0:N_vis-1)' * (2*pi/N_vis);
        coords_vis = [cos(t_vis), sin(t_vis)];
        graph_title = 'Ring graph';

    case 'path'
        N_vis = 30;
        B_vis = ones(N_vis,2);
        A_vis = spdiags(B_vis,[1 -1],N_vis,N_vis);
        A_vis(1,end) = 0;
        A_vis(end,1) = 0;
        A_vis = sparse(A_vis);
        coords_vis = [(1:N_vis).', zeros(N_vis,1)];
        graph_title = 'Path graph';

    case 'swiss_roll'
        A_vis = A;
        N_vis = N;
        graph_title = 'Swiss roll graph';
        if evalin('base','exist(''swiss_roll_coords'',''var'')')
            coords3d = evalin('base','swiss_roll_coords');
        end
        is_swiss3d = true;

    case 'bunny'
        A_vis = A;
        N_vis = N;
        graph_title = 'Bunny graph';
        if evalin('base','exist(''bunny_coords'',''var'')')
            coords_vis = evalin('base','bunny_coords');
        end

    case 'minnesota'
        A_vis = A;
        N_vis = N;
        graph_title = 'Minnesota graph';
        S = load('minnesota.mat');
        if isfield(S,'xy')
            coords_vis = S.xy;
        elseif isfield(S,'coords')
            coords_vis = S.coords;
        end
end

figure(2); clf;
%set(gcf,'Position', default_pos);

if is_swiss3d && ~isempty(coords3d)
    scatter3(coords3d(:,1), coords3d(:,2), coords3d(:,3), ...
        8, coords3d(:,2), 'filled', 'MarkerEdgeAlpha', 0.3);
    colormap('jet');
    grid on;
    axis equal;
    view(45, 30);
    set(gca,'FontSize',16);
    xlabel('x'); ylabel('height'); zlabel('z');
    title(graph_title,'FontSize',16);
    axis tight
else
    G_vis = graph(A_vis);
    if ~isempty(coords_vis)
        plot(G_vis,'XData',coords_vis(:,1),'YData',coords_vis(:,2), ...
            'Marker','.', 'MarkerSize',8, 'NodeColor',[0 0 0], ...
            'EdgeAlpha',0.3, 'LineWidth',1, 'NodeLabel',{});
        axis equal off;
    else
        plot(G_vis,'Layout','force','Iterations',100, ...
            'Marker','.', 'MarkerSize',6, 'NodeColor',[0 0 0], ...
            'EdgeAlpha',0.25, 'LineWidth',0.8, 'NodeLabel',{});
        axis off;
    end
    set(gca,'FontSize',16);
    title(graph_title,'FontSize',16);
end

fh_graph = gcf;

%% Thresholds and exact counting
L_vals = linspace(1e-6, lmx, nL);
N_exact = arrayfun(@(Lt) sum(Lam_exact <= Lt), L_vals);


%% Plot of the selected filter for the user

Lthr_plot = 0.7;


% Build the exact same filter handle used in the stochastic scan
g_plot = select_filter(Lthr_plot, sigma_rel, filter_type);

% Plot g(lambda) on the spectral interval
npts = 2000;
xg = linspace(0, lmx, npts);
yg = g_plot(xg);

fh_filter = figure; clf;


plot(xg, yg, 'LineWidth', 2);
grid on;
set(gca,'FontSize',16);
xlabel('\lambda','FontSize',16);
ylabel('g(\lambda)','FontSize',16);
title(sprintf('Filter used: %s (L_{thr}=%.3f, \\sigma_{rel}=%.3f)', ...
    filter_type, Lthr_plot, sigma_rel), ...
    'FontSize',16, 'Interpreter','none');
ylim([-0.05, 1.05]);
xlim([0, lmx]);






%% Filter sanity check
test_L = 0.5;
g_test = select_filter(test_L, sigma_rel, filter_type);
fw = g_test(Lam_exact);
ab = sum(Lam_exact <= test_L);

%% Stochastic scan
Nq = 2*(m_trace + 1);
NL_approx = zeros(size(L_vals));

for i = 1:numel(L_vals)
    Lthr = L_vals(i);
    g = select_filter(Lthr, sigma_rel, filter_type);
    c = sgwt_cheby_coeff(g, m_trace, Nq, [0, lmx]);
    c(1) = c(1)/2;

    acc = 0;
    for k = 1:q_probes
        z = 2*(rand(N,1) > 0.5) - 1;
        T0 = z; T1 = Lhat*z;
        yv = c(1)*T0;
        if m_trace >= 1, yv = yv + c(2)*T1; end
        for j = 2:m_trace
            Tn = 2*(Lhat*T1) - T0;
            yv = yv + c(j+1)*Tn;
            T0 = T1; T1 = Tn;
        end
        acc = acc + (z'*yv);
    end
    NL_approx(i) = acc / q_probes;
end

%% Weyl fit
% eignvalue threshold (adjust according to your data set)
alpha = 0.04;
beta = 0.4;
frac_approx = NL_approx / N;
mask_final = (frac_approx >= alpha) & (frac_approx <= beta) & (NL_approx > 0);

x_fit_final = log(L_vals(mask_final));
y_fit_final = log(NL_approx(mask_final) + eps);
valid_final = isfinite(x_fit_final) & isfinite(y_fit_final);

p_final = polyfit(x_fit_final, y_fit_final, 1);
d_final = 2*p_final(1);
C_final = exp(p_final(2));




%% Results Summary
fprintf('\n results \n');
if isnan(C_theory)
    fprintf('Estimated (stochastic): d = %.3f, C = %.1f\n', d_final, C_final);
else
    fprintf('Theoretical: d = 1.000, C = %.1f\n', C_theory);
    fprintf('Estimated (stochastic): d = %.3f, C = %.1f\n', d_final, C_final);
    fprintf('Relative errors: d: %.1f%%, C: %.1f%%\n', ...
        100*abs(d_final - 1.0)/1.0, 100*abs(C_final - C_theory)/C_theory);
end



%% Standalone N_G(ℓ) vs ℓ

x_exact = log(L_vals);
y_exact = log(N_exact + eps);
valid_exact = isfinite(x_exact) & isfinite(y_exact) & (N_exact > 0);
fh_NL = figure(4); clf;
%set(gcf,'Position',default_pos);
hold on; grid on; set(gca,'FontSize',16);

plot(L_vals, N_exact, 'k-', 'LineWidth', 2.5, ...
    'DisplayName', 'Exact count');

plot(L_vals, NL_approx, 'r-', 'LineWidth', 1.5, ...
    'DisplayName', 'Stochastic approx.');

xlabel('$\ell$','FontSize',16,'Interpreter','latex');
ylabel('$N_G(\ell)$','FontSize',16,'Interpreter','latex');
title(sprintf('%s (N=%d): Eigenvalue Counting', ds_pretty, N), ...
    'FontSize',16, 'Interpreter', 'none');
legend('Location','northwest','FontSize',14,'Interpreter','latex');
xlim([0, 2]);
ylim([0, N]);
set(gca,'TickLabelInterpreter','latex');

%% Standalone Weyl log–log fit
fh_WeylSingle = figure(5); clf;
%set(gcf,'Position',default_pos);
hold on; grid on; set(gca,'FontSize',16);

plot(x_exact(valid_exact), y_exact(valid_exact), ...
    'k-', 'LineWidth', 2, ...
    'DisplayName', '$\log N_G(\ell)$');

xx = linspace(min(x_fit_final(valid_final)), max(x_fit_final(valid_final)), 200);
plot(xx, polyval(p_final, xx), ...
    'b-', 'LineWidth', 2.5, ...
    'DisplayName', sprintf('Fit from approx. $N_G(\\ell)$: $d=%.2f$, $C=%.0f$', d_final, C_final));

xlabel('$\log(\ell)$','FontSize',16,'Interpreter','latex');
ylabel('$\log(N_G(\ell))$','FontSize',16,'Interpreter','latex');
title('Weyl Law: log-log Scale','FontSize',16);
legend('Location','northwest','FontSize',14,'Box','on','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex');


%% Auto-save figures (EPS + PDF, per-dataset subfolder; overwrite)
basedir = 'G:\.shortcut-targets-by-id\1fhidSzHlhR0M4Z2b4M-Q-llfyvf71Agj\Isidoros Iakovidis\Paper Graph Spectrum Dimension\Codes';
results_root = fullfile(basedir, 'Results');
if ~exist(results_root, 'dir')
    mkdir(results_root);
end

ds_name = lower(dataset);
outdir  = fullfile(results_root, ds_name);
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

% Ensure vector export renderer
set(fh_hist,   'Renderer', 'painters');
set(fh_graph,  'Renderer', 'painters');
set(fh_filter,   'Renderer', 'painters');
set(fh_NL,     'Renderer', 'painters');
set(fh_WeylSingle, 'Renderer', 'painters');

% Base names inside the dataset subfolder
b_hist  = 'eigenvalue_histogram';
b_graph = 'graph';
b_whole = 'NL_and_Weylfit';   % combined 2x1 figure (optional)
b_ax1   = 'NL_vs_L';          % standalone N(L) vs L
b_ax2   = 'Weyl_loglog_fit';  % standalone Weyl log-log

% Save all as full figures (no axes-only export)
save_eps_pdf(fh_hist, fullfile(outdir, b_hist));
save_eps_pdf(fh_graph, fullfile(outdir, b_graph));
save_eps_pdf(fh_filter, fullfile(outdir, b_whole));
save_eps_pdf(fh_NL, fullfile(outdir, b_ax1));
save_eps_pdf(fh_WeylSingle,  fullfile(outdir, b_ax2));

if exist('fh_swiss_3d', 'var')
    save_eps_pdf(fh_swiss_3d, fullfile(outdir, 'swiss_roll_3D'));
end

fprintf('\nSaved EPS and PDF to: %s\n', outdir);

%% Local helper function for export
function save_eps_pdf(hfig, basepath)
ok1 = false;
ok2 = false;
try
    exportgraphics(hfig, [basepath '.eps'], 'ContentType','vector');
    ok1 = true;
catch
end
try
    exportgraphics(hfig, [basepath '.pdf'], 'ContentType','vector');
    ok2 = true;
catch
end
if ~(ok1 && ok2)
    try
        saveas(hfig, basepath, 'epsc');
    end
    try
        saveas(hfig, basepath, 'pdf');
    end
end
end
