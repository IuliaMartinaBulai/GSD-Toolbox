%--------------------------------------------------------------------------
% File: filters_plot.m
%
% Goal:
%   Plot and visualize the three filter functions used in the paper
%   (exponential, polynomial, and tanh-based filters), and automatically
%   save the figures in EPS and PDF formats.
%
% Use:
%   Run this script directly to generate the filter plots.
%   The parameter L can be adjusted to change the filter scale.
%
% Inputs:
%   None (all parameters are set inside the script).
%
% Outputs:
%   - EPS and PDF figures saved in the folder specified by `outdir`.
%   - One figure per filter:
%       * f_exp.eps / f_exp.pdf
%       * G_poly.eps / G_poly.pdf
%       * F_tanh.eps / F_tanh.pdf
%
% Recalls:
%   None 
%
% Authors:
%   Iulia Martina Bulai,
%   Isidoros Iakovidis,
%   Tim Steger
%
% Date last modified:
%   13-12-2025
%%-------------------------------------------------------------------------

clear;
clc;

L = 2;                               % set your L
x = linspace(0, 1.5*L, 4000);        % wide enough domain

% Functions
y_f = exp(- (x./(0.6*L)).^4);
y_G = max(0, 1 - (x./L).^5);
y_F = 0.5 * (1 - tanh((x - L)./0.4));

% Where to save
outdir = 'G:\.shortcut-targets-by-id\1fhidSzHlhR0M4Z2b4M-Q-llfyvf71Agj\Isidoros Iakovidis\Paper Graph Spectrum Dimension\Codes\results';
if ~exist(outdir, 'dir'); mkdir(outdir); end

names  = {'f_exp','G_poly','F_tanh'};   % distinct filenames
titles = { ...
    'f(x) = exp(-(x/(0.6L))^4)', ...
    'G(x) = max(0, 1 - (x/L)^5)', ...
    'F(x) = 0.5(1 - tanh((x - L)/0.4))' ...
    };

ys = {y_f, y_G, y_F};

for i = 1:3
    y = ys{i};
    idx_last = find(y > 1e-4, 1, 'last');
    if isempty(idx_last), idx_last = numel(y); end
    xx = x(1:idx_last);
    yy = y(1:idx_last);

    figure('Color','w');
    plot(xx, yy, 'LineWidth', 2);
    grid on; box on;
    xlabel('x');
    title(titles{i}, 'Interpreter','none');
    ylim([0 1.05]);
    xlim([0 max(xx)]);
    set(gca,'FontSize',12);
    set(gcf,'PaperPositionMode','auto');
    drawnow;

    base = fullfile(outdir, names{i});
    print(gcf, [base '.eps'], '-depsc', '-painters');
    print(gcf, [base '.pdf'], '-dpdf',  '-painters');
end

