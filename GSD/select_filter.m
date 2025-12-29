%--------------------------------------------------------------------------
% File: select_filter.m
%
% Goal:
%   Return a filter function handle g(x) used to approximate the spectral
%   projector in the stochastic eigenvalue counting procedure.
%
% Use:
%   g = select_filter(Lthr, sigma_rel, type)
%
% Inputs:
%   Lthr      : threshold eigenvalue (scalar)
%   sigma_rel : relative width parameter (used for tanh filter)
%   type      : filter type, one of {'tanh'}
%
% Outputs:
%   g : function handle defining the selected filter
%
% Definitions:
%   - Tanh filter:
%       F(x) = 0.5 * (1 - tanh((x - Lthr) / sigma))
%
% Recalls:
%   None (standalone utility function).
%
% Authors:
%   Iulia Martina Bulai,
%   Isidoros Iakovidis,
%   Tim Steger
%
% Date last modified:
%   13-12-2025
%%-------------------------------------------------------------------------

function g = select_filter(Lthr, sigma_rel, type)


    if nargin < 3
        error('Usage: g = select_filter(Lthr, sigma_rel, type)');
    end

    switch lower(type)
        case 'tanh'
            sigma = sigma_rel * Lthr;
            g = @(x) 0.5 * (1 - tanh((x - Lthr) ./ sigma));

        otherwise
            error('Unknown filter type: %s (choose tanh)', type);
    end
end
