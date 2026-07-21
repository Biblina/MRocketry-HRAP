% HRAP Simulation Environment - CoolProp Properties Used
%
% Program:  NOX (CoolProp-backed replacement)
%
% Purpose:  Drop-in replacement for HRAP's empirical saturated nitrous oxide
%           correlations. Returns the same struct fields and units, but the
%           values are derived from CoolProp's Helmholtz-energy equation of
%           state (Lemmon & Span) rather than independent curve fits.
%
%           Properties are tabulated from CoolProp on first call and cached
%           to disk, so the simulation loop performs fast interpolation
%           instead of thousands of MATLAB->Python calls.
%
% Requires: Python environment configured through pyenv, with CoolProp installed
%           (pip install coolprop)
%
% Returned fields (units match the original HRAP NOX.m):
%   op.Pv - vapor pressure [Pa]
%   op.rho_l - saturated liquid density [kg/m^3]
%   op.rho_v - saturated vapor density [kg/m^3]
%   op.Hv - latent heat of vaporization [kJ/kg]
%   op.Cp - saturated liquid specific heat [kJ/(kg*K)]
%   op.Z - saturated vapor compressibility [-]
function [op] = NOX(T)

persistent tbl

if isempty(tbl)
    cachefile = fullfile(fileparts(mfilename('fullpath')),'NOX_coolprop_cache.mat');
    if isfile(cachefile)
        S = load(cachefile);
        tbl = S.tbl;
    else
        fprintf('Building CoolProp property tables for nitrous oxide...\n');
        tbl = build_tables('NitrousOxide');
        save(cachefile,'tbl');
        fprintf('Done. Cached to %s\n', cachefile);
    end
end

% Warn if called outside the tabulated range (silent extrapolation is how
% the original correlations hid its errors)
if T < tbl.T(1) || T > tbl.T(end)
    warning('NOX:OutOfRange','T = %.2f K is outside tabulated range [%.2f, %.2f] K.', T, tbl.T(1), tbl.T(end));
end

op.Pv = interp1(tbl.T, tbl.Pv, T, 'pchip');
op.rho_l = interp1(tbl.T, tbl.rho_l, T, 'pchip');
op.rho_v = interp1(tbl.T, tbl.rho_v, T, 'pchip');
op.Hv = interp1(tbl.T, tbl.Hv, T, 'pchip');
op.Cp = interp1(tbl.T, tbl.Cp, T, 'pchip');
op.Z = interp1(tbl.T, tbl.Z, T, 'pchip');

end

% Table generation - queries CoolProp once over a temperature grid
function tbl = build_tables(fluid)

CP = @(out,T,Q) double(py.CoolProp.CoolProp.PropsSI(out,'T',T,'Q',Q,fluid));

Tc = double(py.CoolProp.CoolProp.PropsSI('Tcrit',fluid)); % ~309.52 K

% Range: ~-90 C (matching the original fits' lower bound) up to just short
% of the critical point, where saturation properties become ill-conditioned
tbl.T = linspace(183, Tc - 0.15, 600);
n = numel(tbl.T);

[tbl.Pv, tbl.rho_l, tbl.rho_v, tbl.Hv, tbl.Cp, tbl.Z] = deal(zeros(1,n));

for i = 1:n
    Ti = tbl.T(i);
    tbl.Pv(i) = CP('P',Ti,0); % Pa
    tbl.rho_l(i) = CP('D',Ti,0); % kg/m^3
    tbl.rho_v(i) = CP('D',Ti,1); % kg/m^3
    tbl.Hv(i) = (CP('H',Ti,1) - CP('H',Ti,0))/1000; % J/kg -> kJ/kg
    tbl.Cp(i) = CP('C',Ti,0)/1000; % J/(kg*K) -> kJ/(kg*K)
    tbl.Z(i) = CP('Z',Ti,1); % dimensionless
end

tbl.fluid = fluid;
tbl.Tcrit = Tc;
tbl.created = datetime('now');

end