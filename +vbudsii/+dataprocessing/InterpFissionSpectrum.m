function [fissionSpectrum] = InterpFissionSpectrum(p)
%
%VBUDSII/DATAPROCESSING/INTERPFISSIONSPECTRUM Returns the fission spectrum at
% the energy group edges given in p.fineGroupDef.
%
% INPUT
%   groupDef       Fine group definition.
%
% OUTPUT
%   fissionSpectrum     [p.nFineGroups x 1 double] Fission neutron spectrum.

% NOTES
%   This data is needed to determine spectral neutron sources.
%
% MAJOR REVISIONS
%   date        handle      description
%   20111030    cld2469     writing comments
%
% TASKLIST
%   1- Inquire as to its incompleteness.
%   2- Find an accurate data source.
%   3- Determine an appropriate interpolation scheme.
%   4- INTEGRATED AREA MUST BE 1.

if isfield(p, 'useWatt') && p.useWatt
    avgGroupEnergy = p.fineGroupDef(1:end-1) + 1/2 * diff(p.fineGroupDef);
    notNormalized = wattspectrum(avgGroupEnergy') .* diff(p.fineGroupDef)';
    fissionSpectrum = notNormalized / sum(notNormalized);
    % try using .2 and 4!

    return;
end

if isfield(p, 'doctorKernel')
    disp('DOCTORING KERNEL <---------------');
    fissionSpectrum = zeros(110, 1);
    fissionSpectrum(p.doctorKernel(2)) = 1;
    return;
end

fissionSpectrumData = [       0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
   5.387e-007
  2.2731e-007
  3.2109e-007
  4.5356e-007
  6.4067e-007
  9.0497e-007
  1.2783e-006
  1.8056e-006
  2.5503e-006
  3.6023e-006
  5.0872e-006
  7.1858e-006
   1.015e-005
  1.4335e-005
  2.0247e-005
  2.8594e-005
  4.0381e-005
  5.7022e-005
  8.0515e-005
   0.00011368
   0.00016045
   0.00022648
    0.0003196
   0.00045089
   0.00063589
   0.00089638
    0.0012628
    0.0017777
    0.0024999
    0.0035108
    0.0049203
    0.0068803
    0.0095886
     0.013301
     0.018332
     0.025034
      0.03372
     0.044418
     0.056938
     0.071041
      0.08592
     0.099945
      0.11018
      0.11203
      0.10214
     0.081935
     0.056915
     0.033105
       0.0151
    0.0064333];

energyBounds = 10.^(-4:.1:7);
nGroups = length(energyBounds) - 1;

fissionSpectrum = zeros(p.nFineGroups, 1);

for iGroup = 1:p.nFineGroups
    iBounds = [min(nGroups, ...
               find(energyBounds > p.fineGroupDef(iGroup), 1, 'first')), ...
               min(nGroups, ...
               find(energyBounds <= p.fineGroupDef(iGroup+1), 1, 'last'))];
    fissionSpectrum(iGroup) = mean(fissionSpectrumData(iBounds));
%    fissionSpectrum = fissionSpectrum / sum(fissionSpectrum);
end

end


function f = wattspectrum(E, varargin)
E = E*1e-6;
    if length(varargin) == 2
        a = varargin{1};
        b = varargin{2};
    else
        a = 0.988;
        b = 2.249;
    end
C = sqrt(4/(pi*a^3 * b)) * exp(-a*b/4);
f = C *exp(-E/a) .* sinh((b*E).^(.5));
end
