function [Composition] = makeH2O(densitySet,varargin)
%MAKEH2O Calculate atomic density of the molecule H2O
% Description:  For macroscopic cross sections, calculate atomic density.
%
% USE:  [Composition] = makeH2O(densitySet,varargin)
%
% NOTES: Density set must be chosen, to use room temperature or 600 K
% densities. Called by makeGeometry, or used to define geometry.
%
% EXAMPLES:
%
% MAJOR UPDATES:
%   version  date     NetID   description
%   1.0      20110520 cld72   cleaned up and formatted
%
% FUTURE UPDATES:
%
% DEPENDENCIES:
%   macroconsts_300K.mat
%

isMicToMac = strcmp(varargin{1},'macro');
if ~isMicToMac
    error('Optional argument can only be ''macro''.')
end

load macroconsts_300K.mat

Composition = 1;

if isMicToMac

    % calculate molecular mass
    molecMass = 2*atomicMass.H1 + atomicMass.O16;

    if densitySet == 2
        density.H2O = .720; % VBUDS, or T = 600 K
    end

    % calculate atomic densities
    Composition = density.H2O * Navogadro / molecMass * 1e-24;

end

end
