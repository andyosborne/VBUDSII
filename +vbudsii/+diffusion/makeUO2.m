function [Composition] = makeUO2(densitySet,enrichment,varargin)
%MAKEUO2 Calculate atomic density or atom ratio of each isotope in UO2
% Description:  For macroscopic cross sections, calculate atomic densities.
% For microscopic cross sections, calculate atom ratios (e.g. 2/3 for H in H2O).
%
% USE:  [Composition] = makeUO2(densitySet,varargin)
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

%can output either macroscopic or microscopic
%micmac is optional

% if we want mictomac, Composition atomic densities are in units of
% #/(barn-cm).

isMicToMac = strcmp(varargin{1},'macro');
if ~isMicToMac
    disp('warning: UO2 cross section is not ''macro''.')
end

load macroconsts_300K.mat

Composition = [enrichment; (1-enrichment); 2];

if isMicToMac
    
    % molecular mass of UO2
    molecMass = enrichment*atomicMass.U235 + (1-enrichment)*atomicMass.U238 + 2*atomicMass.O16;
    
    if densitySet == 2
        density.UO2 = 11; % VBUDS, or T = 600 K
    end
    
    % atomic density
    Composition = density.UO2 * Navogadro / molecMass * Composition * 1e-24;
    
end

%% TRASH
% 
% if (nargin < 4) || (nargin > 5)
%     error('Either 4 or 5 inputs are required')
% end
% 
% optionIn = strcmp(varargin{1},'macro');
% 
% if ~optionIn
%     error('Optional argument can only be ''macro''.')
% end
% 
% % just micro combine
% xsOut = enrichment*U235xs + (1-enrichment)*U238xs + 2*O16xs;
% 
% 
% if optionIn
% 
%     xsOut = 
% 
% end
% 
% 
% end