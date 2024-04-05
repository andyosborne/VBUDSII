function xs5 = makeTransport2(xs18,xs102,xs2,xs251)
%MAKETRANSPORT2 Calculate a single transport cross section vector given
%scattering kernel, fission, gamma, and mubar.
% Description:  Calculate a single transport cross section, for any energy group
% structure. 
%
% USE:  xs5 = makeTransport2(xs18,xs102,xs2,xs251)
%
% NOTES: Total cross section is calculated as sum of fission, gamma, and
% scattering Called by makeLibrary, so this is not called at "runtime"
% in the eventual code. That is, transport is treated as an "immutable"
% cross section.
%
% EXAMPLES: 
%
% MAJOR UPDATES:
%   version  date     NetID   description
%   1.0      20110520 cld72   cleaned up and formatted
%              
% FUTURE UPDATES:
%
% DEPENDENCIES: none
%

% total cross section
xs1 = xs18 + xs102 + sum(xs2)';

% transport
xs5 = xs1 - xs251.*sum(xs2)';

end