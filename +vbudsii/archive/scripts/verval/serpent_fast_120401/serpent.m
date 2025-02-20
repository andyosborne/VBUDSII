clear all;
runMCNPX = 1;
runSerpent = 1;
doPDFXS = 1;
plotS0Bounds = 1;
plotVBUDSIflux = 1;


mycd = cd;
cd(fullfile('..','..','..',''));
vbudsiiDir = cd;
cd(mycd); % alternatively, mfilename

[pL, p, g, vbudsiflux] = fast(vbudsiiDir);

%% RUN SERPENT


if runSerpent
    ifname = 'serpentcompare';
    serpentexecpath = fullfile('serpent');
    writeSerpentSpectralInput(p, g, serpentexecpath, ifname);
    cd(serpentexecpath);
    system(['sss ' ifname]);
    cd(mycd);
end

%% RUN MCNPX

ifname = 'serpentcompare';
mcnpxexecpath = fullfile('mcnpx');

if runMCNPX
    % ~exist([mcnpxexecpath filesep 'outp'], 'file')
    writeMCNPXSpectralInput(p, g, mcnpxexecpath, ifname);
    cd(mcnpxexecpath);
    system(['mcnpx i= ' ifname]);
    if ispc
        system('del runtpe');
    else
        system('rm runtpe');
    end
    cd(mycd);
end


[Tally1 TallyP] = TallyPull2([mcnpxexecpath filesep 'outp']);
save XSTallyserpentcompare Tally1 TallyP

%cd(fullfile('..','..',''));
%addpath(fullfile('..','..',''));
addpath(p.vbudsiiDir);

%% RUN VBUDSII
[Results, p , g, L] = main(p,g);
r = Results.Region(1);


%% TEMP FOR SERPENT
%tempserpentplotter


% Must be assigned before dividing out the flux.
mcnpxPower = TallyP('14.100.-6.-8');
mcnpxPowerDensity = mcnpxPower(end);

TallyP0 = TallyP;
MKEYS = keys(TallyP);
for i = 1:length(MKEYS)
    if strcmp(MKEYS{i}(1:2), '14') && ~strcmp(MKEYS{i}, '14.0.0')
        TallyP(MKEYS{i}) = TallyP(MKEYS{i})./TallyP('14.0.0');
        aaa = TallyP(MKEYS{i});
        TallyP(MKEYS{i}) = aaa(1:end-1);
    elseif strcmp(MKEYS{i}(1:2), '24') && ~strcmp(MKEYS{i}, '24.0.0')
        TallyP(MKEYS{i}) = TallyP(MKEYS{i})./TallyP('24.0.0');
        aaa = TallyP(MKEYS{i});
        TallyP(MKEYS{i}) = aaa(1:end-1);
    end
end

tallyxs = Tally1;

astring = 'No MCNPX cross section intervention. Base case.';

%% RUN VBUDSII AGAIN

bstring = ['Using MCNPX cross sections everywhere, ' ...
    ', but no constant scaling of water xs''s.'];
p2 = p;
p2.useMCNPXTallyXS = TallyP;
[Results, p2 , g2, L2] = main(p2,g);
r2 = Results.Region(1);

%cstring = ['Using water xs scaling factor from MCNPX ' ...
%    'cross sections, but no other intervention from MCNPX'];
%p3 = p;
%p3.MCNPXwaterCorrectionFactor = 1.888;
%[Results, p3 , g3, L3] = main(p3,g);
%r3 = Results.Region(1);
%
%dstring = ['Using MCNPX cross sections and the water xs scaling ' ...
%    'correction factor from MCNPX'];
%p4 = p;
%p4.useMCNPXTallyXS = TallyP;
%p4.MCNPXwaterCorrectionFactor = 1.888;
%[Results, p4 , g4, L4] = main(p4,g);
%r4 = Results.Region(1);

%disp('RETURN AT THIS LINE');
%return;
cd(mycd);

%% COMPARE TO SERPENT RESULTS
serpentreportcompare(['home/fitze/Dropbox/UTA/r/code/vbudsii/' ...
    'scripts/verval/serpent_fast_120401/serpentreport'],...
    'serpent_fast_120401', ...
    'serpentcompare_res', ...
    'serpentcompare_xs', ...
    'serpentcompare_det', ...
    Tally1, ...
    TallyP, ...
    mcnpxPowerDensity, ...
    vbudsiflux, ...
                       p, ...
                       g, ...
                       L, ...
                       r, ...
                       p2, ...
                       g2, ...
                       L2, ...
                       r2, ...
                       astring, ...
                       bstring);

% Create report!
% vbudsiireport('/home/fitze/Dropbox/UTA/r/code/vbudsii/scripts/verval/serpentcompare/reportfunc', ...
%                        'serpentcompare', ...
%                        Tally1, ...
%                        TallyP, ...
%                        mcnpxPowerDensity, ...
%                        vbudsiflux, ...
%                        p2, ...
%                        g2, ...
%                        L2, ...
%                        r2)

vbudsiireportcompare('/home/fitze/Dropbox/UTA/r/code/vbudsii/scripts/verval/serpentcompare/compare_ab', ...
                       'serpentcompare', ...
                       Tally1, ...
                       TallyP, ...
                       mcnpxPowerDensity, ...
                       vbudsiflux, ...
                       p, ...
                       g, ...
                       L, ...
                       r, ...
                       p2, ...
                       g2, ...
                       L2, ...
                       r2, ...
                       astring, ...
                       bstring);

vbudsiireportcompare('/home/fitze/Dropbox/UTA/r/code/vbudsii/scripts/verval/serpentcompare/compare_ac', ...
                       'serpentcompare', ...
                       Tally1, ...
                       TallyP, ...
                       mcnpxPowerDensity, ...
                       vbudsiflux, ...
                       p, ...
                       g, ...
                       L, ...
                       r, ...
                       p3, ...
                       g3, ...
                       L3, ...
                       r3, ...
                       astring, ...
                       cstring);

vbudsiireportcompare('/home/fitze/Dropbox/UTA/r/code/vbudsii/scripts/verval/serpentcompare/compare_bc', ...
                       'serpentcompare', ...
                       Tally1, ...
                       TallyP, ...
                       mcnpxPowerDensity, ...
                       vbudsiflux, ...
                       p2, ...
                       g2, ...
                       L2, ...
                       r2, ...
                       p3, ...
                       g3, ...
                       L3, ...
                       r3, ...
                       bstring, ...
                       cstring);

vbudsiireportcompare('/home/fitze/Dropbox/UTA/r/code/vbudsii/scripts/verval/serpentcompare/compare_bd', ...
                       'serpentcompare', ...
                       Tally1, ...
                       TallyP, ...
                       mcnpxPowerDensity, ...
                       vbudsiflux, ...
                       p2, ...
                       g2, ...
                       L2, ...
                       r2, ...
                       p4, ...
                       g4, ...
                       L4, ...
                       r4, ...
                       bstring, ...
                       dstring);

vbudsiireportcompare('/home/fitze/Dropbox/UTA/r/code/vbudsii/scripts/verval/serpentcompare/compare_cd', ...
                       'serpentcompare', ...
                       Tally1, ...
                       TallyP, ...
                       mcnpxPowerDensity, ...
                       vbudsiflux, ...
                       p3, ...
                       g3, ...
                       L3, ...
                       r3, ...
                       p4, ...
                       g4, ...
                       L4, ...
                       r4, ...
                       cstring, ...
                       dstring);
