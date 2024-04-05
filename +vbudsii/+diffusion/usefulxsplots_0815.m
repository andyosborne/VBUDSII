% creates cross section plots for important isotopes and cross sections

% plotxsn(macro,mytitle,legends,xbelow,varargin)

% 1: total
% 2: scatter
% 18: fission
% 102: gamma (capture)

load reducedlib

% H-1
g1 = plotxsn(0,'H-1',{'t','a','s'},'T = 600 BW = inf',...
    H1_3_1_600,...
    H1_3_102_600,...
    sum(H1_0_2_600)');

print('-dpsc','-r300','/home/fitze/Dropbox/UTA/r/doc/usefulxsplots_H1')
print('-dpng','-r300','/home/fitze/Dropbox/UTA/r/doc/usefulxsplots_H1')

% O-16
g2 =plotxsn(0,'O-16',{'t','a','s'},'T = 600 BW = inf',...
    O16_3_1_600,...
    O16_3_102_600,...
    sum(O16_0_2_600)');

print('-dpsc','-r300','/home/fitze/Dropbox/UTA/r/doc/usefulxsplots_O16')
print('-dpng','-r300','/home/fitze/Dropbox/UTA/r/doc/usefulxsplots_O16')

% H2O
g3 = plotxsn(0,'H2O',{'t','a','s'},'T = 600 BW = inf',...
    H2O_3_1_600,...
    H2O_3_102_600,...
    sum(H2O_0_2_600)');

print('-dpsc','-r300','/home/fitze/Dropbox/UTA/r/doc/usefulxsplots_H2O')
print('-dpng','-r300','/home/fitze/Dropbox/UTA/r/doc/usefulxsplots_H2O')

% U-235
g4 = plotxsn(0,'U-235',{'t','a','f','s'},'T = 600 BW = inf',...
    U235_3_1_600,...
    U235_3_102_600+U235_3_18_600,...
    U235_3_18_600,...
    sum(U235_0_2_600)');

print('-dpsc','-r300','/home/fitze/Dropbox/UTA/r/doc/usefulxsplots_U235')
print('-dpng','-r300','/home/fitze/Dropbox/UTA/r/doc/usefulxsplots_U235')

% U-238
g5 = plotxsn(0,'U-238',{'t','a','f','s'},'T = 600 BW = inf',...
    U238_3_1_600,...
    U238_3_102_600+U238_3_102_600,...
    U238_3_18_600,...
    sum(U238_0_2_600)');

print('-dpsc','-r300','/home/fitze/Dropbox/UTA/r/doc/usefulxsplots_U238')
print('-dpng','-r300','/home/fitze/Dropbox/UTA/r/doc/usefulxsplots_U238')
