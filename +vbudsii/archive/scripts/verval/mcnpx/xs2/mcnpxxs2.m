
fineGroupDef = 10.^(-4:.1:7)*1e-6;
fewGroupDef = [1e-4 1 100e3 1e7];
p = struct('nFineGroups',length(fineGroupDef)-1,... % leave blank.
           'nFewGroups',length(fewGroupDef)-1,...
           'fineGroupDef',fineGroupDef,...
           'fewGroupDef',fewGroupDef,...
           'temp',1200,...
           'powerDensity',50,... % MW/m^3
           'nTimeSteps',1,...
           'makeLibrary',1,...
           'makeLibraryTempFlag',2,...
           'XSLibraryMAT',fullfile('..','..','data','XSLibrary.mat'),...
           'verbose',1,...
           'resolveXS',1,...
           'S0iterthresh',.00001); % this largely affects the validation effort.


uc.pinDiam = 1; % [1,1];
uc.pinPitch = 2; % [2,2];
uc.f = 1; % some weighting of the fuel regions
uc.g = 1; % some weighting of the moderator regions
uc.sauerConst.mod = 2.35;
uc.sauerConst.fuel = 5.00;

enrichment = .032;
density_H2O = 0.72; % g/cm^3
density_UO2 = 11; %g/cm^3
[ao, wo, N_UO2] = matl([92235 enrichment;
                        92238 1-enrichment;
                        8016 2], 1, density_UO2);

p.uc = uc;


fid = fopen('mcnpxxs2','w');
fprintf(fid,'mcnpxxs2 111107\n');

temperature = p.temp*8.617e-11;

fprintf(fid,'C cell cards\n');
% density, material density stuff, pitch, diameter, etc.
fprintf(fid,'10 100 -%.4f -1 imp:n=1 vol=%.4e tmp=%.4e\n',...
    density_UO2, .25*pi*p.uc.pinDiam^2, temperature);
fprintf(fid,'20 0 1 imp:n=0\n');

fprintf(fid,'\nC surface cards\n');
fprintf(fid,'1 cz %.4f\n',p.uc.pinDiam/2);

fprintf(fid,'\nC data cards\n');
fprintf(fid,'kcode 1000 1.0 10 40\n');
fprintf(fid,'ksrc 0.0 0.0 0.0\n');
fprintf(fid,'m100 92235 1\n');
% tallies
fprintf(fid,'print 40 128 140\n');
fprintf(fid,'E0');
for i = 1:22
    fprintf(fid,'     %.4e  %.4e  %.4e  %.4e  %.4e\n',...
                fineGroupDef((i-1)*5+2),fineGroupDef((i-1)*5+3),...
                fineGroupDef((i-1)*5+4),fineGroupDef((i-1)*5+5),...
                fineGroupDef((i-1)*5+6));
end

fprintf(fid,'FC14 fluxcapacitor14\n');
fprintf(fid,'F14:n 10\n');

fprintf(fid,'FC34 fluxcapacitor34\n');
fprintf(fid,'F34:n 10\n');
fprintf(fid,'FM34 1 100 (1) (2) (18) (102) (-6 -7)\n');

fclose(fid);
