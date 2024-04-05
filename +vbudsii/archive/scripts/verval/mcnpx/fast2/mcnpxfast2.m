
fineGroupDef = 10.^(-4:.1:7)*1e-6;
fewGroupDef = [1e-4 1 100e3 1e7];
p = struct('nFineGroups',length(fineGroupDef)-1,... % leave blank.
           'nFewGroups',length(fewGroupDef)-1,...
           'fineGroupDef',fineGroupDef,...
           'fewGroupDef',fewGroupDef,...
           'powerDensity',50,... % MW/m^3
           'nTimeSteps',1,...
           'makeLibrary',1,...
           'makeLibraryTempFlag',2,...
           'XSLibraryMAT',fullfile('..','..','data','XSLibrary.mat'),...
           'verbose',1,...
           'resolveXS',1,...
           'S0iterthresh',.00001); % this largely affects the validation effort.

% Determining pin diameter.
uc.pinPitch = 2; % [2,2];

fuelVolFraction = 0.37;
uc.pinDiam = sqrt(uc.pinPitch^2 * fuelVolFraction * 4 / pi);

uc.f = 1; % some weighting of the fuel regions
uc.g = 1; % some weighting of the moderator regions
uc.sauerConst.mod = 2.35;
uc.sauerConst.fuel = 5.00;

p.uc = uc;

enrichment = .256;
density_Na = 0.882; % g/cm^3
density_UO2 = 11; %g/cm^3
[ao, wo, N_Na] = matl([11023 1], 1, density_Na);
[ao, wo, N_UO2] = matl([92235 enrichment;
                        92238 1-enrichment;
                        8016 2], 1, density_UO2);



fid = fopen('mcnpxfast2','w');
fprintf(fid,'mcnpxfast2 111128\n');

do2D = 0;

tempconversion = 8.617e-11;
tempFuel = (900 + 273)*tempconversion;
tempCoolant = (300 + 273)*tempconversion;

fprintf(fid,'C cell cards\n');
fprintf(fid,'C --- universe 1 --- \n');
% density, material density stuff, pitch, diameter, etc.
fprintf(fid,'10 100 -%.4f -1 u=1 imp:n=1 vol=%.4e tmp=%.4e\n',...
    density_UO2, .25*pi*p.uc.pinDiam^2, tempFuel);
fprintf(fid,'20 200 -%.4f 1 u=1 imp:n=1 vol=%.4e tmp=%.4e\n',...
    density_Na, p.uc.pinPitch^2 - .25*pi*p.uc.pinDiam^2, tempCoolant);

fprintf(fid,'C --- real world ---\n');
if do2D
    fprintf(fid,'30 0 -2 fill=1 u=2 imp:n=1\n');
    fprintf(fid,'40 0 -3 : 4 u=2 imp:n=0\n');
    % unsure about how this line below actually works
    fprintf(fid,'50 0 -2 fill=2 lat=1 imp:n=1\n');
else
    fprintf(fid,'60 0 -2 fill=1 lat=1 imp:n=1\n');
end

fprintf(fid,'\nC surface cards\n');
fprintf(fid,'1 cz %.4f\n',p.uc.pinDiam/2);
fprintf(fid,'2 rpp -%.4f %.4f -%.4f %.4f 0. 0.\n',...
    .5*p.uc.pinPitch, .5*p.uc.pinPitch, .5*p.uc.pinPitch, .5*p.uc.pinPitch);

rxrheight = 100;
if do2D
    fprintf(fid,'3 pz -%.4f\',.5*rxrheight);
    fprintf(fid,'4 pz %.4f\',.5*rxrheight);
end

fprintf(fid,'\nC data cards\n');
fprintf(fid,'kcode 1000 1.0 10 40\n');
fprintf(fid,'ksrc 0.0 0.0 0.0\n');
fprintf(fid,'m100 92235 %.4e   92238 %.4e\n',N_UO2(1),N_UO2(2));
fprintf(fid,'     8016  %.4e\n',N_UO2(3));
fprintf(fid,'m200 11023 1\n');

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

fprintf(fid,'FC24 fluxcapacitor24\n');
fprintf(fid,'F24:n 20\n');

fprintf(fid,'FC34 fluxcapacitor34\n');
fprintf(fid,'F34:n 10\n');
fprintf(fid,'FM34 -1 100 (1) (2) (18) (102) (-6 -7)\n');

fprintf(fid,'FC44 fluxcapacitor44\n');
fprintf(fid,'F44:n 20\n');
fprintf(fid,'FM44 -1 200 (1) (2) (18) (102) (-6 -7)\n');

fclose(fid);
