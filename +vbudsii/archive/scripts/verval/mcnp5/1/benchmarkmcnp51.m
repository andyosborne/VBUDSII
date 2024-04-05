
fineGroupDef = 10.^(-4:.1:7)*1e-6;
fewGroupDef = [1e-4 1 100e3 1e7];
p = struct('nFineGroups',length(fineGroupDef)-1,... % leave blank.
           'nFewGroups',length(fewGroupDef)-1,...
           'fineGroupDef',fineGroupDef,...
           'fewGroupDef',fewGroupDef,...
           'temp',600,...
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
[ao, wo, N_H2O] = matl([1001 2;
                        8016 1], 1, density_H2O);
[ao, wo, N_UO2] = matl([92235 enrichment;
                        92238 1-enrichment;
                        8016 2], 1, density_UO2);

p.uc = uc;


fid = fopen('benchmarkmcnp51','w');
fprintf(fid,'validate vbudsii 111016\n');

do2D = 0;

temperature = p.temp*8.617e-11;

fprintf(fid,'C cell cards\n');
fprintf(fid,'C --- universe 1 --- \n');
% density, material density stuff, pitch, diameter, etc.
fprintf(fid,'10 100 -%.4f -1 u=1 imp:n=1 vol=%.4e tmp=%.4e\n',...
    density_UO2, .25*pi*p.uc.pinDiam^2, temperature);
fprintf(fid,'20 200 -%.4f 1 u=1 imp:n=1 vol=%.4e tmp=%.4e\n',...
    density_H2O, p.uc.pinPitch^2 - .25*pi*p.uc.pinDiam^2, temperature);

fprintf(fid,'C --- real world ---\n');
if do2D
    fprintf(fid,'30 0 -2 fill=1 u=2 imp:n=1\n');
    fprintf(fid,'40 0 -3 : 4 u=2 imp:n=0\n');
    % unsure about how this line below actually works
    fprintf(fid,'50 0 -2fill=2 lat=1 imp:n=1\n');
else
    fprintf(fid,'60 0 -2 fill=1 lat=1 imp:n=1\n');
end

fprintf(fid,'\nC surface cards\n');
fprintf(fid,'1 cz %.4f\n',p.uc.pinDiam);
fprintf(fid,'2 rpp -%.4f %.4f -%.4f %.4f 0. 0.\n',...
    .5*p.uc.pinPitch, .5*p.uc.pinPitch, .5*p.uc.pinPitch, .5*p.uc.pinPitch);

rxrheight = 100;
if do2D
    fprintf(fid,'3 pz -%.4f\',.5*rxrheight);
    fprintf(fid,'4 pz %.4f\',.5*rxrheight);
end

fprintf(fid,'\nC data cards\n');
fprintf(fid,'kcode 10000 1.0 150 800\n');
fprintf(fid,'ksrc 0.0 0.0 0.0\n');
fprintf(fid,'m100 92235 %.4e   92238 %.4e\n',N_UO2(1),N_UO2(2));
fprintf(fid,'     8016  %.4e\n',N_UO2(3));
fprintf(fid,'m200 1001 2 8016 1\n');
fprintf(fid,'mt200 lwtr\n');
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
fprintf(fid,'FM34 -1 100 (-6 -8)\n');

fprintf(fid,'FC44 fluxcapacitor44\n');
fprintf(fid,'F44:n 20\n');
fprintf(fid,'FM44 -1 200 (-6 -8)\n');

fclose(fid);
