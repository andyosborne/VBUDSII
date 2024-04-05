function writeMCNPXSpectralInput(p, g, mcnpxexecpath, ifname);
mycd = cd;
if ~exist(mcnpxexecpath, 'dir')
    mkdir(mcnpxexecpath)
end
cd(mcnpxexecpath);

if exist('outp', 'file')
    if ispc
        system('del outp');
    else
        system('rm outp');
    end
end

if exist('runtpe', 'file')
    if ispc
        system('del runtpe');
    else
        system('rm runtpe');
    end
end

if exist('srctp', 'file')
    if ispc
        system('del scrtp');
    else
        system('rm srctp');
    end
end

fid = fopen(ifname,'w');
fprintf(fid, '%s %s \n', ifname, datestr(now, 30));

do2D = 0;

% see section 5.4.4 of the MCNPX 2.7 manual for this conversion factor.
tempconversion = 8.617e-11;

fprintf(fid,'C cell cards\n');
fprintf(fid,'C --- universe 1 --- \n');
% density, material density stuff, pitch, diameter, etc.
% fuel
fprintf(fid,'10 100 -%.4f -1 u=1 imp:n=1 vol=%.4e tmp=%.4e\n',...
    g.regionDef(1).cellDef(1).initDensity,...
    .25 * pi * p.uc.pinDiam^2, ...
    g.regionDef(1).cellDef(1).initTemp * tempconversion...
    );
% coolant/moderator
fprintf(fid,'20 200 -%.4f 1 u=1 imp:n=1 vol=%.4e tmp=%.4e\n',...
    g.regionDef(1).cellDef(2).initDensity,...
    p.uc.pinPitch^2 - .25*pi*p.uc.pinDiam^2, ...
    g.regionDef(1).cellDef(2).initTemp * tempconversion...
    );

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
fprintf(fid,'kcode 2000 1.0 50 100\n');
fprintf(fid,'ksrc 0.0 0.0 0.0\n');
% pin cell
fprintf(fid,'m100\n');
for zaididx = 1:length(g.regionDef(1).cellDef(1).initZAIDs)
    fprintf(fid,'     %i%s %.4e\n',...
        g.regionDef(1).cellDef(1).initZAIDs(zaididx),...
        MCNPX_temp(g.regionDef(1).cellDef(1).initTemp, 'k2l'), ...
        g.regionDef(1).cellDef(1).initNumDensities(zaididx));
end
for zaididx = 1:length(g.regionDef(1).cellDef(1).initZAIDs)
    fprintf(fid,'m%i  %i%s 1\n',100+zaididx,g.regionDef(1).cellDef(1).initZAIDs(zaididx), MCNPX_temp(g.regionDef(1).cellDef(1).initTemp, 'k2l'));
end
% annular/anti-pin cell
fprintf(fid,'m200\n');
if g.regionDef(1).cellDef(2).initZAIDs == 222;
    g.regionDef(1).cellDef(2).initZAIDs = [1001 8016];
    N = g.regionDef(1).cellDef(2).initNumDensities;
    g.regionDef(1).cellDef(2).initNumDensities = [2*N N];
end
for zaididx = 1:length(g.regionDef(1).cellDef(2).initZAIDs)
    fprintf(fid,'     %i%s %.4e\n',...
        g.regionDef(1).cellDef(2).initZAIDs(zaididx),...
        MCNPX_temp(g.regionDef(1).cellDef(2).initTemp, 'k2l'), ...
        g.regionDef(1).cellDef(2).initNumDensities(zaididx));
end
if sum(g.regionDef(1).cellDef(2).initZAIDs == 222) || ...
    sum(g.regionDef(1).cellDef(2).initZAIDs == 1001)
    fprintf(fid, 'mt200 lwtr%s\n', ...
        MCNPX_temp(g.regionDef(1).cellDef(2).initTemp, 'k2l', 'lwtr'));
end
for zaididx = 1:length(g.regionDef(1).cellDef(2).initZAIDs)
    fprintf(fid,'m%i  %i%s 1\n',200+zaididx, ...
        g.regionDef(1).cellDef(2).initZAIDs(zaididx), ...
        MCNPX_temp(g.regionDef(1).cellDef(2).initTemp, 'k2l'));
    if g.regionDef(1).cellDef(2).initZAIDs(zaididx) == 1001
        fprintf(fid, 'mt%i lwtr%s\n',200+zaididx, ...
            MCNPX_temp(g.regionDef(1).cellDef(2).initTemp, 'k2l', 'lwtr'));
    end
end

% tallies
fprintf(fid,'print 40 128 140\n');
fprintf(fid,'E0\n');
p.fineGroupDef = 1e-6*p.fineGroupDef;
for i = 1:22
    fprintf(fid,'     %.4e  %.4e  %.4e  %.4e  %.4e\n',...
                p.fineGroupDef((i-1)*5+2), p.fineGroupDef((i-1)*5+3),...
                p.fineGroupDef((i-1)*5+4), p.fineGroupDef((i-1)*5+5),...
                p.fineGroupDef((i-1)*5+6));
end

fprintf(fid,'F14:n 10\n');

fprintf(fid,'F24:n 20\n');

fprintf(fid,'FM14 (1)\n');
fprintf(fid,'     (-1 100 (1) (2) (4) (18) (102) (-8) (-6 -7) (-6 -8))\n');

for zaididx = 1:length(g.regionDef(1).cellDef(1).initZAIDs)
    fprintf(fid,'     (1 %i (1) (2) (4) (18) (102) (-8) (-6 -7) (-6 -8))\n',100+zaididx);
end

fprintf(fid,'FM24 (1)\n');
fprintf(fid,'     (-1 200 (1) (2) (4) (18) (102) (-8) (-6 -7) (-6 -8))\n');

for zaididx = 1:length(g.regionDef(1).cellDef(2).initZAIDs)
    fprintf(fid,'     (1 %i (1) (2) (4) (18) (102) (-8) (-6 -7) (-6 -8))\n',200+zaididx);
end

fclose(fid);

cd(mycd);
