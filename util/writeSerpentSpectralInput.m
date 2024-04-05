function writeSerpentSpectralInput(p, g, serpentexecpath, ifname);
mycd = cd;
cd(serpentexecpath);

% TODO:
% 1- manage temperature input
% 2- therm lwtr lwj3.11t
% 3- set xsplot group def
% 4- group def stuff in general
%


fid = fopen(ifname,'w');

% Set overview settings
fprintf(fid, 'set title "%s, created %s"\n', ifname, datestr(now, 30));
fprintf(fid, 'set acelib "/home/fitze/LANL/SERPENT/sss_jeff311u.xsdata"\n');

%% Materials

% fuel
fprintf(fid, 'mat %s -%.6e\n', g.regionDef(1).cellDef(2).name, ...
    g.regionDef(1).cellDef(2).initDensity);
for iZAID = 1:length(g.regionDef(1).cellDef(2).initZAIDs)
    fprintf(fid, '%i.09c %.6e\n', ...
        g.regionDef(1).cellDef(2).initZAIDs(iZAID), ...
        g.regionDef(1).cellDef(2).initNumDensities(iZAID));
end
fprintf(fid, '\n');

% coolant
% Redefine g for defining water as a combination of 1001 and 8016.
if g.regionDef(1).cellDef(1).initZAIDs == 222
    g.regionDef(1).cellDef(1).initZAIDs = [1001 8016];
    N = g.regionDef(1).cellDef(1).initNumDensities;
    g.regionDef(1).cellDef(1).initNumDensities = [2*N N];
end

% Print ratios, etc.
if g.regionDef(1).cellDef(1).initZAIDs == 1001
    fprintf(fid, 'mat %s -%.6e moder lwtr 1001\n', ...
        g.regionDef(1).cellDef(1).name, ...
        g.regionDef(1).cellDef(1).initDensity);
else
    fprintf(fid, 'mat %s -%.6e\n', ...
        g.regionDef(1).cellDef(1).name, ...
        g.regionDef(1).cellDef(1).initDensity);
end

for iZAID = 1:length(g.regionDef(1).cellDef(1).initZAIDs)
    fprintf(fid, '%i.06c %.6e\n', ...
        g.regionDef(1).cellDef(1).initZAIDs(iZAID), ...
        g.regionDef(1).cellDef(1).initNumDensities(iZAID));
end

% Thermal scattering consideration.
if g.regionDef(1).cellDef(1).initZAIDs == 1001
    fprintf(fid, 'therm lwtr lwj3.11t\n');
end

%% Geometry

% Periodic boundary condition for infinite reactor.
fprintf(fid, 'set bc 3\n');

% Group constant generation (default, kinda)
fprintf(fid, 'set gcu 0\n');
fprintf(fid, 'set sym 8\n');
%fprintf(fid, 'set nfg 2 0.625\n');
fprintf(fid, 'set nfg 110 ');
fprintf(fid, '%.5e ', 10.^((-3.9:0.1:6.9) - 6));
fprintf(fid, '\n');

% Neutron source definition (stolen from PWR MOX example).
fprintf(fid, 'set pop 2000 100 20\n');

% Pin
fprintf(fid, 'pin 1\n');
fprintf(fid, '%s %.4e\n', g.regionDef(1).cellDef(2).name, p.uc.pinDiam / 2);
fprintf(fid, '%s\n', g.regionDef(1).cellDef(1).name);

% Surface for the pin.
fprintf(fid, 'surf 1000 sqc 0.0  0.0  %.4e\n', p.uc.pinPitch / 2);
fprintf(fid, '\n');
% Universe 0
fprintf(fid, 'cell 100 0 fill 1 -1000\n');
fprintf(fid, 'cell 101 0 outside 1000\n');

% Plot the reator geometry.
fprintf(fid, 'plot 3 500 500\n');

% Plot the results.
fprintf(fid, 'mesh 3 500 500\n');

% Assign cross section output data.
fprintf(fid, 'set xsplot 110 1e-10 10\n');

%% Detectors
fprintf(fid, 'ene egrid 3 110 1e-10 10\n');
fprintf(fid, 'det fuel dm %s de egrid\n', ...
    g.regionDef(1).cellDef(2).name);
fprintf(fid, 'det coolant dm %s de egrid\n', ...
    g.regionDef(1).cellDef(1).name);

% Detector materials
tempmodifier = {'.06c' '.09c'};
used_materials = [];
for iCell = 1:2
    for iZAID = 1:length(g.regionDef(1).cellDef(iCell).initZAIDs)
        if sum(g.regionDef(1).cellDef(iCell).initZAIDs(iZAID) == ...
            used_materials) == 0
            fprintf(fid, 'mat %i 1.0 %i%s 1.0\n', ...
                g.regionDef(1).cellDef(iCell).initZAIDs(iZAID), ...
                g.regionDef(1).cellDef(iCell).initZAIDs(iZAID), ...
                tempmodifier{iCell} ...
                );
            used_materials = [used_materials ...
                g.regionDef(1).cellDef(iCell).initZAIDs(iZAID)];
        end
    end
end


% fuel
cell_det_names = {'coolant', 'fuel'};
MTs = [1 2 18 51:91 102];
iDetector = 0;
% First attempt to obtain cell-level cross sections.
% for iCell = 1:2
%     for iMT = 1:length(MTs)
%         iDetector = iDetector + 1;
%         fprintf(fid, 'det c%sm%i dm %s de egrid dr %i %s dt 3 %s\n', ...
%             g.regionDef(1).cellDef(iCell).name, ...
%             MTs(iMT), ...
%             g.regionDef(1).cellDef(iCell).name, ...
%             MTs(iMT), ...
%             g.regionDef(1).cellDef(iCell).name, ...
%             cell_det_names{iCell});
%     end
% end
for iCell = 1:2
    for iZAID = 1:length(g.regionDef(1).cellDef(iCell).initZAIDs)
        for iMT = 1:length(MTs)
            if MTs(iMT) == 18 && (iCell == 1 || ...
                    g.regionDef(1).cellDef(iCell).initZAIDs(iZAID) == 8016)
                % Do not try to get fission for nonfissionable isotopes.
            elseif MTs(iMT) >= 51 && MTs(iMT) <= 91 && ...
                    g.regionDef(1).cellDef(iCell).initZAIDs(iZAID) == 1001
            elseif MTs(iMT) >= 58 && MTs(iMT) < 91 && ...
                    g.regionDef(1).cellDef(iCell).initZAIDs(iZAID) == 8016
            elseif MTs(iMT) >= 78 && MTs(iMT) < 91 && iCell == 1
            elseif MTs(iMT) >= 85 && MTs(iMT) < 91 && ...
                    g.regionDef(1).cellDef(iCell).initZAIDs(iZAID) == 92235
            else
                iDetector = iDetector + 1;
                fprintf(fid, 'det c%sz%im%i dm %s de egrid dr %i %i dt 3 %s\n', ...
                    g.regionDef(1).cellDef(iCell).name, ...
                    g.regionDef(1).cellDef(iCell).initZAIDs(iZAID), ...
                    MTs(iMT), ...
                    g.regionDef(1).cellDef(iCell).name, ...
                    MTs(iMT), ...
                    g.regionDef(1).cellDef(iCell).initZAIDs(iZAID), ...
                    cell_det_names{iCell});
            end
        end
    end
end


fclose(fid);

cd(mycd);
