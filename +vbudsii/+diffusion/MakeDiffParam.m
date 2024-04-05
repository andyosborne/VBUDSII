function DiffParam = MakeDiffParam(L, p, g, Region)
%MAKEDIFFPARAM Turn the structure Region into diffusion parameters ready
%for diffusion.
% Description:  From Region, scattering, absorption, transport, nu-fission,
% and fission are used to output the diffusion parameters scattering,
% xs-sink, nu-fission, diffusion coefficient, fission, and diffusion
% length. The input is ordered by increasing energy, but diffusion
% parameters are ordered as increasing lethargy (in group order). User has
% the option of displaying the six parameters and k_{inf} to the command
% window. Diffparam is a cell structure, where the number of cells is the
% number of regions, and each cell contains a structure whose fields are
% the diffusion parameters listed.
%    
% USE:  DiffParam = makeDiffParam(L,Region,g,p)
%
% NOTES: When producing the six parameters, a value for nu is required.
% This value as of now is hardcoded as 2.4, and not taken from Region.
%
%    DiffParam{regidx}.scatker
%                     .xssink
%                     .vfission
%                     .diffco
%                     .fission
%                     .difflength
%
% EXAMPLES: 
%
% MAJOR UPDATES:
%   version  date     NetID   description
%   1.0      20110519 cld72   cleaned up and formatted
%              
% FUTURE UPDATES:
%   1- grab nu instead of using a hard-coded value
%
% DEPENDENCIES: none
%

DiffParam = cell(1, g.nRegions);

for regidx = 1:g.nRegions

    % grab information from Region
    xsIn.scatker = Region(regidx).few(L.MT(2)).value;
    xsIn.absorb = Region(regidx).few(L.MT(6)).value;
    xsIn.transport = Region(regidx).few(L.MT(8)).value;
    xsIn.vfission = Region(regidx).few(L.MT(9)).value;
    xsIn.fission = Region(regidx).few(L.MT(18)).value;
    
    % (E',E) arrangement, switch scattering from increasing-energy to
    % increasing-group
    scatker = xsIn.scatker;
    groupscatker = zeros(p.nGroups);
    for i = 1:p.nGroups
        for j = 1:p.nGroups
            i2 = p.nGroups + 1 - i;
            j2 = p.nGroups + 1 - j;
            groupscatker(i,j) = scatker(j2,i2);
        end
    end
    
    % assemble xssink in increasing energy order first
    xssink = xsIn.absorb;
    
    for nrgidx = 2:p.nGroups % don't do this for the most thermal group
        
        xssink(nrgidx) = xssink(nrgidx) + sum(xsIn.scatker(1:nrgidx-1,nrgidx));
        
    end
    
    % CALCULATION OF DIFFUSION COEFFICIENT
    diffco = 1/3 ./ xsIn.transport(end:-1:1);
    
    % all are in INCREASING-GROUP order
    dp = struct('scatker',groupscatker,...
        'xssink',xssink(end:-1:1),...
        'vfission',xsIn.vfission(end:-1:1),...
        'diffco',diffco,...
        'fission',xsIn.fission(end:-1:1),...
        'difflength',sqrt(diffco./xsIn.absorb(end:-1:1)));
    
    DiffParam{regidx} = dp;
    
    %% show 6 parameters maybe
    
    if p.nGroups == 3 % 3-group, lets do our old stuff
        
        macroxs12 = xsIn.scatker(2,3);
        macroxs13 = xsIn.scatker(1,3);
        macroxs23 = xsIn.scatker(1,2);
        macroxsf1 = xsIn.fission(3);
        macroxsf2 = xsIn.fission(2);
        macroxsf3 = xsIn.fission(1);
        macroxsa1 = xsIn.absorb(3);
        macroxsa2 = xsIn.absorb(2);
        macroxsa3 = xsIn.absorb(1);
        diffco1 = 1/3/xsIn.transport(3);
        diffco2 = 1/3/xsIn.transport(2);
        diffco3 = 1/3/xsIn.transport(1);
                
        nu = 2.4;
        
        denom1 = macroxs12 + macroxs13 + macroxsa1;
        epsilon1 = nu*macroxsf1/denom1;
        epsilon2 = macroxs12/denom1;
        epsilon3 = macroxs13/denom1;
        
        denom2 = macroxs23 + macroxsa2;
        r        = nu*macroxsf2/denom2;
        pp       = macroxs23/denom2;
        
        etaf = nu*macroxsf3/macroxsa3;
        
        kinf1 = epsilon1 + epsilon2*r + etaf*(epsilon2*pp + epsilon3);
        
        %% display
        
        if ( p.printDiffParam )
            
            fprintf('\n3-Group Macroscopic Cross sections (1/cm) \n')
            fprintf('xstype     region-level\n')
            fprintf('----------------------------------------------\n')
            fprintf('12         %.6f\n',macroxs12)
            fprintf('13         %.6f\n',macroxs13)
            fprintf('23         %.6f\n',macroxs23)
            fprintf('1f         %.6f\n',macroxsf1)
            fprintf('2f         %.6f\n',macroxsf2)
            fprintf('3f         %.6f\n',macroxsf3)
            fprintf('1a         %.6f\n',macroxsa1)
            fprintf('2a         %.6f\n',macroxsa2)
            fprintf('3a         %.6f\n',macroxsa3)
            fprintf('1d         %.6f\n',diffco1)
            fprintf('2d         %.6f\n',diffco2)
            fprintf('3d         %.6f\n',diffco3)
            
            fprintf('\nCalculation for CADY method, first formulation\n')
            fprintf('----------------------------------------------\n')
            fprintf('epsilon1:  %.6f\n',epsilon1)
            fprintf('epsilon2:  %.6f\n',epsilon2)
            fprintf('epsilon3:  %.6f\n',epsilon3)
            fprintf('r:         %.6f\n',r)
            fprintf('p:         %.6f\n',pp)
            fprintf('eta f:     %.6f\n',etaf)
            fprintf('k inf:     %.6f\n\n',kinf1)
        end
        
        % store diffusion parameters to a file
        fid = fopen('makeDiffParamOutput.txt','a'); % append data to end of file
        fprintf(fid,['~~~~~~~ Run ' datestr(now,30) ' ~~~~~~~~~~~~~~\r']);
        fprintf(fid,'Calculation for CADY method, first formulation\r');
        fprintf(fid,'----------------------------------------------\r');
        fprintf(fid,'epsilon1:  %.6f\r',epsilon1);
        fprintf(fid,'epsilon2:  %.6f\r',epsilon2);
        fprintf(fid,'epsilon3:  %.6f\r',epsilon3);
        fprintf(fid,'r:         %.6f\r',r);
        fprintf(fid,'p:         %.6f\r',pp);
        fprintf(fid,'eta f:     %.6f\r',etaf);
        fprintf(fid,'k inf:     %.6f\r',kinf1);
        fclose(fid);
    end
end

end


