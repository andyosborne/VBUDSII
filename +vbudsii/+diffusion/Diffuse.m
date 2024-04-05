function [R,fluxOut, powerFrac] = Diffusion(DiffParamIn,DiffParamRefl)

global N R Delta

% 3-GROUP-ONLY CODE
numGroups = 3; % THIS HAS to be ascertained from inputs , like length(DiffParam{1}.vfission)

% see email on 101019 from recktenwald for the layout of this ThreeGroup
% code

% last edited (by and date): cld72@cornell.edu M101102
tic
plotCore = 1;

%  \                           \             \       \      \
%   |         campaign III      |  camp II    |   I   | refl | 
%  /                           /             /       /      /


        
  %% GEOMETRY INPUT
numCampaigns = 3;

% campaignMap, length

% WARNING: 3-REGION ONLY CODE (using region 4 as the moderator, with
% moderator defined as camp III's water/moderator)
% campaignMap = [ 3 3 3 3 3 3 3 3 3 3 3 2 2 2 2 2 2 2 2 1 1 1 1 1 1 4 4 4 4 4 4 4 4]; % USER INPUT OF GEOMETRY
% campaignMap = [ 3 3 2 2 1 1 4 4 4];
% campaignMap = [ 3*ones(1,30) 2*ones(1,30) 1*ones(1,30) 4*ones(1,100) ];
% initialize loop parameters / geometry
    
% gp = 9;

myMapcase = 2;
gp = 100;

switch myMapcase
    case 1
        campaignMap = [ 3 3 2 2 1 1 4 4 4];
    case 2
        campaignMap = [ 3*ones(1,gp) 2*ones(1,gp) 1*ones(1,gp) 4*ones(1,gp) ];
end

N = length(campaignMap); % number of grid points
k = 0:N;


% CoreRadius = 13.34; %14.53; % number of square fuel elements of size a

% 200 is radius of core in centimeters
%
% UNITS OF CM
% R = sqrt(CoreRadius^2 / CoreSpaces * k);

% reflR = sqrt(CoreRadius^2 / CoreSpaces * (N));


%% Process DiffParamIn

DiffParam = DiffParamIn;
DiffParam{4} = DiffParamRefl; % THIS IS VIOLATING STUFF, HARDCODED 4

% form of DiffParam
% DiffParam{campidx}.
%                  .scatker
%                  .xssink
%                  .vfission
%                  .diffco

Scatker = zeros(N,numGroups,numGroups);
Xssink = zeros(N,numGroups);
Vfission = zeros(N,numGroups);
Diffco = zeros(N,numGroups);
    
for spaceidx = 1:N
    
    campidx = campaignMap(spaceidx);
    
    Scatker(spaceidx,:,:) = DiffParam{campidx}.scatker;
    
    Xssink(spaceidx,:) = DiffParam{campidx}.xssink;
    Vfission(spaceidx,:) = DiffParam{campidx}.vfission;
    Diffco(spaceidx,:) = DiffParam{campidx}.diffco;
end

% actual geometry

CoreRadius = 200; % number of square fuel elements of size a
% CoreSpaces = N-sum(campaignMap==4); % side length of fuel elements (cm)
% 200 is radius of core in centimeters
%
% UNITS OF CM


CoreSpaces = N-sum(campaignMap==4); % side length of fuel elements (cm)
reflR = 1.3*CoreRadius;
campR = sqrt(CoreRadius^2 / 3 * (0:3));
%         R = sqrt(CoreRadius^2 / CoreSpaces * k);
campaign3R = linspace(campR(1),campR(2),sum(campaignMap==3)+1);
campaign2R = linspace(campR(2),campR(3),sum(campaignMap==2)+1);
campaign1R = linspace(campR(3),campR(4),sum(campaignMap==1)+1);
reflR = linspace(campR(4),reflR,sum(campaignMap==4)+1);
R = [campaign3R campaign2R(2:end) campaign1R(2:end) reflR(2:end)];
R = R';


% plot the core
close all
if plotCore
    h = figure(1);
    for i = 1:N
        hold on;
        thet = 0:.01:2*pi;
        x = R(i)*cos(thet);
        y = R(i)*sin(thet);
        plot(x,y,'Color',[0 .5 .25])
        axis equal
    end
end

% thickenss of each annular space,
Delta = diff(R);

%% BOUNDARY CONDITIONS
betaBC = [1 0];
extCurrent = [0 0];


%% ITERATION
% first fast source guess
avgdiffco = 0;
for i = 1:numCampaigns
    avgdiffco = avgdiffco + DiffParam{i}.diffco(1);
end
avgdiffco = avgdiffco / numCampaigns;
d = 2.13 * avgdiffco;
sourceFast = cos( R(2:end) / ( R (end) + d ) ); % space property
% load diffusesourceFast.mat
% prepare flux plot
% hflux = figure(2);
% gflux = axes;
% hold on
% set(hflux,'Visible','off')
% 
% % prepare sourceFast plot
% hsource = figure(3);
% gsource = axes;
% hold on
% set(hsource,'Visible','off')

% error
threshold = 1e-5;
isErrorful = 1;
maxIterations = 2000;
counter = 1;
ErrrVector = []; % used to plot error as a function of iteration
while isErrorful && counter <= maxIterations
    
    fluxOut = zeros(N+1,numGroups);
    
    %% Fast group diffusion
    groupidx = 1;
    
    % fast
    fluxOut(:,groupidx) = SingleGroup( sourceFast, Xssink(:,groupidx), Diffco(:,groupidx), betaBC, extCurrent );
    
    %% Resonance diffusion!
    groupidx = 2;
    
    % 3-GROUP-ONLY.
    % this calculation should really be done by the subfunction
    sourceRes = zeros(N,1);
    
    for spaceidx = 1:N
        nodalidx = spaceidx + 1;
        
        sourceRes(spaceidx) = sourceRes(spaceidx) + Scatker(spaceidx,groupidx-1,groupidx)...
                            * mean(fluxOut(nodalidx-1:nodalidx,1));
    end

    fluxOut(:,groupidx) = SingleGroup( sourceRes, Xssink(:,groupidx), Diffco(:,groupidx), betaBC, extCurrent);

    %% tHERMAL Diffusion!
    groupidx = 3;
    % 3group only!!!
    sourceTher = zeros(N,1);
    
    % EDIT BIG TIME , note that this can be done using matrix algebra
    % instead, across spaces
    
    for spaceidx = 1:N
        nodalidx = spaceidx + 1;
        
        for gidx = 1:groupidx-1
            
            sourceTher(spaceidx) = sourceTher(spaceidx) + Scatker(spaceidx,gidx,groupidx)...
                                * mean(fluxOut(nodalidx-1:nodalidx,gidx));
            
        end
    end
    
    fluxOut(:,groupidx) = SingleGroup( sourceTher, Xssink(:,groupidx), Diffco(:,groupidx), betaBC, extCurrent);

    %% Compute fast sources
    
    sourceFastNew = zeros(N,1);
    
    
    for spaceidx = 1:N
        for groupidx = 1:numGroups
            
            nodalidx = spaceidx + 1;
            
            sourceFastNew(spaceidx) = sourceFastNew(spaceidx) + Vfission(spaceidx,groupidx)...
                * mean(fluxOut(nodalidx-1:nodalidx,groupidx));
            
        end
    end
    
    Keff = norm(sourceFastNew)/norm(sourceFast);
    Keff2 = sourceFastNew(1)/sourceFast(1);
   
    sourceFastNew = sourceFastNew/sourceFastNew(1);

    errr = norm(sourceFastNew-sourceFast)/norm(sourceFastNew);
    
    ErrrVector = [ErrrVector errr];
    %% Plot fluxes
    
%     figure(hflux);
%     plot(R*[1 1 1],fluxOut)
%     set(hflux,'Visible','off')
%     
%     figure(hsource);
%     plot(mean([R(2:end) R(1:end-1)],2),sourceFast)
%     set(hsource,'Visible','off')
    
%     fprintf('Iteration Number %d: Error is %.4e \n',counter,errr)
    counter = counter+1;
    isErrorful = errr > threshold;
    sourceFast = sourceFastNew;
    
end

%% we have convergence

runtime = toc;
% flux
% figure(hflux);
% set(hflux,'Visible','on');
% xlabel('R (cm)')
% ylabel('\phi (# cm^{-2} s^{-1})')
% legend('fast','resonance','thermal')
% 
% figure(hsource);
% set(hsource,'Visible','on');
% xlabel('Radial distance to center of space (cm)')
% ylabel('Fast source')

%% FINAL FLUX PLOT
hfinalflux = figure(4);
gfinalflux = axes;
hold on;
plot(R*[1 1 1],fluxOut)
legend('fast','resonation','thermal') % EDIT resonation
xlabel({'R (cm)','',['Number of iterations: ' num2str(counter-1)],['Runtime: ' num2str(runtime) ' seconds']})
ylabel('\phi (# cm^{-2} s^{-1})')

bottoM = min(min(fluxOut));
toP = max(max(fluxOut));
for i = 1:N+1
    plot(R(i)*[1 1],[min(bottoM,0) 1.1*toP],':','Color',[.8 .8 .8])
end

for i = 2:4
    plot(campR(i)*[1 1],[min(bottoM,0) 1.1*toP],'k')%'LineWidth',1.5)
end

set(gfinalflux,'XLim',[0 R(end)],'YLim',[0 1.1*toP],'Box','on')


h5 = figure(5);
g5 = axes;
semilogy(ErrrVector)
xlabel('iteration (-)')
ylabel('Error (-)')

fluxOut;
fprintf('keff is: %.4f\n',Keff);
fprintf('keff2 is: %.4f\n',Keff2);

%% Output flux

% fid = fopen('DiffusionOutput.txt','a');
% % 3-GROUP ONLY
% fprintf(fid,['\r~~~~~~~ Run ' datestr(now,30) ' ~~~~~~~~~~~~~~\r']);
% for i = 1:N+1
%     fprintf(fid,'%7.5f   %7.5f   %7.5f\r',fluxOut(i,:));
% end
% 
% fclose(fid);

%% Calculate power fractions

powerSpace = zeros(N,1);

for spaceidx = 1:N % can 
    
    nodalidx = spaceidx + 1;
    
    for groupidx = 1:numGroups % can do away with the energy loop using matrix algebra
        
        powerSpace(spaceidx) = powerSpace(spaceidx) + Vfission(spaceidx,groupidx)...
                              * mean(fluxOut(nodalidx-1:nodalidx,groupidx));
        
    end
    
end

powerTotal = sum(powerSpace);

powerFrac = zeros(numCampaigns,1);

campaign3idxs = find(campaignMap==3);
campaign2idxs = find(campaignMap==2);
campaign1idxs = find(campaignMap==1);
powerFrac(3) = sum( powerSpace( campaign3idxs ) ) / powerTotal; % campaign III
powerFrac(2) = sum( powerSpace( campaign2idxs ) ) / powerTotal; % campaign II
powerFrac(1) = sum( powerSpace( campaign1idxs ) ) / powerTotal; % campaign I

disp('Power Fractions!')
for i = numCampaigns:-1:1
    fprintf('Campaign %d: %.5f \n',i,powerFrac(i))
end

fprintf('Number of iterations: %d\nRuntime: %.2f seconds\n',counter-1,runtime);

% save the fast source

save diffusesourceFast.mat sourceFast
end






















function [fluxOut] = SingleGroup(sourceIn,xssink,diffco,betaBC,extCurrent)

global N R Delta 

% from cady packet 1 and from everything else too (typed writeup of Cady,
% Recktenwald's PDF writing on 101019...
%% GRIDSPACE LOOP!

% e, d, k defined for space
e = zeros(N,1);
d = zeros(N,1);
t = zeros(N,1);

% a, c, b, sigma defined for grid points
a = zeros(N,1);
c = zeros(N,1);

for spaceidx = 1:N % loop through spaces, get a, b, c, sigma
    
    % defined such that nodalidx = 1 when k = 0 (centerline) and
    % nodalidx = N+1 when k = N (boundary)
    nodalidx = spaceidx + 1;
    
    e(spaceidx) = 0.5 * xssink(spaceidx) * Delta(spaceidx);
    d(spaceidx) = 2.0 * diffco(spaceidx) / Delta(spaceidx);
    % for cylinders. t = 1/2 for slabs, r(nodalidx-1)^2 / (r(nodalidx)^2 +
    % r(nodalidx-1)^2) for spheres
    
    p = 1;
    
    t(spaceidx) = R(nodalidx-1)^p / ( R(nodalidx)^p + R(nodalidx-1)^p ); % R is defined nodally
    
    a(spaceidx) = (1 - t(spaceidx)) * ( e(spaceidx) - d(spaceidx) );
    c(spaceidx) = t(spaceidx)*( e(spaceidx) - d(spaceidx) );
    
    % reconsider the use of such long indices
    
end

% deal with b and sigma now
b = zeros(N+1,1);
sigma = zeros(N+1,1);

% centerline
nodalidx = 1; % k = 0
spaceidx = 1; % k = 1

b(nodalidx) = t(spaceidx) * e(spaceidx) + ( 1 - t(spaceidx) ) * d(spaceidx)...
             + 0.5 * (1 - betaBC(1)) / (1 + betaBC(1));

sigma(nodalidx) = 0.5 * sourceIn(spaceidx)*Delta(spaceidx) + 2 * extCurrent(1) / (1 + betaBC(1));

% edge of reactor (end of reflector...)
nodalidx = N+1;
spaceidx = N;

b(nodalidx) = (1 - t(spaceidx))*e(spaceidx) +  (t(spaceidx) * d(spaceidx))...
             + 0.5 * (1 - betaBC(2)) / (1 + betaBC(2));

sigma(nodalidx) = 0.5 * sourceIn(spaceidx)*Delta(spaceidx) + 2 * extCurrent(2) / (1 + betaBC(2));

for nodalidx = 2:N % nodalidx = 1 (k = 0) and nodalidx = N+1 (k = N) taken care of
    
    spaceidx = nodalidx-1;
    
    % b(k) = (1-t(k))e(k) + t(k)*d(k) + t(k+1)*e(k+1) + (1-t(k+1)*d(k+1)
    b(nodalidx) = (1 - t(spaceidx))*e(spaceidx) +  (t(spaceidx) * d(spaceidx))...
        + (t(spaceidx+1) * e(spaceidx+1)) + (1 - t(spaceidx+1))*d(spaceidx+1);
    
    sigma(nodalidx) = 0.5 * (sourceIn(spaceidx)*Delta(spaceidx) + sourceIn(spaceidx+1)*Delta(spaceidx+1));
    % again, source is defined on grid point nodes, and Delta is defined on
    % spaces
        
end

%% Solve the matrix!

    fluxOut = thomas(b, a, c, sigma);

end


%% usefull

%         set(h,'NextPlot','replacechildren')

% This code performs the following operations:
%
% calculates diffusion parameters
% plugs stuff into a diffusion equation ??

% INPUTS
% DiffParam cell array, each cell holding like 12 values for each region
% MUST BE FLEXIBLE OUTPUT FOR # of groups
% homogenized; all it sees are the 3-group numbers for a single homogenized
% medium
% output:
%
% last edited: cld72 A101022


% explain nodal indices (since there's a 0th node, we must shift over 1
% not for the space indices
% DiffParam{campidx}.
%                  .absorb
%                  .vfission
%                  .scatker
%                  .diffco

% Nodal Parameters:
% Flux
% 
% Space Parameters:
% Delta
% Source



% [b0 a1                  ]
% [c1 b1 a2               ]
% [   c2 b2 a3            ]
% [      c3 b3            ]
% [             .         ]
% [               .       ]
% [                 .     ]
% [                     aN]
% [                  cN bN]

% phsyical meaning of a,b,c
      

% this initial flux guess does not belong here. where do we use our
% input flux?

%% trash

%     if mod(counter,10)==1
%         figure(h3);
%         plot(R*[1 1 1],fluxOut)
%         %         set(h3,'Visible','off')
%         xlabel(['iteration number ' num2str(counter) ', error is ' num2str(errr)])
%         legend('fast','res','thermal')
%         pause(.2)
%         
%     end
% 	% plot fast source
%     figure(h3);
%     plot(mean([R(2:end) R(1:end-1)],2),sourceFast)
%     xlabel('radial distance to center of space (cm)')
%     ylabel('Fast source')
%     set(h3,'Visible','off')


    
%     e
%     t
% %     d%     
%     i = [2:N+1,   1:N,   1:N+1];
%     j = [1:N,   2:N+1,   1:N+1];
%     s = [a; c; b];
%     
%     
%     s = [c; a; b];
%     
%     K = sparse(i,j,s);
%     
%     K = full(K);
%     K;
%     det(K);
%     fluxOut;
%     sigma;
