% this will test my Library (MakeLibrary) against Andy's OO code, when we do
% not include any Bondarenko considerations.


mycd = cd;
cd(fullfile('..','..',''));
vbudsiiDir = cd;
cd(mycd); % alternatively, mfilename

%% DEFINE PARAMETER STRUCTURE
fineGroupDef = 10.^(-4:.1:7);
fewGroupDef = [1e-4 1 100e3 1e7];
p = struct('nFineGroups',length(fineGroupDef)-1,... % leave blank.
           'nFewGroups',length(fewGroupDef)-1,...
           'fineGroupDef',fineGroupDef,...
           'fewGroupDef',fewGroupDef,...
           'temp',600,...
           'powerDensity',50,... % MW/m^3
           'nTimeSteps',1,...
           'vbudsiiDir',vbudsiiDir,...
           'makeLibrary',1,...
           'makeLibraryFromAndy',1,...
           'XSLibraryMAT',fullfile('..','..','data','XSLibrary.mat'),...
           'verbose',1);

addpath(fullfile(p.vbudsiiDir,'data',''));
addpath(fullfile(p.vbudsiiDir,'preprocessor',''));
addpath(fullfile(p.vbudsiiDir,'library',''));
addpath(fullfile(p.vbudsiiDir,'multicell',''));
addpath(fullfile(p.vbudsiiDir,'diffusion',''));
addpath(fullfile(p.vbudsiiDir,'bateman',''));
addpath(fullfile(p.vbudsiiDir,'postprocessor',''));
% I don't want to use this structure for the multicell source code.
addpath(fullfile(p.vbudsiiDir,'multicell','modules','CollisionProbability',''));
addpath(fullfile(p.vbudsiiDir,'multicell','modules','CrossSectionLibrary',''));
addpath(fullfile(p.vbudsiiDir,'multicell','modules','Material',''));
addpath(fullfile(p.vbudsiiDir,'multicell','modules','NjoyParameter',''));
addpath(fullfile(p.vbudsiiDir,'multicell','modules','NuclearPhysics',''));

%% RUN VBUDSII
%cd(fullfile('..','..',''));
%addpath(fullfile('..','..',''));
addpath(p.vbudsiiDir);

lib = CrossSectionLibrary.soleInstance();

L = MakeLibrary(p);
% for all zaids, all mts, temps, and bw
%loglog(p.fineGroupDef(1:p.nFineGroups),L.z(L.ZAID(92235)).m(L.MT(18)).xs(:,L.T(600),L.BW(1)),p.fineGroupDef(1:p.nFineGroups),lib.xsSpectrum('U235','fission',600));
%legend('L','lib')
%L.z(L.ZAID(92235)).m(L.MT(18)).xs(:,L.T(600),L.BW(10)) - lib.xsSpectrum('U235','fission',600);

z = {222, 'H2O';
     8016, 'O16';
     92235, 'U235';
     92238, 'U238'};

m = {4,'inelastic';
     18,'fission';
     102,'nGamma'};

close all;
for i = 1:4
    for j = [1 2 3]; %1:3
        for t = 600
            a = 0;
            for k = 1; %L.BWs
                a = a + 1;
                New = L.z(L.ZAID(z{i,1})).m(L.MT(m{j,1})).xs(:,L.T(t),L.BW(k));
                Old = lib.xsSpectrum(z{i,2}, m{j,2}, t);
                dif = New - Old;
                err = sum(abs(dif))/length(dif);
                fprintf('z: %d, m: %d, bw: %d ::   %.10f \n',z{i,1},m{j,1},k,err);
                if sum(abs(New)) == 0&& sum(abs(Old)) ~= 0
                    fprintf('z : %d, m: %d, bw: %d :: FLAG\n',z{i,1},m{j,1},k);
                end
                if a == 1 && sum((size(New) == size(Old))) == 2
%                subplot(length(L.BWs),1,a);
                    figure
                loglog(p.fineGroupDef(1:end-1)'*[1,1],[Old New]);
%                figure
 %               semilogx(p.fineGroupDef(1:end-1),dif);
                legend('old','new')
                end
             end
        end
    end
end

disp('SCAT');
for i = 1:4
    for t = 600
        for k = L.BWs
            New = L.z(L.ZAID(z{i,1})).m(L.MT(2)).xs(:,:,L.T(t),L.BW(k));
            Old = lib.elScatKernel(z{i,2}, t);
            dif = New - Old;
            err = sum(sum(abs(dif)))/length(dif);
            fprintf('z: %d, m: %d, bw: %d :: %.10f  \n',z{i,1},2,k,err);
            if sum(sum(abs(New))) == 0 && sum(sum(abs(Old))) ~= 0
                fprintf('z : %d, m: %d, bw: %d :: FLAG\n',z{i,1},2,k);
            end
        end
    end
end










%% CHECK RESULTS
%{
if 1...
    sum(r.Cell(2).spectralFlux ~= phi(:,2))
    disp('multicell3 test failed');
else
    disp('multicell3 test passed');
%    [Region(1).Cell(1).spectralFlux == phi(:,1), ...
%    Region(1).Cell(2).spectralFlux == phi(:,2)]
end
%}
