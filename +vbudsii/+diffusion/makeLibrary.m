function Library = makeLibrary(p)
%MAKELIBRARY Reformat ENDF library into appropriate VBUDS format.
% Description:  This file will take a material matrix and return the atomic
% percentages, weight percentages, atoms/b-cm, and density.
%
% USE:  Library = makeLibrary(p)
%
% NOTES: Input p is not actually used in this function.
%
% EXAMPLES:
%
% MAJOR UPDATES:
%   version  date     NetID   description
%   1.0      20110327 cld72   created
%   2.0      20110501 cld72   implementation after forgetting everything
%   3.0      20110519 cld72   submitting some final work this week
%
% FUTURE UPDATES:
%   1- use more appropriate BW values.
%
% DEPENDENCIES:
%   reducedlib.mat
%   makeTransport2
%

load('reducedlib.mat');

% Library = struct('ZAID',[ 1001,
%                           8016,
%                          10018,
%                          92235,
%                          92238],...
%                  'MT',sparse([2 5 18 102 251 452]),... %nufission?
%                  'Temp',[300,600,900,1200,1500],... % kelvin
%                  'BW',[0 1 2 3 5 10],...
%                  'Bin',[1 3 110 1100],...
%                  'xs',[]);

% list the ZAIDs, MTs, Ts, BWs, and Bins that can be used
ZAIDs = [222 1001 8016 92235 92238];
MTs = [2 5 7 16 18 102 251 452]; % 16 is n2n, 7 is nu fission (mcnp table3.5)
Ts = [300 600 900 1200 1500];
BWs = [0 1 2 3 5 10];
Bins = [1 3 110 1100];
% bb=sparse([0]); want to incude -1 for BW

% initialize the library
L = struct('ZAIDs',ZAIDs,...
           'MTs',MTs,...
           'Ts',Ts,...
           'BWs',BWs,...
           'Bins',Bins,...
           'ZAID',[],...
           'MT',[],... %nufission?
           'T',[],... % kelvin
           'BW',[],...
           'Bin',[],...
           'xss',struct('NBins',struct('xs',[])));

% create mapping of values to an array index.
L.ZAID(ZAIDs) = 1:length(ZAIDs);
L.MT(MTs) = 1:length(MTs);
L.T(Ts) = 1:length(Ts);
% L.BW(BWs) = 1:length(BWs);
L.BW = containers.Map({0 1 2 3 5 10},{1 2 3 4 5 6});
L.Bin(Bins) = 1:length(Bins);

% make mapping fields sparse to save space
L.ZAID = sparse(L.ZAID);
L.MT = sparse(L.MT);
L.T = sparse(L.T);
L.Bin = sparse(L.Bin);

nTs = length(Ts);
nBWs = length(BWs);
nEThree = 3;
nECoarse = 110;
nEFine = 1100;

% initialize the size of the library.
for i = ZAIDs
    for m = Bins
        if m == 110
            for j = MTs
                if j == 2
                    L.xss(L.ZAID(i),L.MT(j)).NBins(L.Bin(m)).xs = zeros(nTs,nBWs,m,m);
                else
                    L.xss(L.ZAID(i),L.MT(j)).NBins(L.Bin(m)).xs = zeros(nTs,nBWs,m);
                end
            end
        end
    end
end
%% ZAID   222: WATER

%   222 scatter kernel
for i = Ts
    L.xss(L.ZAID(222),L.MT(2)).NBins(L.Bin(110)).xs(L.T(i),L.BW(0),:,:)   = H2O_0_2_600_0;
    L.xss(L.ZAID(222),L.MT(2)).NBins(L.Bin(110)).xs(L.T(i),L.BW(1),:,:)   = H2O_0_2_600_1;
    L.xss(L.ZAID(222),L.MT(2)).NBins(L.Bin(110)).xs(L.T(i),L.BW(2),:,:)   = H2O_0_2_600_2;
    L.xss(L.ZAID(222),L.MT(2)).NBins(L.Bin(110)).xs(L.T(i),L.BW(3),:,:)   = H2O_0_2_600_3;
    L.xss(L.ZAID(222),L.MT(2)).NBins(L.Bin(110)).xs(L.T(i),L.BW(5),:,:)   = H2O_0_2_600_5;
    L.xss(L.ZAID(222),L.MT(2)).NBins(L.Bin(110)).xs(L.T(i),L.BW(10),:,:)   = H2O_0_2_600_10;
end

%   222 transport

%   222 fission
% for k = Bins
%     L.xss(L.ZAID(1001),L.MT(18)).NBins(L.Bin(k)).xs = zeros(nTs,nBWs,k);
% end

%   222 gamma
for i = Ts
    L.xss(L.ZAID(222),L.MT(102)).NBins(L.Bin(110)).xs(L.T(i),L.BW(0),:,:)   = H2O_3_102_600_0;
    L.xss(L.ZAID(222),L.MT(102)).NBins(L.Bin(110)).xs(L.T(i),L.BW(1),:,:)   = H2O_3_102_600_1;
    L.xss(L.ZAID(222),L.MT(102)).NBins(L.Bin(110)).xs(L.T(i),L.BW(2),:,:)   = H2O_3_102_600_2;
    L.xss(L.ZAID(222),L.MT(102)).NBins(L.Bin(110)).xs(L.T(i),L.BW(3),:,:)   = H2O_3_102_600_3;
    L.xss(L.ZAID(222),L.MT(102)).NBins(L.Bin(110)).xs(L.T(i),L.BW(5),:,:)   = H2O_3_102_600_5;
    L.xss(L.ZAID(222),L.MT(102)).NBins(L.Bin(110)).xs(L.T(i),L.BW(10),:,:)   = H2O_3_102_600_10;
end

%   222 mubar
for i = Ts
    L.xss(L.ZAID(222),L.MT(251)).NBins(L.Bin(110)).xs(L.T(i),L.BW(0),:,:)   = H2O_3_251_600_0;
    L.xss(L.ZAID(222),L.MT(251)).NBins(L.Bin(110)).xs(L.T(i),L.BW(1),:,:)   = H2O_3_251_600_1;
    L.xss(L.ZAID(222),L.MT(251)).NBins(L.Bin(110)).xs(L.T(i),L.BW(2),:,:)   = H2O_3_251_600_2;
    L.xss(L.ZAID(222),L.MT(251)).NBins(L.Bin(110)).xs(L.T(i),L.BW(3),:,:)   = H2O_3_251_600_3;
    L.xss(L.ZAID(222),L.MT(251)).NBins(L.Bin(110)).xs(L.T(i),L.BW(5),:,:)   = H2O_3_251_600_5;
    L.xss(L.ZAID(222),L.MT(251)).NBins(L.Bin(110)).xs(L.T(i),L.BW(10),:,:)   = H2O_3_251_600_10;
end

%  8016 nu. only relevant for fissile isotopes

%% ZAID  1001: where is H1w dealt with?
%  1001 scatter kernel
for i = BWs
    L.xss(L.ZAID(1001),L.MT(2)).NBins(L.Bin(110)).xs(L.T(300),L.BW(i),:,:) = H1_0_2_300;
    L.xss(L.ZAID(1001),L.MT(2)).NBins(L.Bin(110)).xs(L.T(300),L.BW(i),1:45,1:45) = H1w_0_2_296(1:45,1:45);
    L.xss(L.ZAID(1001),L.MT(2)).NBins(L.Bin(110)).xs(L.T(600),L.BW(i),:,:) = H1_0_2_600;
    L.xss(L.ZAID(1001),L.MT(2)).NBins(L.Bin(110)).xs(L.T(600),L.BW(i),1:45,1:45) = H1w_0_2_600(1:45,1:45);
    L.xss(L.ZAID(1001),L.MT(2)).NBins(L.Bin(110)).xs(L.T(900),L.BW(i),:,:) = H1_0_2_900;
    L.xss(L.ZAID(1001),L.MT(2)).NBins(L.Bin(110)).xs(L.T(900),L.BW(i),1:45,1:45) = H1w_0_2_1000(1:45,1:45);
    L.xss(L.ZAID(1001),L.MT(2)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(i),:,:) = H1_0_2_1200;
    L.xss(L.ZAID(1001),L.MT(2)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(i),:,:) = H1_0_2_1500;
end

%  1001 transport

%  1001 fission

%  1001 gamma
L.xss(L.ZAID(1001),L.MT(102)).NBins(L.Bin(110)).xs(L.T(300),L.BW(0),:)   = H1_3_102_300_0;
L.xss(L.ZAID(1001),L.MT(102)).NBins(L.Bin(110)).xs(L.T(300),L.BW(1),:)   = H1_3_102_300_1;
L.xss(L.ZAID(1001),L.MT(102)).NBins(L.Bin(110)).xs(L.T(300),L.BW(2),:)   = H1_3_102_300_2;
L.xss(L.ZAID(1001),L.MT(102)).NBins(L.Bin(110)).xs(L.T(300),L.BW(3),:)   = H1_3_102_300_3;
L.xss(L.ZAID(1001),L.MT(102)).NBins(L.Bin(110)).xs(L.T(300),L.BW(5),:)   = H1_3_102_300_5;
L.xss(L.ZAID(1001),L.MT(102)).NBins(L.Bin(110)).xs(L.T(300),L.BW(10),:)   = H1_3_102_300_10;

L.xss(L.ZAID(1001),L.MT(102)).NBins(L.Bin(110)).xs(L.T(600),L.BW(0),:)   = H1_3_102_600_0;
L.xss(L.ZAID(1001),L.MT(102)).NBins(L.Bin(110)).xs(L.T(600),L.BW(1),:)   = H1_3_102_600_1;
L.xss(L.ZAID(1001),L.MT(102)).NBins(L.Bin(110)).xs(L.T(600),L.BW(2),:)   = H1_3_102_600_2;
L.xss(L.ZAID(1001),L.MT(102)).NBins(L.Bin(110)).xs(L.T(600),L.BW(3),:)   = H1_3_102_600_3;
L.xss(L.ZAID(1001),L.MT(102)).NBins(L.Bin(110)).xs(L.T(600),L.BW(5),:)   = H1_3_102_600_5;
L.xss(L.ZAID(1001),L.MT(102)).NBins(L.Bin(110)).xs(L.T(600),L.BW(10),:)   = H1_3_102_600_10;

L.xss(L.ZAID(1001),L.MT(102)).NBins(L.Bin(110)).xs(L.T(900),L.BW(0),:)   = H1_3_102_900_0;
L.xss(L.ZAID(1001),L.MT(102)).NBins(L.Bin(110)).xs(L.T(900),L.BW(1),:)   = H1_3_102_900_1;
L.xss(L.ZAID(1001),L.MT(102)).NBins(L.Bin(110)).xs(L.T(900),L.BW(2),:)   = H1_3_102_900_2;
L.xss(L.ZAID(1001),L.MT(102)).NBins(L.Bin(110)).xs(L.T(900),L.BW(3),:)   = H1_3_102_900_3;
L.xss(L.ZAID(1001),L.MT(102)).NBins(L.Bin(110)).xs(L.T(900),L.BW(5),:)   = H1_3_102_900_5;
L.xss(L.ZAID(1001),L.MT(102)).NBins(L.Bin(110)).xs(L.T(900),L.BW(10),:)   = H1_3_102_900_10;

L.xss(L.ZAID(1001),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(0),:)   = H1_3_102_1200_0;
L.xss(L.ZAID(1001),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(1),:)   = H1_3_102_1200_1;
L.xss(L.ZAID(1001),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(2),:)   = H1_3_102_1200_2;
L.xss(L.ZAID(1001),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(3),:)   = H1_3_102_1200_3;
L.xss(L.ZAID(1001),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(5),:)   = H1_3_102_1200_5;
L.xss(L.ZAID(1001),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(10),:)   = H1_3_102_1200_10;

L.xss(L.ZAID(1001),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(0),:)   = H1_3_102_1500_0;
L.xss(L.ZAID(1001),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(1),:)   = H1_3_102_1500_1;
L.xss(L.ZAID(1001),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(2),:)   = H1_3_102_1500_2;
L.xss(L.ZAID(1001),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(3),:)   = H1_3_102_1500_3;
L.xss(L.ZAID(1001),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(5),:)   = H1_3_102_1500_5;
L.xss(L.ZAID(1001),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(10),:)   = H1_3_102_1500_10;

%  1001 mubar
for i = BWs(1:end-1)
    L.xss(L.ZAID(1001),L.MT(251)).NBins(L.Bin(110)).xs(L.T(300),L.BW(i),:)   = H1_3_251_300;
    L.xss(L.ZAID(1001),L.MT(251)).NBins(L.Bin(110)).xs(L.T(600),L.BW(i),:)   = H1_3_251_600;
    L.xss(L.ZAID(1001),L.MT(251)).NBins(L.Bin(110)).xs(L.T(900),L.BW(i),:)   = H1_3_251_900;
    L.xss(L.ZAID(1001),L.MT(251)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(i),:)   = H1_3_251_1200;
    L.xss(L.ZAID(1001),L.MT(251)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(i),:)   = H1_3_251_1500;
end
L.xss(L.ZAID(1001),L.MT(251)).NBins(L.Bin(110)).xs(L.T(300),L.BW(10),:)   = H1_3_251_300_10;
L.xss(L.ZAID(1001),L.MT(251)).NBins(L.Bin(110)).xs(L.T(600),L.BW(10),:)   = H1_3_251_600_10;
L.xss(L.ZAID(1001),L.MT(251)).NBins(L.Bin(110)).xs(L.T(900),L.BW(10),:)   = H1_3_251_900_10;
L.xss(L.ZAID(1001),L.MT(251)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(10),:)   = H1_3_251_1200_10;
L.xss(L.ZAID(1001),L.MT(251)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(10),:)   = H1_3_251_1500_10;

%  1001 nu. only relevant for fissile isotopes.

%% ZAID  8016
%  8016 scatter kernel
for i = BWs
    L.xss(L.ZAID(8016),L.MT(2)).NBins(L.Bin(110)).xs(L.T(300),L.BW(i),:,:) = O16_0_2_300;
    L.xss(L.ZAID(8016),L.MT(2)).NBins(L.Bin(110)).xs(L.T(600),L.BW(i),:,:) = O16_0_2_600;
    L.xss(L.ZAID(8016),L.MT(2)).NBins(L.Bin(110)).xs(L.T(900),L.BW(i),:,:) = O16_0_2_900;
    L.xss(L.ZAID(8016),L.MT(2)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(i),:,:) = O16_0_2_1200;
    L.xss(L.ZAID(8016),L.MT(2)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(i),:,:) = O16_0_2_1500;
end

%  8016 transport

%  8016 fission

%  8016 gamma
L.xss(L.ZAID(8016),L.MT(102)).NBins(L.Bin(110)).xs(L.T(300),L.BW(0),:)   = O16_3_102_300_0;
L.xss(L.ZAID(8016),L.MT(102)).NBins(L.Bin(110)).xs(L.T(300),L.BW(1),:)   = O16_3_102_300_1;
L.xss(L.ZAID(8016),L.MT(102)).NBins(L.Bin(110)).xs(L.T(300),L.BW(2),:)   = O16_3_102_300_2;
L.xss(L.ZAID(8016),L.MT(102)).NBins(L.Bin(110)).xs(L.T(300),L.BW(3),:)   = O16_3_102_300_3;
L.xss(L.ZAID(8016),L.MT(102)).NBins(L.Bin(110)).xs(L.T(300),L.BW(5),:)   = O16_3_102_300_5;
L.xss(L.ZAID(8016),L.MT(102)).NBins(L.Bin(110)).xs(L.T(300),L.BW(10),:)   = O16_3_102_300_10;

L.xss(L.ZAID(8016),L.MT(102)).NBins(L.Bin(110)).xs(L.T(600),L.BW(0),:)   = O16_3_102_600_0;
L.xss(L.ZAID(8016),L.MT(102)).NBins(L.Bin(110)).xs(L.T(600),L.BW(1),:)   = O16_3_102_600_1;
L.xss(L.ZAID(8016),L.MT(102)).NBins(L.Bin(110)).xs(L.T(600),L.BW(2),:)   = O16_3_102_600_2;
L.xss(L.ZAID(8016),L.MT(102)).NBins(L.Bin(110)).xs(L.T(600),L.BW(3),:)   = O16_3_102_600_3;
L.xss(L.ZAID(8016),L.MT(102)).NBins(L.Bin(110)).xs(L.T(600),L.BW(5),:)   = O16_3_102_600_5;
L.xss(L.ZAID(8016),L.MT(102)).NBins(L.Bin(110)).xs(L.T(600),L.BW(10),:)   = O16_3_102_600_10;

L.xss(L.ZAID(8016),L.MT(102)).NBins(L.Bin(110)).xs(L.T(900),L.BW(0),:)   = O16_3_102_900_0;
L.xss(L.ZAID(8016),L.MT(102)).NBins(L.Bin(110)).xs(L.T(900),L.BW(1),:)   = O16_3_102_900_1;
L.xss(L.ZAID(8016),L.MT(102)).NBins(L.Bin(110)).xs(L.T(900),L.BW(2),:)   = O16_3_102_900_2;
L.xss(L.ZAID(8016),L.MT(102)).NBins(L.Bin(110)).xs(L.T(900),L.BW(3),:)   = O16_3_102_900_3;
L.xss(L.ZAID(8016),L.MT(102)).NBins(L.Bin(110)).xs(L.T(900),L.BW(5),:)   = O16_3_102_900_5;
L.xss(L.ZAID(8016),L.MT(102)).NBins(L.Bin(110)).xs(L.T(900),L.BW(10),:)   = O16_3_102_900_10;

L.xss(L.ZAID(8016),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(0),:)   = O16_3_102_1200_0;
L.xss(L.ZAID(8016),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(1),:)   = O16_3_102_1200_1;
L.xss(L.ZAID(8016),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(2),:)   = O16_3_102_1200_2;
L.xss(L.ZAID(8016),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(3),:)   = O16_3_102_1200_3;
L.xss(L.ZAID(8016),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(5),:)   = O16_3_102_1200_5;
L.xss(L.ZAID(8016),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(10),:)   = O16_3_102_1200_10;

L.xss(L.ZAID(8016),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(0),:)   = O16_3_102_1500_0;
L.xss(L.ZAID(8016),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(1),:)   = O16_3_102_1500_1;
L.xss(L.ZAID(8016),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(2),:)   = O16_3_102_1500_2;
L.xss(L.ZAID(8016),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(3),:)   = O16_3_102_1500_3;
L.xss(L.ZAID(8016),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(5),:)   = O16_3_102_1500_5;
L.xss(L.ZAID(8016),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(10),:)   = O16_3_102_1500_10;

%  8016 mubar
for i = BWs(1:end-1)
    L.xss(L.ZAID(8016),L.MT(251)).NBins(L.Bin(110)).xs(L.T(300),L.BW(i),:)   = O16_3_251_300;
    L.xss(L.ZAID(8016),L.MT(251)).NBins(L.Bin(110)).xs(L.T(600),L.BW(i),:)   = O16_3_251_600;
    L.xss(L.ZAID(8016),L.MT(251)).NBins(L.Bin(110)).xs(L.T(900),L.BW(i),:)   = O16_3_251_900;
    L.xss(L.ZAID(8016),L.MT(251)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(i),:)   = O16_3_251_1200;
    L.xss(L.ZAID(8016),L.MT(251)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(i),:)   = O16_3_251_1500;

end
L.xss(L.ZAID(8016),L.MT(251)).NBins(L.Bin(110)).xs(L.T(300),L.BW(10),:)   = O16_3_251_300_10;
L.xss(L.ZAID(8016),L.MT(251)).NBins(L.Bin(110)).xs(L.T(600),L.BW(10),:)   = O16_3_251_600_10;
L.xss(L.ZAID(8016),L.MT(251)).NBins(L.Bin(110)).xs(L.T(900),L.BW(10),:)   = O16_3_251_900_10;
L.xss(L.ZAID(8016),L.MT(251)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(10),:)   = O16_3_251_1200_10;
L.xss(L.ZAID(8016),L.MT(251)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(10),:)   = O16_3_251_1500_10;

%  8016 nu. only relevant for fissile isotopes

%% ZAID 92235

% 92235 scatter kernel
for i = BWs
    L.xss(L.ZAID(92235),L.MT(2)).NBins(L.Bin(110)).xs(L.T(300),L.BW(i),:,:) = U235_0_2_300;
    L.xss(L.ZAID(92235),L.MT(2)).NBins(L.Bin(110)).xs(L.T(600),L.BW(i),:,:) = U235_0_2_600;
    L.xss(L.ZAID(92235),L.MT(2)).NBins(L.Bin(110)).xs(L.T(900),L.BW(i),:,:) = U235_0_2_900;
    L.xss(L.ZAID(92235),L.MT(2)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(i),:,:) = U235_0_2_1200;
    L.xss(L.ZAID(92235),L.MT(2)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(i),:,:) = U235_0_2_1500;
end

% 92235 transport
% for i = BWs
%     L.xss(L.ZAID(92235),L.MT(5)).NBins(L.Bin(110)).xs(L.T(300),L.BW(i),:)   = U235_3_5_300;
%     L.xss(L.ZAID(92235),L.MT(5)).NBins(L.Bin(110)).xs(L.T(600),L.BW(i),:)   = U235_3_5_600;
%     L.xss(L.ZAID(92235),L.MT(5)).NBins(L.Bin(110)).xs(L.T(900),L.BW(i),:)   = U235_3_5_900;
%     L.xss(L.ZAID(92235),L.MT(5)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(i),:)   = U235_3_5_1200;
%     L.xss(L.ZAID(92235),L.MT(5)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(i),:)   = U235_3_5_1500;
% end

%  92235 fission
L.xss(L.ZAID(92235),L.MT(18)).NBins(L.Bin(110)).xs(L.T(300),L.BW(0),:)   = U235_3_18_300_0;
L.xss(L.ZAID(92235),L.MT(18)).NBins(L.Bin(110)).xs(L.T(300),L.BW(1),:)   = U235_3_18_300_1;
L.xss(L.ZAID(92235),L.MT(18)).NBins(L.Bin(110)).xs(L.T(300),L.BW(2),:)   = U235_3_18_300_2;
L.xss(L.ZAID(92235),L.MT(18)).NBins(L.Bin(110)).xs(L.T(300),L.BW(3),:)   = U235_3_18_300_3;
L.xss(L.ZAID(92235),L.MT(18)).NBins(L.Bin(110)).xs(L.T(300),L.BW(5),:)   = U235_3_18_300_5;
L.xss(L.ZAID(92235),L.MT(18)).NBins(L.Bin(110)).xs(L.T(300),L.BW(10),:)   = U235_3_18_300_10;

L.xss(L.ZAID(92235),L.MT(18)).NBins(L.Bin(110)).xs(L.T(600),L.BW(0),:)   = U235_3_18_600_0;
L.xss(L.ZAID(92235),L.MT(18)).NBins(L.Bin(110)).xs(L.T(600),L.BW(1),:)   = U235_3_18_600_1;
L.xss(L.ZAID(92235),L.MT(18)).NBins(L.Bin(110)).xs(L.T(600),L.BW(2),:)   = U235_3_18_600_2;
L.xss(L.ZAID(92235),L.MT(18)).NBins(L.Bin(110)).xs(L.T(600),L.BW(3),:)   = U235_3_18_600_3;
L.xss(L.ZAID(92235),L.MT(18)).NBins(L.Bin(110)).xs(L.T(600),L.BW(5),:)   = U235_3_18_600_5;
L.xss(L.ZAID(92235),L.MT(18)).NBins(L.Bin(110)).xs(L.T(600),L.BW(10),:)   = U235_3_18_600_10;

L.xss(L.ZAID(92235),L.MT(18)).NBins(L.Bin(110)).xs(L.T(900),L.BW(0),:)   = U235_3_18_900_0;
L.xss(L.ZAID(92235),L.MT(18)).NBins(L.Bin(110)).xs(L.T(900),L.BW(1),:)   = U235_3_18_900_1;
L.xss(L.ZAID(92235),L.MT(18)).NBins(L.Bin(110)).xs(L.T(900),L.BW(2),:)   = U235_3_18_900_2;
L.xss(L.ZAID(92235),L.MT(18)).NBins(L.Bin(110)).xs(L.T(900),L.BW(3),:)   = U235_3_18_900_3;
L.xss(L.ZAID(92235),L.MT(18)).NBins(L.Bin(110)).xs(L.T(900),L.BW(5),:)   = U235_3_18_900_5;
L.xss(L.ZAID(92235),L.MT(18)).NBins(L.Bin(110)).xs(L.T(900),L.BW(10),:)   = U235_3_18_900_10;

L.xss(L.ZAID(92235),L.MT(18)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(0),:)   = U235_3_18_1200_0;
L.xss(L.ZAID(92235),L.MT(18)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(1),:)   = U235_3_18_1200_1;
L.xss(L.ZAID(92235),L.MT(18)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(2),:)   = U235_3_18_1200_2;
L.xss(L.ZAID(92235),L.MT(18)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(3),:)   = U235_3_18_1200_3;
L.xss(L.ZAID(92235),L.MT(18)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(5),:)   = U235_3_18_1200_5;
L.xss(L.ZAID(92235),L.MT(18)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(10),:)   = U235_3_18_1200_10;

L.xss(L.ZAID(92235),L.MT(18)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(0),:)   = U235_3_18_1500_0;
L.xss(L.ZAID(92235),L.MT(18)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(1),:)   = U235_3_18_1500_1;
L.xss(L.ZAID(92235),L.MT(18)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(2),:)   = U235_3_18_1500_2;
L.xss(L.ZAID(92235),L.MT(18)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(3),:)   = U235_3_18_1500_3;
L.xss(L.ZAID(92235),L.MT(18)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(5),:)   = U235_3_18_1500_5;
L.xss(L.ZAID(92235),L.MT(18)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(10),:)   = U235_3_18_1500_10;

% 92235 gamma
L.xss(L.ZAID(92235),L.MT(102)).NBins(L.Bin(110)).xs(L.T(300),L.BW(0),:)   = U235_3_102_300_0;
L.xss(L.ZAID(92235),L.MT(102)).NBins(L.Bin(110)).xs(L.T(300),L.BW(1),:)   = U235_3_102_300_1;
L.xss(L.ZAID(92235),L.MT(102)).NBins(L.Bin(110)).xs(L.T(300),L.BW(2),:)   = U235_3_102_300_2;
L.xss(L.ZAID(92235),L.MT(102)).NBins(L.Bin(110)).xs(L.T(300),L.BW(3),:)   = U235_3_102_300_3;
L.xss(L.ZAID(92235),L.MT(102)).NBins(L.Bin(110)).xs(L.T(300),L.BW(5),:)   = U235_3_102_300_5;
L.xss(L.ZAID(92235),L.MT(102)).NBins(L.Bin(110)).xs(L.T(300),L.BW(10),:)   = U235_3_102_300_10;

L.xss(L.ZAID(92235),L.MT(102)).NBins(L.Bin(110)).xs(L.T(600),L.BW(0),:)   = U235_3_102_600_0;
L.xss(L.ZAID(92235),L.MT(102)).NBins(L.Bin(110)).xs(L.T(600),L.BW(1),:)   = U235_3_102_600_1;
L.xss(L.ZAID(92235),L.MT(102)).NBins(L.Bin(110)).xs(L.T(600),L.BW(2),:)   = U235_3_102_600_2;
L.xss(L.ZAID(92235),L.MT(102)).NBins(L.Bin(110)).xs(L.T(600),L.BW(3),:)   = U235_3_102_600_3;
L.xss(L.ZAID(92235),L.MT(102)).NBins(L.Bin(110)).xs(L.T(600),L.BW(5),:)   = U235_3_102_600_5;
L.xss(L.ZAID(92235),L.MT(102)).NBins(L.Bin(110)).xs(L.T(600),L.BW(10),:)   = U235_3_102_600_10;

L.xss(L.ZAID(92235),L.MT(102)).NBins(L.Bin(110)).xs(L.T(900),L.BW(0),:)   = U235_3_102_900_0;
L.xss(L.ZAID(92235),L.MT(102)).NBins(L.Bin(110)).xs(L.T(900),L.BW(1),:)   = U235_3_102_900_1;
L.xss(L.ZAID(92235),L.MT(102)).NBins(L.Bin(110)).xs(L.T(900),L.BW(2),:)   = U235_3_102_900_2;
L.xss(L.ZAID(92235),L.MT(102)).NBins(L.Bin(110)).xs(L.T(900),L.BW(3),:)   = U235_3_102_900_3;
L.xss(L.ZAID(92235),L.MT(102)).NBins(L.Bin(110)).xs(L.T(900),L.BW(5),:)   = U235_3_102_900_5;
L.xss(L.ZAID(92235),L.MT(102)).NBins(L.Bin(110)).xs(L.T(900),L.BW(10),:)   = U235_3_102_900_10;

L.xss(L.ZAID(92235),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(0),:)   = U235_3_102_1200_0;
L.xss(L.ZAID(92235),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(1),:)   = U235_3_102_1200_1;
L.xss(L.ZAID(92235),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(2),:)   = U235_3_102_1200_2;
L.xss(L.ZAID(92235),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(3),:)   = U235_3_102_1200_3;
L.xss(L.ZAID(92235),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(5),:)   = U235_3_102_1200_5;
L.xss(L.ZAID(92235),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(10),:)   = U235_3_102_1200_10;

L.xss(L.ZAID(92235),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(0),:)   = U235_3_102_1500_0;
L.xss(L.ZAID(92235),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(1),:)   = U235_3_102_1500_1;
L.xss(L.ZAID(92235),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(2),:)   = U235_3_102_1500_2;
L.xss(L.ZAID(92235),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(3),:)   = U235_3_102_1500_3;
L.xss(L.ZAID(92235),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(5),:)   = U235_3_102_1500_5;
L.xss(L.ZAID(92235),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(10),:)   = U235_3_102_1500_10;

% 92235 mubar
for i = BWs(1:end-1)
    L.xss(L.ZAID(92235),L.MT(251)).NBins(L.Bin(110)).xs(L.T(300),L.BW(i),:)   = U235_3_251_300;
    L.xss(L.ZAID(92235),L.MT(251)).NBins(L.Bin(110)).xs(L.T(600),L.BW(i),:)   = U235_3_251_600;
    L.xss(L.ZAID(92235),L.MT(251)).NBins(L.Bin(110)).xs(L.T(900),L.BW(i),:)   = U235_3_251_900;
    L.xss(L.ZAID(92235),L.MT(251)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(i),:)   = U235_3_251_1200;
    L.xss(L.ZAID(92235),L.MT(251)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(i),:)   = U235_3_251_1500;
end
L.xss(L.ZAID(92235),L.MT(251)).NBins(L.Bin(110)).xs(L.T(300),L.BW(10),:)   = U235_3_251_300_10;
L.xss(L.ZAID(92235),L.MT(251)).NBins(L.Bin(110)).xs(L.T(600),L.BW(10),:)   = U235_3_251_600_10;
L.xss(L.ZAID(92235),L.MT(251)).NBins(L.Bin(110)).xs(L.T(900),L.BW(10),:)   = U235_3_251_900_10;
L.xss(L.ZAID(92235),L.MT(251)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(10),:)   = U235_3_251_1200_10;
L.xss(L.ZAID(92235),L.MT(251)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(10),:)   = U235_3_251_1500_10;

% 92235 nu
for i = BWs
    L.xss(L.ZAID(92235),L.MT(452)).NBins(L.Bin(110)).xs(L.T(300),L.BW(i),:) = U235_3_452_300;
    L.xss(L.ZAID(92235),L.MT(452)).NBins(L.Bin(110)).xs(L.T(600),L.BW(i),:) = U235_3_452_600;
    L.xss(L.ZAID(92235),L.MT(452)).NBins(L.Bin(110)).xs(L.T(900),L.BW(i),:) = U235_3_452_900;
    L.xss(L.ZAID(92235),L.MT(452)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(i),:) = U235_3_452_1200;
    L.xss(L.ZAID(92235),L.MT(452)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(i),:) = U235_3_452_1500;
end

% 92235 nu fission
for j = Ts
    for i = BWs
        L.xss(L.ZAID(92235),L.MT(7)).NBins(L.Bin(110)).xs(L.T(j),L.BW(i),:,:) = ...
            L.xss(L.ZAID(92235),L.MT(452)).NBins(L.Bin(110)).xs(L.T(j),L.BW(i),:).*...
            L.xss(L.ZAID(92235),L.MT(18)).NBins(L.Bin(110)).xs(L.T(j),L.BW(i),:);
    end
end


%% ZAID 92238

% 92238 scatter kernel
for i = BWs
    L.xss(L.ZAID(92238),L.MT(2)).NBins(L.Bin(110)).xs(L.T(300),L.BW(i),:,:) = U238_0_2_300;
    L.xss(L.ZAID(92238),L.MT(2)).NBins(L.Bin(110)).xs(L.T(600),L.BW(i),:,:) = U238_0_2_600;
    L.xss(L.ZAID(92238),L.MT(2)).NBins(L.Bin(110)).xs(L.T(900),L.BW(i),:,:) = U238_0_2_900;
    L.xss(L.ZAID(92238),L.MT(2)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(i),:,:) = U238_0_2_1200;
    L.xss(L.ZAID(92238),L.MT(2)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(i),:,:) = U238_0_2_1500;
end

% 92238 transport

% 92238 fission
L.xss(L.ZAID(92238),L.MT(18)).NBins(L.Bin(110)).xs(L.T(300),L.BW(0),:)   = U238_3_18_300_0;
L.xss(L.ZAID(92238),L.MT(18)).NBins(L.Bin(110)).xs(L.T(300),L.BW(1),:)   = U238_3_18_300_1;
L.xss(L.ZAID(92238),L.MT(18)).NBins(L.Bin(110)).xs(L.T(300),L.BW(2),:)   = U238_3_18_300_2;
L.xss(L.ZAID(92238),L.MT(18)).NBins(L.Bin(110)).xs(L.T(300),L.BW(3),:)   = U238_3_18_300_3;
L.xss(L.ZAID(92238),L.MT(18)).NBins(L.Bin(110)).xs(L.T(300),L.BW(5),:)   = U238_3_18_300_5;
L.xss(L.ZAID(92238),L.MT(18)).NBins(L.Bin(110)).xs(L.T(300),L.BW(10),:)   = U238_3_18_300_10;

L.xss(L.ZAID(92238),L.MT(18)).NBins(L.Bin(110)).xs(L.T(600),L.BW(0),:)   = U238_3_18_600_0;
L.xss(L.ZAID(92238),L.MT(18)).NBins(L.Bin(110)).xs(L.T(600),L.BW(1),:)   = U238_3_18_600_1;
L.xss(L.ZAID(92238),L.MT(18)).NBins(L.Bin(110)).xs(L.T(600),L.BW(2),:)   = U238_3_18_600_2;
L.xss(L.ZAID(92238),L.MT(18)).NBins(L.Bin(110)).xs(L.T(600),L.BW(3),:)   = U238_3_18_600_3;
L.xss(L.ZAID(92238),L.MT(18)).NBins(L.Bin(110)).xs(L.T(600),L.BW(5),:)   = U238_3_18_600_5;
L.xss(L.ZAID(92238),L.MT(18)).NBins(L.Bin(110)).xs(L.T(600),L.BW(10),:)   = U238_3_18_600_10;

L.xss(L.ZAID(92238),L.MT(18)).NBins(L.Bin(110)).xs(L.T(900),L.BW(0),:)   = U238_3_18_900_0;
L.xss(L.ZAID(92238),L.MT(18)).NBins(L.Bin(110)).xs(L.T(900),L.BW(1),:)   = U238_3_18_900_1;
L.xss(L.ZAID(92238),L.MT(18)).NBins(L.Bin(110)).xs(L.T(900),L.BW(2),:)   = U238_3_18_900_2;
L.xss(L.ZAID(92238),L.MT(18)).NBins(L.Bin(110)).xs(L.T(900),L.BW(3),:)   = U238_3_18_900_3;
L.xss(L.ZAID(92238),L.MT(18)).NBins(L.Bin(110)).xs(L.T(900),L.BW(5),:)   = U238_3_18_900_5;
L.xss(L.ZAID(92238),L.MT(18)).NBins(L.Bin(110)).xs(L.T(900),L.BW(10),:)   = U238_3_18_900_10;

L.xss(L.ZAID(92238),L.MT(18)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(0),:)   = U238_3_18_1200_0;
L.xss(L.ZAID(92238),L.MT(18)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(1),:)   = U238_3_18_1200_1;
L.xss(L.ZAID(92238),L.MT(18)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(2),:)   = U238_3_18_1200_2;
L.xss(L.ZAID(92238),L.MT(18)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(3),:)   = U238_3_18_1200_3;
L.xss(L.ZAID(92238),L.MT(18)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(5),:)   = U238_3_18_1200_5;
L.xss(L.ZAID(92238),L.MT(18)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(10),:)   = U238_3_18_1200_10;

L.xss(L.ZAID(92238),L.MT(18)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(0),:)   = U238_3_18_1500_0;
L.xss(L.ZAID(92238),L.MT(18)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(1),:)   = U238_3_18_1500_1;
L.xss(L.ZAID(92238),L.MT(18)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(2),:)   = U238_3_18_1500_2;
L.xss(L.ZAID(92238),L.MT(18)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(3),:)   = U238_3_18_1500_3;
L.xss(L.ZAID(92238),L.MT(18)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(5),:)   = U238_3_18_1500_5;
L.xss(L.ZAID(92238),L.MT(18)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(10),:)   = U238_3_18_1500_10;

% 92238 gamma
L.xss(L.ZAID(92238),L.MT(102)).NBins(L.Bin(110)).xs(L.T(300),L.BW(0),:)   = U238_3_102_300_0;
L.xss(L.ZAID(92238),L.MT(102)).NBins(L.Bin(110)).xs(L.T(300),L.BW(1),:)   = U238_3_102_300_1;
L.xss(L.ZAID(92238),L.MT(102)).NBins(L.Bin(110)).xs(L.T(300),L.BW(2),:)   = U238_3_102_300_2;
L.xss(L.ZAID(92238),L.MT(102)).NBins(L.Bin(110)).xs(L.T(300),L.BW(3),:)   = U238_3_102_300_3;
L.xss(L.ZAID(92238),L.MT(102)).NBins(L.Bin(110)).xs(L.T(300),L.BW(5),:)   = U238_3_102_300_5;
L.xss(L.ZAID(92238),L.MT(102)).NBins(L.Bin(110)).xs(L.T(300),L.BW(10),:)   = U238_3_102_300_10;

L.xss(L.ZAID(92238),L.MT(102)).NBins(L.Bin(110)).xs(L.T(600),L.BW(0),:)   = U238_3_102_600_0;
L.xss(L.ZAID(92238),L.MT(102)).NBins(L.Bin(110)).xs(L.T(600),L.BW(1),:)   = U238_3_102_600_1;
L.xss(L.ZAID(92238),L.MT(102)).NBins(L.Bin(110)).xs(L.T(600),L.BW(2),:)   = U238_3_102_600_2;
L.xss(L.ZAID(92238),L.MT(102)).NBins(L.Bin(110)).xs(L.T(600),L.BW(3),:)   = U238_3_102_600_3;
L.xss(L.ZAID(92238),L.MT(102)).NBins(L.Bin(110)).xs(L.T(600),L.BW(5),:)   = U238_3_102_600_5;
L.xss(L.ZAID(92238),L.MT(102)).NBins(L.Bin(110)).xs(L.T(600),L.BW(10),:)   = U238_3_102_600_10;

L.xss(L.ZAID(92238),L.MT(102)).NBins(L.Bin(110)).xs(L.T(900),L.BW(0),:)   = U238_3_102_900_0;
L.xss(L.ZAID(92238),L.MT(102)).NBins(L.Bin(110)).xs(L.T(900),L.BW(1),:)   = U238_3_102_900_1;
L.xss(L.ZAID(92238),L.MT(102)).NBins(L.Bin(110)).xs(L.T(900),L.BW(2),:)   = U238_3_102_900_2;
L.xss(L.ZAID(92238),L.MT(102)).NBins(L.Bin(110)).xs(L.T(900),L.BW(3),:)   = U238_3_102_900_3;
L.xss(L.ZAID(92238),L.MT(102)).NBins(L.Bin(110)).xs(L.T(900),L.BW(5),:)   = U238_3_102_900_5;
L.xss(L.ZAID(92238),L.MT(102)).NBins(L.Bin(110)).xs(L.T(900),L.BW(10),:)   = U238_3_102_900_10;

L.xss(L.ZAID(92238),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(0),:)   = U238_3_102_1200_0;
L.xss(L.ZAID(92238),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(1),:)   = U238_3_102_1200_1;
L.xss(L.ZAID(92238),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(2),:)   = U238_3_102_1200_2;
L.xss(L.ZAID(92238),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(3),:)   = U238_3_102_1200_3;
L.xss(L.ZAID(92238),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(5),:)   = U238_3_102_1200_5;
L.xss(L.ZAID(92238),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(10),:)   = U238_3_102_1200_10;

L.xss(L.ZAID(92238),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(0),:)   = U238_3_102_1500_0;
L.xss(L.ZAID(92238),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(1),:)   = U238_3_102_1500_1;
L.xss(L.ZAID(92238),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(2),:)   = U238_3_102_1500_2;
L.xss(L.ZAID(92238),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(3),:)   = U238_3_102_1500_3;
L.xss(L.ZAID(92238),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(5),:)   = U238_3_102_1500_5;
L.xss(L.ZAID(92238),L.MT(102)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(10),:)   = U238_3_102_1500_10;

% 92238 mubar
for i = BWs(1:end-1)
    L.xss(L.ZAID(92238),L.MT(251)).NBins(L.Bin(110)).xs(L.T(300),L.BW(i),:)   = U238_3_251_300;
    L.xss(L.ZAID(92238),L.MT(251)).NBins(L.Bin(110)).xs(L.T(600),L.BW(i),:)   = U238_3_251_600;
    L.xss(L.ZAID(92238),L.MT(251)).NBins(L.Bin(110)).xs(L.T(900),L.BW(i),:)   = U238_3_251_900;
    L.xss(L.ZAID(92238),L.MT(251)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(i),:)   = U238_3_251_1200;
    L.xss(L.ZAID(92238),L.MT(251)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(i),:)   = U238_3_251_1500;
end
L.xss(L.ZAID(92238),L.MT(251)).NBins(L.Bin(110)).xs(L.T(300),L.BW(10),:)   = U238_3_251_300_10;
L.xss(L.ZAID(92238),L.MT(251)).NBins(L.Bin(110)).xs(L.T(600),L.BW(10),:)   = U238_3_251_600_10;
L.xss(L.ZAID(92238),L.MT(251)).NBins(L.Bin(110)).xs(L.T(900),L.BW(10),:)   = U238_3_251_900_10;
L.xss(L.ZAID(92238),L.MT(251)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(10),:)   = U238_3_251_1200_10;
L.xss(L.ZAID(92238),L.MT(251)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(10),:)   = U238_3_251_1500_10;

% 92238 nu
for i = BWs
    L.xss(L.ZAID(92238),L.MT(452)).NBins(L.Bin(110)).xs(L.T(300),L.BW(i),:) = U238_3_452_300;
    L.xss(L.ZAID(92238),L.MT(452)).NBins(L.Bin(110)).xs(L.T(600),L.BW(i),:) = U238_3_452_600;
    L.xss(L.ZAID(92238),L.MT(452)).NBins(L.Bin(110)).xs(L.T(900),L.BW(i),:) = U238_3_452_900;
    L.xss(L.ZAID(92238),L.MT(452)).NBins(L.Bin(110)).xs(L.T(1200),L.BW(i),:) = U238_3_452_1200;
    L.xss(L.ZAID(92238),L.MT(452)).NBins(L.Bin(110)).xs(L.T(1500),L.BW(i),:) = U238_3_452_1500;
end

% 92238 nu fission
for j = Ts
    for i = BWs
        L.xss(L.ZAID(92238),L.MT(7)).NBins(L.Bin(110)).xs(L.T(j),L.BW(i),:,:) = ...
            L.xss(L.ZAID(92238),L.MT(452)).NBins(L.Bin(110)).xs(L.T(j),L.BW(i),:).*...
            L.xss(L.ZAID(92238),L.MT(18)).NBins(L.Bin(110)).xs(L.T(j),L.BW(i),:);
    end
end

% TRANSPORT FOR ALL. MAKE SURE WE'RE CAREFUL ABOUT WATER. MUST DEAL WITH
% H1W STILL
for i = ZAIDs
    if i ~= 222 % no no for water
        for j = Ts
            for k = BWs
                L.xss(L.ZAID(i),L.MT(5)).NBins(L.Bin(110)).xs(L.T(j),L.BW(k),:) = ...
                    makeTransport2(squeeze(L.xss(L.ZAID(i),L.MT(18)).NBins(L.Bin(110)).xs(L.T(j),L.BW(k),:)),...
                    squeeze(L.xss(L.ZAID(i),L.MT(102)).NBins(L.Bin(110)).xs(L.T(j),L.BW(k),:)),...
                    squeeze(L.xss(L.ZAID(i),L.MT(2)).NBins(L.Bin(110)).xs(L.T(j),L.BW(k),:,:)),...
                    squeeze(L.xss(L.ZAID(i),L.MT(251)).NBins(L.Bin(110)).xs(L.T(j),L.BW(k),:)));
            end
        end
    end
end

% make 222 transport
L.xss(L.ZAID(222),L.MT(5)).NBins(L.Bin(110)).xs = 2*L.xss(L.ZAID(1001),L.MT(5)).NBins(L.Bin(110)).xs + L.xss(L.ZAID(8016),L.MT(5)).NBins(L.Bin(110)).xs;

Library = L;

end
