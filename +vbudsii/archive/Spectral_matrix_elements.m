function [NuSigmaF,Sigma,Mubar,Bateman_matrix]=Spectral_matrix_elements(MATs,N,Isfissile,p,L)
%% File Description
% This function will use the library to create the crosssection library for
% a given region of the reactor.
%
%  Standard Syntax:
%     [NuSigmaF,Sigma,Mubar,Bateman_matrix]=spectral_matrix_elements(MATs,N,Reactor_Temp);
%
%
% Input :
%         MATs = ['U235','U238',' O16', ' H2O', '  H1'];
%         Materials as defined by their 4 char code (Spaces for shorter
%         characters.  This is just a list of the possible materials in the
%         region.
%
%         N = [0.0003867,0.007347,0.015470,0.014500,0];
%         Atomic Density of each material as measured in units of
%         [atoms/b/cm].  It is worth noting that some bound molecules are
%         treated like atoms (H2O,D2O,ZrO, etc) and thus the units would be
%         [Molecules/b/cm].
%
%         Reactor_Temp = 615;
%         Temperature of the region being analyzed.  The Data is collected
%         in specific temperatures (300,600,900,1200,1500) and linear
%         interpolation is used to find the data at that point.
%
% Output:
%         NuSigmaF [# of neutrons / Fission * Fission / b]
%               Fission cross section * Nu
%
%         Sigma.&  [1/cm]
%               Cross sections (Macroscopic)
%                   T-Total
%                   E-Elastic
%                   I-Inelastic
%                   G-Absobtion (N,gamma)
%                   F-Fission
%                   R-Transport
%               NOTE:  Sigma.T= E+I+G+F
%                      Sigma.R= T-Mubar*E
%
%          Mubar   [unitless]
%               Average Cosine for the elastic scattering of neutrons.
%
%          Sigma.S [1/cm]
%               Elastic Scattering Kernel, contains du', the bin size of
%               the sources.  sum(Sigma.S,1)=Sigma.E
%
% Operations:
%          Collection of crosssections
%          Microscopic Xsection -> Macroscopic Xsection
%          Average material Xsections
%          Average Xsections based on Temperature
%
%          Sigma.T Calculated from others
%          Sigma.R Calculated from others
%
%
% STILL NEEDED - OUTPUT: Inelastic Scattering Cross section.
%                OUTPUT: BATEMAN MATRIX
%                MODIFY: Bonderenko weighting needs to be considered
%                MODFIY: H2O and other molecules need to be taken at
%                        similar tempatures (otherwise we error).
%
%
%                         ENDF LIBRARY
%                             |
%                             V
%                     ___________________
%                     |                 |
%             MATs -->|                 | --> Sigma.&      [1/cm]
%         (materials) |                 |     Crosssections
%                     |    Spectral     |
%             N ----->|     Matrix      |
%         (atomic     |    Elements.m   | --> NuSigmaF     [#/cm]
%           density)  |                 |
%                     |                 | --> Sigma.S      [1/cm]
%      Reactor_Temp ->|                 |     ScatteringKernel * du'
%                     |                 |
%                     |_________________|
%
%
%
% These crosssections will be compiled by Spectral_matrix.m into a large
% generalized eigenvalue equation and then solved by Matlab's eigenvalue
% solver.

%These are just starters for testing the runs.  They will be removed.

if p.verbose
    disp('Entering Spectral_matrix_elements');
end

%% Input data from Library
lib = CrossSectionLibrary.soleInstance();
load water %contains all data. %  WARNING: AT THIS POINT WATER CONTAINS SOME EXTRA DATA VALUES, FOR INSTANCE, TIME

[regionnum,matnum]=size(N);
%% function GET EQN for a specific  region

%Determine how temperature selection fits into data
if p.temp<=L.Ts(1),
    temp=L.Ts(1);
    %['Warning: Temperature selection is below ' num2str(L.Ts(1)) 'k, using ' num2str(L.Ts(1)) 'k instead.']
elseif p.temp>=L.Ts(end),
    temp=L.Ts(2);
    %['Warning: Temperature selection is above ' num2str(L.Ts(end)) 'k, using ' num2str(L.Ts(end)) 'k instead.']
elseif find(abs(L.Ts-p.temp)<=5),
    temp=L.Ts( find(L.Ts > p.temp-5,1) );
    %['Proximity notice: Since the temperature is within 5k of ' num2str(temp) 'k we will use ' num2str(temp) 'k.']
else
    hi_temp_index=find(L.Ts > p.temp,1);
    temp=[L.Ts(hi_temp_index-1),L.Ts(hi_temp_index)];
end;

for l=1:2, %linear interpolation between multiple temperatures.
    Temp=temp(l);
    
    %Create Fission Cross section
    i=0; %start counter
    for j=1:1:matnum 
        %Collect Fission cross sections from each material in the region.
        if N(j)==0 || Isfissile(j)==0, %then do nothing
        else %create Fission cross section
            i=i+1;   %sigmaF is actually nu*sigmaF
            Sf(:,i) = N(j) * lib.xsSpectrum(MATs{j}, 'fission', Temp);
        end
    end
    if i==0,
        NuSigmaF=zeros(p.nFineGroups,1);
        Sigma.F=zeros(p.nFineGroups,1);
    else
        NuSigmaF= sum(Sf,2); %add cross sections for different elements.
        Sigma.F=sum(Sf,2);
    end
    
    
    %Sie(:,1)=zeros(str_length,1); %Some regions will have no inelastic scattering.
    
    %Create other Cross sections
    i=0; %create a counter
    %preallocate
    St = zeros(p.nFineGroups,matnum);
    Ss = zeros(p.nFineGroups,p.nFineGroups,matnum);
    Sg = zeros(p.nFineGroups,matnum);
    Mu = zeros(p.nFineGroups,matnum);
    Sie = zeros(p.nFineGroups,matnum);
    for j=1:1:matnum   %Collect other material properties from each material in the region.
        if N(j)==0, %then do nothing
        else
            i=i+1; %create cross section
            St(:,i) = N(j)*lib.xsSpectrum(MATs{j}, 'total', Temp);
            Ss(:,:,i) = N(j)*lib.elScatKernel(MATs{j}, Temp);
            Sg(:,i) = N(j)*lib.xsSpectrum(MATs{j}, 'nGamma', Temp);
            Mu(:,j) = N(j)*lib.muScat(MATs{j}, Temp);
            Sie(:,i) = N(j)*lib.xsSpectrum(MATs{j}, 'inelastic', Temp);
        end
    end
    
    
    %Combine elements.
    %Sigma.T=sum(St,2);  %Total Cross Section CHRIS DETERMINED THIS UNNEC.
    Sigma.S=sum(Ss,3);  %Scattering Kernel
    Sigma.G=sum(Sg,2);  %(n,Gamma) cross section
    Sigma.I=sum(Sie,2); %Inelastic Scattering cross section  %need correct for 10 group poorly formed lib.
    Sigma.E=sum(Sigma.S,1)'; %Elastic Scattering cross section (from Kernel).
    Sigma.T=Sigma.G+Sigma.I+Sigma.E+Sigma.F;
    Mubar=sum(Mu,2);   %Mubar - Average cosine of the scattering angle
    Sigma.R=Sigma.T-Mubar.*Sigma.E;
    
    
    %[Sigma.R,Sigma.T,Sigma.E,Sigma.G,Sigma.I,Mubar] %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %[Sigma.T, Sigma.E, Sigma.R, Mubar]
    
    %SigmaError_percent=(SigmaT-(SigmaG+SigmaI+SigmaE+SigmaF))./SigmaT
    
    
    if length(temp)==1, break,
    elseif l==1,
        SigmaLO=Sigma;
        MubarLO=Mubar;
        NuSigmaFLO=NuSigmaF;
        
    elseif l==2,
        SigmaHI=Sigma;
        MubarHI=Mubar;
        NuSigmaFHI=NuSigmaF;
        
        %Linear interpolation
        Factor=(p.temp-temp(1))/(temp(2)-temp(1));  %1 lo, 2 hi
        Mubar=MubarLO*(1+Factor)-Factor*MubarHI;
        for j=char([fieldnames(Sigma)])',
            eval([ 'Sigma.' j '=SigmaLO.' j '*(1+Factor)-Factor*SigmaHI.' j ';' ]);
        end
        
        %Sigma=SigmaLO.ii*(1+Factor)-Factor*SigmaHI;
        NuSigmaF=NuSigmaFLO*(1+Factor)-Factor*NuSigmaFHI;
    end;  %
    
end

Bateman_matrix='this code needs to be written'; %This needs to be completed

%save Library1

