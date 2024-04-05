function [outputs]=frontend(ratios,per,inputs)
% This function interacts with the GUI

MATs=  ['U235','U238',' O16', ' H2O', '  H1'];
density = [10.3, 10.3, 10.3, 1.0, 0.01];  %[g/cm3];
lambda = [7.038e8*Year,4.468e9*Year, Stable, Stable, Stable];
Atomic_Mass=[235,238,16,18,1];
%Halflife time scales (in seconds)
uSec =1e-6;
mSec =1e-3;
Sec  =1;
Min  =60;
Hour =3600;
Day  =86400;
Week =604000;
Month=2629743.83;
Year =31556926;
Stable=inf;

lambda = []; %halflife in [seconds]
%This list should have some common molecules
%Molecules=[' UO2',' H20',' ZO2', ' etc'];
%natural=['Uran','carb','


f=[1]; %Cell fraction moderator region
g=[1]; %Cell fraction Fissile region

f=[0.7002, 0.2998]; %Cell fraction moderator region
g=[0.6473, 0.2739, 0.0498, 0.0249, 0.0031, 0.0010]; %Cell fraction Fissile region



% Calculate N, the (# of atoms)/b/cm



% Atomic_density(i)=0.6022*density(i)/Atomic_Mass(i);  %[# atoms / b·cm]





if per=2 %User input a ratio: MAT% by volume
    N(i)=ratios(i)*0.6022*density(i)/Atomic_Mass(i);  %[# atoms / b·cm]
elseif per=1  %User input a ratio: MAT% by mass
    N(i)=ratios(i)*0.6022/Atomic_Mass(i);             %[# atoms / b·cm]   
elseif per=0  %User input a ratio: Atomic%
    
    
    
    
    N(i)=0.6022*density(i)/atomicN(i);  %[# atoms / b·cm]
    