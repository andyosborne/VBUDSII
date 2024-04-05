function [PI_Transmission_Probabilities]=Create_Pi(Sigma, p)
%% File Description
% This file will make use of the Geometry, and Transmission crosssections
% to formulate the Pi Matricies.  At this point we are using Sauer's
% formula.
%
%  Assumes
%     Sigma.R(i,u)
%
%   Geometries needed (Assume pin geometry )
%       p - pin pitch
%       d - pin diameter
%       f -
%       g -
%
%   Transmission cross sections for each region
%
%   Updates needed
%       1. Deal with single transmission cross sections
%       2. Deal with alternate geometries
%            It might be worth making a module that takes the geometry and
%            returns the appropriate x and c values for Sauer's.
%            The lines that depend on pin geometry are marked with ---->  %(PIN)
%
%

if p.verbose
    disp('Entering Create_Pi');
end

d=[1,1]; %[0.2,0.2]; %pin diameter
pp=[2,2]; %[0.3,0.3]; %pin pitch
f=1; %[1 2 3 4 5]/15;
n_mod_cells=length(f);
gg=1; %[4 3 5 6]/18;
n_fuel_cells=length(gg);

n_tot_cells=n_fuel_cells+n_mod_cells; %length(Sigma.R(:,1));

if length(d)==1,
    d=d*ones(1,n_mod_cells);
    pp=pp*ones(1,n_mod_cells); 'warning if # fuel > # mods';
end

%% Get Average T values

%Moderator regions (ie anti-cylinders)
xm=zeros(n_mod_cells,p.nFineGroups); %length(Sigma.R));
Tm=xm; Pm=xm;
xf=zeros(n_fuel_cells,p.nFineGroups);      %length(Sigma.R));
Tf=xf; Pf=xf;


c=2.35; %Constant for Sauer's formula for a moderator                     %(PIN)
for i=1:1:n_mod_cells
    xm(i,:)=Sigma.R(i,:)*d(i)*(4*pp(i)^2/(pi*d(i)^2)-1);                       %(PIN)
    Tm(i,:)=(1./(1+xm(i,:)/c).^c);  %columns are energy
    Pm(i,:)=(1-Tm(i,:))./xm(i,:);
end
Tavg_m=(f*Tm);  %returns the average transmission coeff in the moderator
%as a row vector in energy.

%Fuel regions (cylinder regions)
c=5.00; %Constant for Sauer's formula for a cylinder fuel pin             %(PIN)
for i=1:1:n_fuel_cells
    index=i+n_mod_cells;  %moderator regions are first, then Fuel regions
    xf(i,:)=Sigma.R(index,:)*d(index);                                        %(PIN)
    Tf(i,:)=(1./(1+xf(i,:)/c).^c);  %columns are energy
    Pf(i,:)=(1-Tf(i,:))./xf(i,:);
end
Tavg_f=(gg*Tf);  %returns the average transmission coeff in the fuel
%as a row vector in energy.

tau_m=(1-Tavg_m)./(1-Tavg_m.*Tavg_f);   %Moderator collision probability
tau_f=(1-Tavg_f)./(1-Tavg_m.*Tavg_f);   %Fuel collision probability

%Checks
if abs(1-sum(gg))>=0.005, 'ERROR: g values do not sum to 1 - in Create_Pi', end
if abs(1-sum(f))>=0.005, 'ERROR: f values do not sum to 1 - in Create_Pi', end

%% Create Pi Matrix terms

%   [       |          ]   This matrix has the form PI(i,j,u), where the u
%   [  M1   |    M4    ]   is not visible.  It should be noted that M1 & M2
%   [_______|__________]   are square. M3&M4 represent the transmission
%   [       |          ]   between mod->fuel, and fuel->mod respectively.
%   [  M3   |    M2    ]
%   [       |          ]   The convention on PI is PI( i <- j , u ).  Thus,
%   [       |          ]   sums are to be on i and u & contractions on j.

PI=ones(n_tot_cells, ...
    n_tot_cells, ...
    p.nFineGroups);

for i=1:1:p.nFineGroups
    M1=diag(1-Pm(:,i)) + diag(f)*ones(length(f))*diag(Pm(:,i))*Tavg_f(i)*tau_m(i);
    M2=diag(1-Pf(:,i)) + diag(gg)*ones(length(gg))*diag(Pf(:,i))*Tavg_m(i)*tau_f(i);
    M3=diag(gg)*ones(length(gg),length(f))*diag(Pm(:,i))*tau_f(i);
    M4=diag(f)*ones(length(f),length(gg))*diag(Pf(:,i))*tau_m(i);
    PI(:,:,i)=[M1, M4 ; M3, M2];
end

PI_Transmission_Probabilities=PI;

return

% P11(1,:)=PI(1,1,:)
% P21(1,:)=PI(2,1,:)
% P12(1,:)=PI(1,2,:)
% P22(1,:)=PI(2,2,:)
% 
% plot(P11,'k')
% hold
% plot(P12,'gg')
% plot(P22,'b')
% plot(P21,'r')
% pause
% 
% [Pm' , Pf']
% [Tavg_f', Tavg_m', tau_m', tau_f']
% 
% pause
% 
