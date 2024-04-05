function [Phiout,Kout]=Bench1(MATs,N,phi)
% This function will analyze a very simple crosssection.
% The only cross section is Fission, and it only Fissions at 0.5MeV.



MATs=['U235','U238',' O16', ' H2O', '  H1'];
N=[0.0003867,0.007347,0.015470,0.014500,0];
%Chi=[zeros(6,1);
%     1.23227794e-005;
%     0.000377192647;
%     0.0117409149;
%     0.284172008;
%     0.703697986];



%  [NuSigmaF,SigmaT,SigmaS]=geteqntest1(MATs,N);
%    Chi=[zeros(1,9),1,0]'
%    [Phiout,Kout]=solveeqneig(NuSigmaF,SigmaT,SigmaS,Chi)
%      plot(Phiout(:,10));
%      pause
%      [Phiout,Kout]=solveeqnBING(NuSigmaF,SigmaT,SigmaS,Chi);   
%      
%      
%pause
%        [NuSigmaF,SigmaT,SigmaS]=geteqntest2(MATs,N);
%    Chi=[zeros(1,9),1,0]'
%    [Phiout,Kout]=solveeqneig(NuSigmaF,SigmaT,SigmaS,Chi);
%      plot(Phiout(:,10));
%pause   
        [NuSigmaF,SigmaT,SigmaS]=geteqntest2(MATs,N);
    Chi=[zeros(1,9),1,0]'
    [Phiout,Kout]=solveeqneig(NuSigmaF,SigmaT,SigmaS,Chi);
      plot(Phiout(:,10));
pause
    [Phiout,Kout]=solveeqnBING(NuSigmaF,SigmaT,SigmaS,Chi);   
      %plot(Phiout(:,1)/(Phiout(end,1)))
%% function SOLVE with generalized eigenvectors


function [V,D]=solveeqneig(NuSigmaF,SigmaT,SigmaS,Chi,bins);

%edges=10.^(-4:1:7);   %These two are for energy
%du=diff(edges)';  
du=2.303*ones(11,1);

A=diag(SigmaT)-SigmaS;
B=diag(Chi)*ones(length(Chi))*diag(du.*NuSigmaF);

[V,D] = eig(A,B);

phi=abs(V(:,1));
k=D(1,1);

%phi=rand(11,1);
T1=SigmaT .* phi;
T2=SigmaS * phi;
T3= 2.303*Chi*(NuSigmaF' *phi);
T4=T2+T3*k;  
T5=T4-T1;  %should be 0 if phi is an eigenvector

T6=A*phi-(T1-T2); %should be 0 for any phi
T7=B*phi-T3;  %should be 0 for any phi
Eigentest=[V(:,1), T1,T2,T3,T4,T5,T6, T7]

D(1,1)
%% function GET EQN test 1
function [NuSigmaF,SigmaT,SigmaS]=geteqntest1(MATs,N)
%this is a crappy model

SigmaF=[zeros(1,9),1000,0]';   %SigmaF=[0.3,0,zeros(1,9)]';
SigmaS=diag(ones(1,11));      %SigmaS=diag(0.5*ones(1,11))+ [zeros(10,1), diag(0.5*ones(1,10)); zeros(1,11)];
                               %SigmaS(1,1)=1;
SigmaA=zeros(11,1);            %ones(11,1);
SigmaT=sum(SigmaS,1)'+SigmaF; %+0.0252*SigmaA;
NuSigmaF=2.5*SigmaF;

%% function GET EQN test 2
function [NuSigmaF,SigmaT,SigmaS]=geteqntest2(MATs,N)
%this is a crappy model

SigmaF=[zeros(1,9),1000,0]';   %SigmaF=[0.3,0,zeros(1,9)]';
SigmaS=diag(ones(1,11));       %SigmaS=diag(0.5*ones(1,11))+ [zeros(10,1), diag(0.5*ones(1,10)); zeros(1,11)];
SigmaS(3,10)=0.5 %2e-6;
SigmaS(10,10)=0.5;                               %SigmaS(1,1)=1;
SigmaA=zeros(11,1);            %ones(11,1);
SigmaT=sum(SigmaS,1)'+SigmaF; %+0.0252*SigmaA;
NuSigmaF=2.5*SigmaF;


%% Bing's solver
function [V,D]=solveeqnBING(NuSigmaF,SigmaT,SigmaS,Chi);

phi=rand(11,1)
kinf=1;

n=700; tic;
for i=1:1:n
    phinew=( SigmaS * phi + Chi/kinf*((NuSigmaF'*phi)*2.303))./SigmaT;
    phi=phinew;%/sum(phinew);
    if mod(i,10)==1 & i>=40
    kinf=mean((Chi*(NuSigmaF'*phi)*2.303))/mean( SigmaT.*phi-SigmaS*phi );
    end
    
 %   if max(phinew-phi)<=1e-5, ['loop done in ', int2str(i), ' iterations: ', int2str(toc*1000), ' milliseconds' ],  break, end
%end
    if i==1, growth(i)=sum(phinew);
    else, growth(i)=sum(phinew)*growth(i-1);
    end;
%test=SigmaT.*phinew-( SigmaS * phinew + Chi/kinf*(NuSigmaF'*phinew))
if mod(i,1)==1
plot(phi,'r');
hold on
plot(phinew);
hold off
pause(0.5)
kout=kinf;
end
end
V=phi;
D=kout;
kinf
figure, plot(growth)

%% function GET EQN
function [NuSigmaF,SigmaT,SigmaS]=geteqn(MATs,N)

MATs=['U235','U238',' O16', ' H2O', '  H1'];
%fissile first.
%load data
load water
%load falsewater


%% Fission
  %nu*sigmaF is actually nu*sigmaF
  for i=1:1:2 %2 fissile materials
     eval([ 'Sf(:,i)=N(i)*' MATs(4*(i-1)+1:4*i) '_3_18_600' '.*' MATs(4*(i-1)+1:4*i) '_3_452_600;' ]);
  end
      NuSigmaF= sum(Sf,2); 

%% Total
  %sigmaF is actually nu*sigmaF
  for i=1:1:5 %5 materials
     eval([ 'St(:,i)=N(i)*' MATs(4*(i-1)+1:4*i) '_3_1_600;' ]);
  end
      SigmaT= sum(St,2); 

        
%% Scatter
     %SigmaS is actually SigmaS*binwidth
  for i=1:1:5 %2 fissile materials
     eval([ 'Ss(:,:,i)=N(i)*' MATs(4*(i-1)+1:4*i) '_0_2_600;' ]);
  end
      SigmaS= sum(Ss,3); 

%% SIGMAT error analysis
  for i=1:1:5 %5 materials
     eval([ 'Sg(:,i)=N(i)*' MATs(4*(i-1)+1:4*i) '_3_102_600;' ]);
  end
  for i=1:1:3 %3 inelastic scattering materials
     eval([ 'Sie(:,i)=N(i)*' MATs(4*(i-1)+1:4*i) '_3_4_600;' ]);
  end
  for i=1:1:2 %2 fissile materials
     eval([ 'Sf(:,i)=N(i)*' MATs(4*(i-1)+1:4*i) '_3_18_600;' ]);
  end

   SigmaF=sum(Sf,2);
   SigmaG=sum(Sg,2);
   SigmaI=sum(Sie(2:end,:),2);   %needs to be modified when I correct library
   SigmaE=sum(SigmaS,1)';
  %SigmaT=SigmaG+SigmaI+SigmaE+SigmaF;
  
  handle1=figure;
   hold off
   semilogy(SigmaT,'k'), hold on
   semilogy(SigmaF,'g')
   semilogy(SigmaG,'c')
   semilogy(SigmaI,'y')
   semilogy(SigmaE,'b')
   %pause,
   semilogy(SigmaG+SigmaI+SigmaE+SigmaF,'r')
   SigmaError_percent=(SigmaT-(SigmaG+SigmaI+SigmaE+SigmaF))./SigmaT
  close(handle1)
%%   
save sigma
