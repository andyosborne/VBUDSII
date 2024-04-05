function [Phiout,Kout]=spectral(MATs,N,phi)
MATs=['U235','U238',' O16', ' H2O', '  H1'];

%handle=figure;
%hold on, n=0;
%for a=0.1:0.1:100,
%    n=n+1;
a=.1; b=1;
N=b*[0.0003867*a,0.007347*a,0.015470*a,0.014500,0];
%N=[0.05,0.95,0,0,0];
%phi=[1 2 3 2 1 1 1 1 2 3 2]';
%phi=rand(1,11)';
Chi=[zeros(6,1);
     1.23227794e-005;
     0.000377192647;
     0.0117409149;
     0.284172008;
     0.703697986];

 load water
Chi=Chi'


[NuSigmaF,SigmaT,SigmaS]=geteqn(MATs,N);

%[NuSigmaF,SigmaT,SigmaS]=geteqntest(MATs,N);
%Chi=[zeros(1,9),0,1]'
    [Phiout,Kout]=solveeqneig(NuSigmaF,SigmaT,SigmaS,Chi);
    %figure(handle)
    %plot(Phiout(:,1)/(Phiout(end,1)),'r')
   % K(n)=Kout(1,1);
%    Ratio(n)=a;
%end
%    figure, plot(Ratio,K,'*-')
%% function SOLVE with generalized eigenvectors


function [V,Kout]=solveeqneig(NuSigmaF,SigmaT,SigmaS,Chi,bins);

%edges=10.^(-4:1:7);   %These two are for energy
%du=diff(edges)';  
du=2.303*ones(length(Chi),1);

A=diag(SigmaT)-SigmaS;
B=diag(Chi)*ones(length(Chi))*diag(du.*NuSigmaF);

[V,D] = eig(A,B);

phi=abs(V(:,1));
k=D(1,1);

%phi=rand(11,1);
T1=SigmaT .* phi;
T2=SigmaS * phi
T3= 2.303*Chi*(NuSigmaF' *phi)
T4=T2+T3*k;  
T5=T4-T1;  %should be 0 if phi is an eigenvector

T6=A*phi-(T1-T2); %should be 0 for any phi
T7=B*phi-T3;  %should be 0 for any phi
Eigentest=[V(:,1), T1,T2,T3,T4,T5,T6, T7];

Kout=1./diag(D);
%% function GET EQN test
function [NuSigmaF,SigmaT,SigmaS]=geteqntest(MATs,N)
%this is a crappy model




% SigmaF=[1000,zeros(1,10)]';
% SigmaS=diag(0.5*ones(1,11))+ [zeros(10,1), diag(0.5*ones(1,10)); zeros(1,11)];
% SigmaS(1,1)=1;
% SigmaA=ones(11,1);
% sum(SigmaS,1)';
% SigmaT=sum(SigmaS,1)'+SigmaF+0.09565*SigmaA;
% NuSigmaF=2.5*SigmaF;

SigmaF=[0.3,0,zeros(1,9)]';
SigmaS=diag(0.5*ones(1,11))+ [zeros(10,1), diag(0.5*ones(1,10)); zeros(1,11)];
SigmaS(1,1)=1;
SigmaA=ones(11,1);
sum(SigmaS,1)';
SigmaT=sum(SigmaS,1)'+SigmaF; %+0.0252*SigmaA;
NuSigmaF=2.5*SigmaF;

%% function GET EQN
function [NuSigmaF,SigmaT,SigmaS]=geteqn(MATs,N)

MATs=['U235','U238',' O16', ' H2O', '  H1'];

load water



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
   SigmaI=sum(Sie,2); 
   SigmaE=sum(SigmaS,1)';
  %SigmaT=SigmaG+SigmaI+SigmaE+SigmaF;
  
 % handle1=figure;
 %  hold off
 %  semilogy(SigmaT,'k'), hold on
 %  semilogy(SigmaF,'g')
 %  semilogy(SigmaG,'c')
 %  semilogy(SigmaI,'y')
 %  semilogy(SigmaE,'b')
 %  %pause,
 %  semilogy(SigmaG+SigmaI+SigmaE+SigmaF,'r')
   SigmaError_percent=(SigmaT-(SigmaG+SigmaI+SigmaE+SigmaF))./SigmaT 
 %  close(handle1)
%%   
save sigma
