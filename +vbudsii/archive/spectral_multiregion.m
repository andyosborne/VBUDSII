function [Neutron_sources, Neutron_sinks, Phi_sink, Scatter_sources,Fission_sources]=spectral(Phi,K)
% This function receives the Volume, Sigma(T,S,F),X,Nu and
% a guess for Phi(i,u).  It returns the Neutron_Sources and Neutron_Sinks
% (right and left side respectively of the spectral equation).  It will
% also return if desired, Phi_sink - the spectra of absorbted neutrons. 

global n_regions m_groups SigmaS SigmaF SigmaT V X PI Region  

    %Variable Initialization
           A=zeros(n_regions,m_groups);

           R=A; Neutron_sources=A; Neutron_sinks=A; Phi_sink=A; 
           Scatter_sources=A; Fission_sources=A; R1=A; R2=A;
     
            for j=1:1:n_regions  %Source Regions (j)  
            ,  Sigmas(:,:)=permute(SigmaS(j,:,:),[3,2,1]);  %set up sigma's for operations
            ,    %Sum across source Lethargy groups (u')
            ,  if Region(j)==1, %Region returns 1 if it contains fissile material       
            ,     R1(j,:)= V(j)*Sigmas*(Phi(j,:)');
            ,     R2(j,:)=V(j)*(X(j,:)')*SigmaF(j,:)*(Phi(j,:)'); % /k  (k is left out until later)
            ,  else
            ,     R1(j,:)= V(j)*Sigmas*(Phi(j,:)');   %R(#,u)=V(#)*Sigmas(u,u',#)*(Phi(u',#)
            ,     R2(j,:)= zeros(m_groups,1);
            ,  end;
            end;
            
            for l=1:1:m_groups  %Sink Lethargy group (u)       
               Scatter_sources(:,l)=PI(:,:,l)*R1(:,l);
               Fission_sources(:,l)=PI(:,:,l)*R2(:,l);  
                 %Neutron_sources(i,#)=PI(i,j,#)*R(j,#)
               Neutron_sources(:,l)=Scatter_sources(:,l)+Fission_sources(:,l)/K;
               Neutron_sinks(:,l)=V(:).*SigmaT(:,l).*Phi(:,l); 
                 %Neutron_sinks(i,#)=V(i).*SigmaT(i,#).*Phi(i,#)  %find the sinks in each region
               Phi_sink(:,l)=Neutron_sources(:,l)./(V(:).*SigmaT(:,l)); 
                 %Phi_source(i,l)=Neutron_sources(i,#)./(V(i).*SigmaT(i,#))
            end;
            
            
            
            
 %All lines below this point are comments on the above code.  They constitute a
 %sample run.
            
    
 %This program will be used to verify the matrix form of the Spectral Equations
T=0; if T==1, %silences this part of the code  
 
 %Global variables 
  global n_regions m_groups SigmaS SigmaF SigmaT V X PI Region  
  n_regions=6;      %cell types
  m_groups=10;      %energy groups 
  SigmaS=rand(n_regions,m_groups,m_groups); %SigmaS(region,lethargy,lethargy')
  SigmaF=rand(n_regions,m_groups);          %SigmaF(region,lethargy')
  SigmaT=rand(n_regions,m_groups);          %SigmaT(region,lethargy)
  V=rand(n_regions,1);                        %volume(region)
  X=rand(n_regions,m_groups);                         %Xi(lethargy)
  PI=rand(n_regions,n_regions,m_groups);    %PI(region, region', lethargy);
  Region=[0,0,1,1,0,0];  % 1 if the region is fissile. %note size(Region)==n_regions  
 
 %Input Variables 
  Phi=rand(n_regions,m_groups);             %Phi(region,lethargy*) 
  K=1;                                      %K_infinity
  
 %Calling routine
  [Neutron_sources, Neutron_sinks, Phi_sink]=spectral(Phi,K); 
 
 %Internal Variables
  A, R, Scatter_sources, Fission_sources, R1, R2 
 
 %Output Variables
  Neutron_sources(region,lethargy);
  Neutron_sinks(region,lethargy);
  Phi_sink(region,lethargy);
 end
 
  
  
  %vector of next neutron interaction
%     Neutron_sinks(region) = V(region)*Sigma(region,lethargy)*Phi(region,lethargy)
%Matrix of creation (fission and scatter) of neutrons 
                                %(do I need to include (n2n)?
%     Neutron_sources(region',lethargy)=V(region)*SigmaS(region,lethargy'->lethargy)*Phi(region,lethargy)
%                     +X(lethargy)/k*V(region)*Nu*SigmaF(region,lethargy')*Phi(region,lethargy)
%     Neutron_to_sinks=Pi(region,region',lethargy)*Neutron_sources(region',
%     lethargy)