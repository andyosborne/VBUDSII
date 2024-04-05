function [Scatter_sources,Fission_sources,Neutron_sinks]=Spectral_matrix_t(V,Sigma,PI,Chi,Region,Phi,test);
     % This program receives the same inputs as the Spectal_matrix.m and is
     % here to verify that the Spectral_matrix.m file returns the correct
     % answers.
     %
     %
% This function receives the Volume, Sigma(T,S,F),X,Nu and
% a guess for Phi(i,u).  It returns the Neutron_Sources and Neutron_Sinks
% (right and left side respectively of the spectral equation).  It will
% also return if desired, Phi_sink - the spectra of absorbted neutrons. 

K=1;
[m_regions,n_bins]=size(Sigma.T); 

    %Variable Initialization
           A=zeros(m_regions,n_bins);

           R=A; Neutron_sources=A; Neutron_sinks=A; Phi_sink=A; 
           Scatter_sources=A; Fission_sources=A; R1=A; R2=A;
     
            for j=1:1:m_regions  %Source Regions (j)  
              Sigmas(:,:)=permute(Sigma.S(j,:,:),[2,3,1]);  %set up sigma's for operations
                 %Sum across source Lethargy groups (u')
              if Region(j)==1, %Region returns 1 if it contains fissile material       
                 R1(j,:)= V(j)*Sigmas*(Phi(j,:)');
                 R2(j,:)=V(j)*(Chi)*Sigma.F(j,:)*(Phi(j,:)'); % /k  (k is left out until later)
              else
                 R1(j,:)= V(j)*Sigmas*(Phi(j,:)');   %R(#,u)=V(#)*Sigmas(u,u',#)*(Phi(u',#)
                 R2(j,:)= zeros(n_bins,1);
              end;
            end;
            
            for l=1:1:n_bins  %Sink Lethargy group (u)       
               Scatter_sources(:,l)=PI(:,:,l)*R1(:,l);
               Fission_sources(:,l)=PI(:,:,l)*R2(:,l);  
                 %Neutron_sources(i,#)=PI(i,j,#)*R(j,#)
               Neutron_sources(:,l)=Scatter_sources(:,l)+Fission_sources(:,l)/K;
               Neutron_sinks(:,l)=V(:).*Sigma.T(:,l).*Phi(:,l); 
                 %Neutron_sinks(i,#)=V(i).*SigmaT(i,#).*Phi(i,#)  %find the sinks in each region
               Phi_sink(:,l)=Neutron_sources(:,l)./(V(:).*Sigma.T(:,l)); 
                 %Phi_source(i,l)=Neutron_sources(i,#)./(V(i).*SigmaT(i,#))
            end;

 %All lines below this point are comments on the above code.  They constitute a
 %sample run.
            
    
 %This program will be used to verify the matrix form of the Spectral Equations
T=0; if T==1, %silences this part of the code  
 
 %Global variables 
 % global m_regions n_bins SigmaS SigmaF SigmaT V X PI Region  
  m_regions=6;      %cell types
  n_bins=10;      %energy groups 
  Sigma.S=rand(m_regions,n_bins,n_bins); %Sigma.S(region,lethargy,lethargy')
  Sigma.F=rand(m_regions,n_bins);          %Sigma.F(region,lethargy')
  SigmaT=rand(m_regions,n_bins);          %SigmaT(region,lethargy)
  V=rand(m_regions,1);                        %volume(region)
  X=rand(m_regions,n_bins);                         %Xi(lethargy)
  PI=rand(m_regions,m_regions,n_bins);    %PI(region, region', lethargy);
  Region=[0,0,1,1,0,0];  % 1 if the region is fissile. %note size(Region)==m_regions  
 
 %Input Variables 
  Phi=rand(m_regions,n_bins);             %Phi(region,lethargy*) 
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
%     Neutron_sources(region',lethargy)=V(region)*Sigma.S(region,lethargy'->lethargy)*Phi(region,lethargy)
%                     +X(lethargy)/k*V(region)*Nu*Sigma.F(region,lethargy')*Phi(region,lethargy)
%     Neutron_to_sinks=Pi(region,region',lethargy)*Neutron_sources(region',
%     lethargy)