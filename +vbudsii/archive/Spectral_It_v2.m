function [kout,phi_it]=Spectral_It(V,Sigma,PI,Chi,Region,phi_it)
    %  This section of code calls Spectral_It.m
    %  Spectral solve takes arguments V,Sigma,PI,phi
    %  It takes the arguments and iterates to find Phi_new and k_inf.
    
[m_regions,n_bins]=size(Sigma.T);
Phi=phi_it;
kinf=1;
du=2.303*ones(1,n_bins);  'what is du?';

n=700;
for i=1:1:n
    
           A=zeros(m_regions,n_bins);
           R=A; Neutron_sources=A; Neutron_sinks=A; Phi_sink=A; 
           Scatter_sources=A; Fission_sources=A; R1=A; R2=A;
     
              [T,S,F]=iterationtest(Sigma,Phi,Chi,V,PI);
            
            for l=1:1:n_bins  %Sink Lethargy group (u)       
               N(:,l)=S(:,l)+F(:,l)/kinf;
                 %Neutron_sinks(i,#)=V(i).*SigmaT(i,#).*Phi(i,#)  %find the sinks in each region
               Phi_sink(:,l)=N(:,l)./(V(:).*Sigma.T(:,l)); 
                 %Phi_source(i,l)=Neutron_sources(i,#)./(V(i).*SigmaT(i,#))
            end;
    
            phinew=Phi_sink; 
            
            Phi=phinew;
          
            if mod(i,10)==1 && i>=40
              kinf=mean(mean(F))/mean(mean( T-S ));
            end    
    
            
            
            
%if max(phinew-phi)<=1e-5, ['loop done in ', int2str(i), ' iterations: ', int2str(toc*1000), ' milliseconds' ],  break, end
%end
%test=SigmaT.*phinew-( SigmaS * phinew + Chi/kinf*(NuSigmaF'*phinew))
if mod(i,1)==1
    'ran a crazy code'
plot(Phi,'r');
hold on
plot(Phinew);
hold off
pause(0.5)
kout=kinf;
end
end


kout=kinf;
phi_it=Phi/max(max(Phi));

%figure
    semilogx(10.^(-4:0.1:6.9),phi_it(1,:),'g');
    hold on
    semilogx(10.^(-4:0.1:6.9),phi_it(2,:),'k');
    xlabel(['kinf_{iteration}=' num2str(kout)])
    %title(['kinf=' num2str(kout)]);
    hold off
    



    
%% Test of the iteration scheme.  
%     T=Neutron_sinks
%     S=Scatter_sources
%     F=Fission_sources
%
function [T,S,F]=iterationtest(Sigma,phi,Chi,V,PI)
        
    for u=1:1:110
        T(:,u)= V(:).*Sigma.T(:,u).*phi(:,u);
    end

    for j=1:1:2
        F1(j,:)=Chi(:) * (V(j) *( Sigma.F(j,:) * phi(j,:)' )* 2.303);
    end
    
    for u=1:1:110  
        for j=1:1:2
            Sigmas(1,:)=Sigma.S(j,u,:);
            S1(j,u)= (V(j) *( Sigmas(1,:)) * phi(j,:)' );
        end
    end

%multiplication by PI
    for u=1:1:110  
        for i=1:1:2
            Pi(1,:)=PI(i,:,u);
            F(i,u)= ( Pi(1,:) * F1(:,u) );
            S(i,u)= ( Pi(1,:) * S1(:,u) );
        end
    end