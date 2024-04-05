function [kout,phi_it,power]=Spectral_It(V,Sigma,PI,Chi,Region,phi_it)
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
     
            for j=1:1:m_regions  %Source Regions (j)  
            ,    Sigmas(:,:)=permute(Sigma.S(j,:,:),[2,3,1]);
                %Sigmas(:,:)=permute(Sigma.S(j,:,:),[3,2,1]);  %set up sigma's for operations
            ,    %Sum across source Lethargy groups (u')
            ,  if Region(j)==1, %Region returns 1 if it contains fissile material       
            ,     R1(j,:)=(V(j)*Sigmas*(Phi(j,:)'));
            ,     R2(j,:)=V(j)*(Chi)*(Sigma.F(j,:)*(du.*Phi(j,:))'); % /k  (k is left out until later)
            ,  else
            ,     R1(j,:)= V(j)*Sigmas*(Phi(j,:)');   %R(#,u)=V(#)*Sigmas(u,u',#)*(Phi(u',#)
            ,     R2(j,:)= zeros(n_bins,1);
            ,  end;
            end;
            
            for l=1:1:n_bins  %Sink Lethargy group (u)       
               Scatter_sources(:,l)=PI(:,:,l)*R1(:,l);
               Fission_sources(:,l)=PI(:,:,l)*R2(:,l);  
                 %Neutron_sources(i,#)=PI(i,j,#)*R(j,#)
               Neutron_sources(:,l)=Scatter_sources(:,l)+Fission_sources(:,l)/kinf;
               Neutron_sinks(:,l)=V(:).*Sigma.T(:,l).*Phi(:,l); 
                 %Neutron_sinks(i,#)=V(i).*SigmaT(i,#).*Phi(i,#)  %find the sinks in each region
               Phi_sink(:,l)=Neutron_sources(:,l)./(V(:).*Sigma.T(:,l)); 
                 %Phi_source(i,l)=Neutron_sources(i,#)./(V(i).*SigmaT(i,#))
            end;
    
            phinew=Phi_sink; 
            
            Phi=phinew;
          
            if mod(i,10)==1 & i>=40
              kinf=mean(mean(Fission_sources))/mean(mean( Neutron_sinks-Scatter_sources ));
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
power=200e3*1.602e-19*sum(sum((Fission_sources)))/max(max(Phi))/sum(V);   % 200,000keV/Fission*Fissions/s
            %this is not the best method as it uses the nth-1 phi value.

figure
    [xx,yy]=stairs(10.^(-4:0.1:6.9),phi_it(1,:));
    loglog(xx,yy,'g');
    hold on
    [xx,yy]=stairs(10.^(-4:0.1:6.9),phi_it(2,:));
    loglog(xx,yy,'k');
    xlabel(['kinf_{iteration}=' num2str(kout)])
    hold off
    
    
    figure
    [xx,yy]=stairs(10.^(-4:0.1:6.9),phi_it(1,:));
    semilogx(xx,yy,'g');
    hold on
    [xx,yy]=stairs(10.^(-4:0.1:6.9),phi_it(2,:));
    semilogx(xx,yy,'k');
    xlabel(['kinf_{iteration}=' num2str(kout)])
    hold off
    
    %semilogx(10.^(-4:0.1:6.9),phi_it(1,:),'g');
    %hold on
    %semilogx(10.^(-4:0.1:6.9),phi_it(2,:),'k');
    %xlabel(['kinf_{iteration}=' num2str(kout)])
    %title(['kinf=' num2str(kout)]);
    %hold off
    
    
    %Test the iteration scheme
[T,S,F,F1,S1]=iterationtest(Sigma,Phi,Chi,V,PI);
%S1-R1;
%F1-R2;

%pause
%[(T-Neutron_sinks)' , (S-Scatter_sources)', (F-Fission_sources)' ]



%pause

    
%% Test of the iteration scheme.  
%     T=Neutron_sinks
%     S=Scatter_sources
%     F=Fission_sources
%
function [T,S,F,F1,S1]=iterationtest(Sigma,phi,Chi,V,PI)
        
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