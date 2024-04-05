function [R,F,L,Rform,Fform,Lform]=Spectral_matrix(p,r,Sigma,PI,Chi,Region,Phi)
% this function will eventually become part of the Spectral Solv function.
%
% Assumed input forms.
%   Sigma.T(i,u)
%   Sigma.S(j,u,u')   (u <- u')
%   Sigma.F(j,u')
%   Chi(u)
%   PI(i,j,u)         (i <- j )
%   r.relVolumes(j)
%
%   j=1:1:r.nCells - spatial
%   u=1:1:p.nFineGroups    - energy
%
%   Regions(j)      - 1 if fissile matrial
%                   - 0 if no fissile material
%
%   phi input can be removed later
%    phi(i,u)
%
%  It returns R F L which look like
%  L*phi=R*phi+1/k_inf*F*phi
%  and can be used in the generalized eigenvector solver.
%
%  It can also return Rform,Fform,Lform where
%   Lform=L*phi;  -> reshaped into a Lform(i,u)
%
%  This line can be erased.
%   global n_regions m_groups SigmaS SigmaF SigmaT V X PI Region

K=1;

% spectral matrix equations

%% L solution (Neutron_sinks)
%This line of code uses SigmaTotal and the volume vector to create a
%matrix that operates on Phi to get the LHS of the spectral equation.
for u=1:1:p.nFineGroups
    LL(:,u)=r.relVolumes(:).*Sigma.T(:,u);
end
[i,u]=size(LL);
L=diag(reshape(LL,(i*u),1));
%this has been checked and evaluates properly.  It assumes that phi is
%going to be written as a column vector in the form phi(1,1) (2,1) ...
%phi(i,1) ...  phi(i,2) ... phi(i, u).

%% R solution (Neutron_sources - due to scatter)
% This line of code uses SigmaScattering Kernel and a volume vector to
% create a matrix that operates on Phi.
% It assumes that Sigma_Scatter is input in the form SigmaS(j,u <- u')

for cellidx = 1:r.nCells
    R3(cellidx,:,:)=r.relVolumes(cellidx)*Sigma.S(cellidx,:,:);
end;

R2=size32(R3,[3,1,2],1);
P1=size32(PI,[2,3,1],2);  %Probably better to stick T*PI into size32

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%I1 will make R(i,u)=>R(u,i)  %can I
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%use the function reshape or rot here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%maybe transpose?
I1=zeros(r.nCells*p.nFineGroups);
for j=1:1:(r.nCells*p.nFineGroups)
    rr=((j-mod(j-1,r.nCells))-1)/r.nCells;
    I1(j,(j-1)*p.nFineGroups+1+rr-rr*r.nCells*p.nFineGroups)=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R=I1*P1*R2;
%R_old=P1*T1;

%% F solution (Neutron_sources - due to Fission)

du=2.303*ones(p.nFineGroups,1);  'what is du?'; % lethargy constant

for j=1:1:r.nCells
    if Region(j)==1,
        for up=1:1:p.nFineGroups %(u')
            F3(j,:,up)=r.relVolumes(j)*Chi(:)*du(up)*Sigma.F(j,up);
        end;
    elseif Region(j)==0,
        F3(j,:,:)=zeros(p.nFineGroups,p.nFineGroups);
    else, 'ERROR bad Region encoding in region j', j
    end,
    
end;
F2=size32(F3,[3,1,2],1);
F=I1*P1*F2;

%   F31(:,:)=F3(1,:,:);
%   F32(:,:)=F3(2,:,:);
%
%   F31-2*F32
%   max(max(ans))
%   pause

%% (n,2n) or other absorption
% This section is left blank and could be used at a future time to include (n,2n) reactions (or other absorption)
% It is commented out, but could easily be used, argin needs to include
% N_2N(r.nCells)
%
%       du=2.303*ones(p.nFineGroups,1);  'what is du?'
%
%            for j=1:1:r.nCells
%             if N_2N(j)==1,
%                for up=1:1:p.nFineGroups %(u')
%                 A3(j,:,up)=r.relVolumes(j)*du(up)*Sigma.?(j,up);
%                end;
%             elseif N_2N(j)==0,
%                 A3(j,:,:)=zeros(p.nFineGroups,p.nFineGroups);
%             else, 'ERROR bad Region encoding in region j', j
%             end,
%            end;
%         A2=size32(A3,[3,1,2],1);
%         A=I1*P1*F2;
%
%         R=R+A;

Rout=R*reshape(Phi,(i*u),1);
Fout=F*reshape(Phi,(i*u),1)*K;
Lout=L*reshape(Phi,(i*u),1);

%plot(Lout-Rout-Fout)

Rform=reshape(Rout,i,u);  %include ' if R_old is used  %Rform_old=reshape(Rout_old,u,i)';
Fform=reshape(Fout,i,u);
Lform=reshape(Lout,i,u);  %This value has been verified


%% Discussion of form (code does not run)
%
t=1; if t==0,
    %The output of the program is such that L*phi=(R+F/k)*phi
    %So to solve for phi and k use the program
    [V,D]=eig(L-R,F)
    
    DD=diag(D); list=1;
    for i=1:1:length(DD)
        if abs(DD(i))~=Inf,
            DDD(list)=DD(i);
            iii(list)=i;
            list=list+1;
        end; end;
    [DDDD,iiii]=max(DDD);
    index=iii(iiii);
    eigval=DDDD;
    eigvec=V(:,index);
end

%% size32 function
function [out]=size32(input,order,type)
%This program scales a tensor of the form
%E(i,j,k) into a tensor of the form
%
%           i=1    i=2    i=3   i=4
%       [ 1     | 1     |               ]
%       [  2    |  2    |               ]
% out=  [   3   | j=3   |               ] k=1
%       [    4  |    4  |               ]
%       [     5 |     5 |               ]
%        ----------------------------
%       [ 1     | 1     |               ] k=2
% etc

%A typical input would take the form
%[M_out]=size32( permute(M_in,[3,2,1]));
%
% Type needs to be defined properly

%Helpful functions
%    permute     - Permute array dimensions.
%    ipermute    - Inverse permute array dimensions.
%    shiftdim    - Shift dimensions.
%    reshape     - Change size


X0=permute(input,order);


[i,j,k]=size(X0);
X1=reshape(X0,i,j*k)';
if type==1,
    X2=[X1(:,1),zeros(j*k,j-1)];
    for n=2:1:i
        X2=[X2,X1(:,n),zeros(j*k,j-1)];
    end
    X3=X2;
    for n=1:1:(j*k)
        r=mod(n-1,j);
        X3(n,:)=circshift(X3(n,:),[0,r]);
    end
    out=X3;
elseif type==2,
    for n=1:1:j*k
        r=mod(n-1,j);
        %[zeros(1,i*r),X1(n,:),zeros(1,i*k-i*(r-1))];
        X2(n,:)=[zeros(1,i*r),X1(n,:),zeros(1,i*(j-r-1))];
    end
    out=X2;
    
end

