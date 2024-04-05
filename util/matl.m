%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% File Created by Susan Cardinali
%                 Cornell University
%                 July 20th, 2010
%                 V1.02
%
% Description:  This file will take a material matrix and return the atomic
% percentages, weight percentages, atoms/b-cm, and density.
%
% Use:  [ao,wo,atom/bcm,density,zaid] = matl( material_vector , type, density, flags)
%
%    material_vector = [ ZAID, ###; ZAID, ###; ...];   FORMAT [2xN] matrix
%    type =  1 if ### is an atomic %                   FORMAT number
%         =  2 if ### is a weight %                    Notes: does not have 
%         =  3 if ### is atoms/b-cm                        to be normalized
%    density = # in [g/cm^3]                           FORMAT number              
%    flags = [ 1 1 ]                                   FORMAT [1x2] logical
%
% Notes: This function can handle non-normalized atomic and weight percents
%        If type = 3, the density is ignored and the calculated density is
%        used.  If you want to use your own density set type = 1;
%        flags allow you to turn off post processing.  (it is not a necessary input)
%              flags(1)=0 - turn off stripzeros 
%              flags(2)=0 - turn off combineduplicates
%
% Examples: matl( [92235 0.05; 92238 0.95; 8016 2] , 1, 11.5);
%           matl( [40000 92; 94239; 8] , 2, 6.10);
%           matl( [40000 92; 94239; 8] , 2, 6.10, [0,1]);
%           matl( [92235 0.001; 92238 0.02] , 3 ); 
%
%   An example where you have volume percents (.44 & .56 below);
%           matl( [11023 0.44*0.925; 92238 0.56*19] , 2, 0.925*0.44+19*0.56) 
%
%
%
% UPDATES:  
%         1.01 - included Zaid list in output
%              - made volume percent [v/o] example
%         1.02 - included stripzeros
%              - included combineduplicates
%              - added input 'flags' to indicate if the user would like to 
%                avoid post processing the data with the 2 new functions.
%              
% File needs to be updated to parse library type ZAIDS like 92235.35c 
% File needs to be updated to handle sparse input. (parse index as Zaid)           
%
% Dependencies: This file requires access to AW.m (AW.m is a lookup table
% for atomic weights).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ao wo X delta zaid]= matl(vecomatls, type, delta, flags)

Na=0.6022E24;        % [atoms/(g*mol)] Avogadro's number
zaid=vecomatls(:,1); % ZAID            pull zaid numbers
atwt=AW(zaid);       % [Amu]          get atomic weights
atwtg=atwt*1.66E-24; % [grams]         atomic weight in grams
Mu=1;                % [g/mol];  molar mass constant keeps units consistent

%% JJB EDIT FIX BOUND HYDROGEN
if ismember([11],vecomatls(:,1))
    fixBoundHydrogen=vecomatls(:,1);
    fixBoundHydrogen(fixBoundHydrogen==11)=1001;
    vecomatls(:,1)=fixBoundHydrogen;
    
end
%%
% Determine and flag inputs.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ((nargin==2) && (exist('type')) && (exist('vecomatls')))
    if (type==1) || (type==2)
        fprintf('\n Warning:Density [g/cm^3] not given! \n');
        fprintf('\n Warning: X cannot be calculated. \n');
    else% (type==3)
        delta=sum(atwtg.*vecomatls(:,2)*10^24);
    end
    flags=[1,1];
elseif (nargin==3)
    if type==3, 
        Delta=sum(atwtg.*vecomatls(:,2)*10^24);
        fprintf('\n Density is ignored for input type [atoms/b-cm]! \n');
        fprintf(' Input density = %0.3g, Calculated density = %0.3g \n', delta,Delta);
        fprintf(' Input density is off by %0.3g% \n ', (delta-Delta)/(Delta) );
        fprintf(' Density has been set to calculated density \n ', (delta-Delta)/(Delta) );
        delta=Delta; clear Delta
    end
    flags=[1,1];
elseif (nargin<2)
    fprintf('\n Too few input parameters. \n');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (type==1) %atom percent%%%%%%%%%%%%%%%%%%%%works%%%%%%%%%
    ao=100*vecomatls(:,2)/sum(vecomatls(:,2));

    wo=100*((ao*Na).*(atwtg)*Mu)/sum((ao*Na).*(atwtg)*Mu);

    X=ao*delta*10^(-24)*(dot(ao,atwtg)).^-1;
    
    
elseif (type==2) %weight percent  %need to fix X
    wo=vecomatls(:,2)./sum(vecomatls(:,2));

    ai=wo.*(sum(atwtg)./atwtg);
    ao=100*ai/sum(ai);
    X=ao*delta*10^-24*(dot(ao,atwtg)).^-1;
    wo=100*wo;
    
    
elseif (type==3) %actual density  %%%%%%%%%%%WORKS%%%%%%%%%%%
    %define actual density
    X=vecomatls(:,2);
    %find total number of atoms (scalar)
    a_tot=sum(X);
    %find atom percent (vector)
    ao=(100/a_tot)*(X);
    %find total weight (scalar)
    %     Xreal=X(1);
    %     X(1)=0;
    w_tot=sum(X.*atwtg);
    %find weight percent (vector)
    wo=(100/w_tot)*(X.*atwtg);
    %     X(1)=Xreal; 
else fprint('\n\n   Type must be an integer 1,2, or 3 \n\n') 
end

%post processing (default is ON)
if flags(1)==1, [ao wo X delta zaid]=stripzeros(ao, wo, X, delta, zaid);        end
if flags(2)==1, [ao wo X delta zaid]=combineduplicates(ao, wo, X, delta, zaid); end

end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ao, wo, X, delta, zaid]=stripzeros(ao, wo, X, delta, zaid)
         if sum(ao==0)~=0,
             index=(ao~=0);
             ao=ao(index);
             wo=wo(index);
             X= X(index);
             zaid=zaid(index);
             delta=delta;
         end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ao, wo, X, delta, zaid]=combineduplicates(ao, wo, X, delta, zaid)
         if length(zaid)~=length(unique(zaid)), 
              [a1,a2]=sort(zaid);
              a3=[1; diff(a1)]~=0;
              index=sort(a2(a3));

              AO=ao;
              ZAID=zaid;
              
              for i=1:1:length(a1),
                  if a3(i)~=0, j=a2(i);
                  else 
                      ao(j)=ao(j)+ao(a2(i));
                      wo(j)=wo(j)+wo(a2(i));
                      X(j)=  X(j)+ X(a2(i));
                  end
              end
              
              ao=ao(index);
              wo=wo(index);
              X= X(index);
              zaid=zaid(index);
              delta=delta;
              
              fprintf('\n \n You have inserted a material(s) more than once. \n \n    The vector \n \n')
              fprintf('       %d \t %2.3f \n',[ZAID,AO]')
              fprintf('    Has been replaced by the vector: \n \n')
              fprintf('       %d \t %2.3f \n',[zaid,ao]')
              
              
         end
end
