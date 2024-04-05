% [t_b,t_c,t_k,libs]=MCNPX_temp(temps,flags,water); flags='c' 'k' 'mev' 'mev2k' etc...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% File Created by Geoffrey Recktenwald
%                 Cornell University 
%                 March 23rd, 2011      
%                 V1.10
%
% Description:  This file will take a vector of temperatures and return a 
% vector of temperatures in different units.  
%
% USE: [t_b,t_c,t_k,lib]=MCNPX_temp(temps,flags)
%
% Input: 
%         temps = a vector of temperatures [must be in same units] (give units in flags)
%         flags = the input temperature and desired output temperature (or all)
%         water = put something here to look for libs for ltwr.
%
%            flags are used to specify the input type.
%               flag='c'     Celcius to all three
%                   ='k'     Kelvin  to all three
%                   ='mev'   Mev     to all three
%                   ='c2mev' Celcius to Mev
%                   ='c2k'   Celcius to Kelvin
%                   ='mev2c' Kelvin  to Celcius
%                   ='k2mev' Kelvin  to Mev
%                   ='mev2k' Mev     to Kelvin
%                   ='mev2c' Mev     to Celcius
%                   ='c2l'   Celcius to library string
%                   ='k2l'   Kelvin  to library string
%                   ='mev2l' Mev     to library string
%                   ='c2ln'  Celcius to library decimal
%                   ='k2ln'  Kelvin  to library decimal
%                   ='mev2ln'Mev     to library decimal
%     Future        ='cNOl'  Celcius skip library <= only to be written if
%                                                    speed is an issue.
% Output:
%          t_b = boltzmann 'temperature'   <== mcnp input temps
%          t_c = temperature in Celcius
%          t_k = temperature in Kelvin
%          libs= MCNP library designator string.  The MCNPX library chosen 
%                will always be the maximum library temp such that 
%                librarytemp < temps
%          libn= MCNP library decimal (0.72).  
%  
% EXCEPTION: Some flags request an output.  In those cases, the first 
%            argout, t_b, will be the temperature in the requested units. 
%
% 
% MCNP Library options        T < 293.6 - 'Unkn' ==> Unknown
%                      293.6< T < 600   - '.70c' 
%                      600  < T < 900   - '.71c'
%                      900  < T < 1200  - '.72c' 
%                     1200  < T < 2500  - '.73c'
%                     2500  < T         - '.74c'
% EXCEPTION:  I use '.70c' for any temp T<600 but issue as warning if T < 292
%
%
% NOTE: if no flag is specified, the code will choose Celcius or Mev
% depending on the magnitude of the temps.
%
% Code: note I use mev and b to denote the same thing (b for Boltzmann)
%
% Examples:
%         [t_b,t_c,t_k]=MCNPX_temp([300, 900],'c');
%         [~,~,t_k]=MCNPX_temp([300, 900],'c');
%         [t_b]=MCNPX_temp([300, 900],'c2mev');
%         [t_k]=MCNPX_temp([4.94e-08 1.011e-07],'mev2k');
%
% Note:   The following are equivalent:
%                          [~,~,t_k]=MCNPX_temp([300, 900],'c');
%                               t_k =MCNPX_temp([300, 900],'c2k');
%
% Versions: 1.00 - Initial version
%           1.10 - Reformatted code and methods 
%                - Added mcnp library codes
%                - Added mcnp library codes for water
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t_b,t_c,t_k,libs]=MCNPX_temp(temps,flags,lwtr)

if ~exist('flags','var'), 
    if max(temps) > 1, flags='c'; 
    else flags='mev'; end
end

c2k=273.15;        % 0 c in kelvin
k2b=8.6173423E-11; % 1 kelvin in Mev (1 k * (Boltzmann constant) / 1 Mev)

switch lower(flags(1))
    case 'c'
       t_c=temps;
       t_k=temps+c2k;
       t_b=t_k*k2b;
    case 'k'
       t_c=temps-c2k;
       t_k=temps;
       t_b=t_k*k2b;
    case {'m','b'}
       t_b=temps;
       t_k=t_b/k2b;
       t_c=t_k-c2k; 
    otherwise
        fprintf('\n\n Please choose a valid flag \n\n')
        t_b=-1; t_c=-1; t_k=-1; lib='Unkn';
end  

%build library strings
if nargin==3  %(use water libraries)
 L.temps=[ 294, 350:50:650, 800,inf]-2;
 L.decim=[0.1:0.01:0.18];
 L.strng={'.10t','.11t','.12t','.13t','.14t','.15t','.16t','.17t','.18t'};
 %fprintf('\n\n Water uses a different scheme  lwtr.14t \n\n');
else
 L.temps=[ 296, 600, 900,1200,2500,inf]-5;
 L.decim=[0.70,0.71,0.72,0.73,0.74];
 L.strng={'.70c','.71c','.72c','.73c','.74c'};
end

if length(t_k)==1,
    i=find(L.temps>t_k,1)-1; i=max(1,i);
    libn=L.decim(i); 
    %fprintf('\n strings end in c.  ex: .70c \n ')
    libs=L.strng{i};
else
   for ii=1:length(t_k),
      i=find(L.temps>t_k(ii),1)-1; i=max(1,i);
      libn(ii)=L.decim(i);
      libs{ii}=L.strng{i};
   end
end

%deal with single output requests.
switch lower(flags)
    case {'c','k','b','mev'}
       %nothing needs to be done
    case {'c2mev','c2b','k2mev','k2b'},
        %nothing needs to be done
    case {'c2k','mev2k','b2k'}
        t_b=t_k;
    case {'k2c','mev2c','b2c'}
        t_b=t_c;
    case {'c2l','k2l','mev2l','c2ls','k2ls','mev2ls'}
        t_b=libs;
    case {'c2ln','k2ln','mev2ln','c2ld','k2ld','mev2ld'}
        t_b=libn;
    otherwise   
         fprintf('\n\n Please choose a valid flag \n\n')
end 





%% IGNORE EVERYTHING BELOW THIS LINE  (legacy code)

% FASTER version if wanted but no libs
%
% switch lower(flags)
%     case 'c'
%        t_c=temps;
%        t_k=temps+c2k;
%        t_b=t_k*k2b;
%     case 'k'
%        t_c=temps-c2k;
%        t_k=temps;
%        t_b=t_k*k2b;
%     case {'mev','b'}
%        t_b=temps;
%        t_k=t_b/k2b;
%        t_c=t_k-c2k;
%     case {'c2mev','c2b'}, t_b=(temps+c2k)*k2b;
%     case 'c2k',           t_b=(temps+c2k)    ;
%     case 'k2c',           t_b=(temps-c2k)    ;
%     case {'k2mev','k2b'}, t_b=(temps    )*k2b;
%     case {'mev2k','b2k'}, t_b=(temps/k2b)    ;
%     case {'mev2c','b2c'}, t_b=(temps/k2b)-c2k;
%     case {'c2l','k2l','mev2l'}
%         t_b=lib;
%     
%     otherwise   
%          fprintf('\n\n Please choose a valid flag \n\n')
%          t_b=-1; t_c=-1; t_k=-1;
% end
