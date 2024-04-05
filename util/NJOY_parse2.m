function XSLib=NJOY_parse(filename,Bin_1ev,verbose,sparsed)
% FORMAT: XSLib=NJOY_parse(filename,Bin_1ev,verbose,sparsed);
% CALL:    >> R=NJOY_parse('output', 1,      0,      1);
%
% Program:
%   Created  : Geoff Recktenwald of UT Austin -- All rights reserved.
%   Version  : 1.00
%   Released : Jan 24nd, 2012
%
% Please send bugs and update requests to gdr7@cornell.edu
%
% This program will take NJOY output and correctly extract the appropriate 
% multigroup data.  Data will be collected into the VBUDS2 format or VBUDS2
% sparse format.
%
%  INPUTS: 
%         Filename : String  File name or FID with data
%         Bin_1ev  : Float   Thermal XS to Res XS transition [in Ev] (default 3) 
%         verbose  : Logical Set to 1 for screen output (default 0, quiet)
%         sparsed  : Logical Set to 0 for full  (default 1, sparse)
%
%  OUTPUT: 
%         XSLib.
%               edges: [111x1 double]            [ev   ] Energy group edges
%             nGroups: 110                       [#    ] number of groups      
%                  MT: [43x1 double]             [MT   ] MT reaction identifier%                 
%                 MTs: [1x452 double]            [     ] sparse map of MT
%                Temp: [8x1 double]              [K    ] temperatures
%               Temps: [1x1000 double]           [     ] sparse map of Temp
%                  s0: [7x1 double]        10^s0 [barns] self shielding values
%                 s0s: [12x1]                    [s0+2 ] shifted sparse map of s0
%                   m. [1x43 struct]
%                          xs: []                [barns] cross section of interest
%                          L1: []                        higher order terms
%                          F0: []                        higher order terms or spec
%                          F1: []                        higher order terms or prod  
%                        kern: [4-D double]      [barns] scattering kernel
%                 kern.
%                        therm: []               thermally treated kernel
%                        nothm: []               kernel with no treatment
%            
%  OUTPUT_sample:   VBUDS2 full   - R.m(R.MTs(91)).xs(:,s0,temp);
%                                                 .kern(:,:,s0,temp);
%
%                   VBUDS2 sparse - R.m(R.MTs(91)).p(s0,temp).xs(:);
%                                                            .kern(:,:);
%
%   A sample input that I've used successfully is
%
%                XSlib=NJOY_parse('output',4);
%
%  WARNINGS:  NJOY has issues with 6-91, particularly for U-235.  In that
%  case it claims to return 8-91 (but it is really 6-91).  To deal with
%  this I have sent mf8 values to mf6 parsers.  If you ever really want
%  MF8, you'll need to either modify the code or pull it out of the .kern.
%
%   NEEDED UPDATES.....
%      Allow the user to request a particular XS
%                   MFMT     : Njoy reactions requested.  [mf1 mt1; mf2 mt2; ... ]
%      Inform the user if a XS has a NAN
% 
%
%  COMMENTS:
%      Version 1.00 - Original
%      Version 1.01 - Updated to include warnings and nan
%                     To create your own warning look for EDIT EDIT
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%SET defaults
if     nargin<2, Bin_1ev=3; verbose=0; sparsed=1;
elseif nargin<3,            verbose=0; sparsed=1;
elseif nargin<4,                       sparsed=1;
end

% Allow the user to scroll through multiple outputfiles.
if strcmp(filename,'all'),
    R=dir;
    for i=1:1:length(R)
        L=fgetl(fopen(R(i).name));
        if L(1:3)=='1**', output1(i)=NJOY_parse(R(1).name,Bin_1ev); end
    end
end %This is not complete... find a way to list the materials on i.
        
    
% Make sure that the filename requested exists.    
if ~exist(filename,'file'), fprintf('Fatal Error: File %s not found.',filename); end

%Open and read the file.
FID=fopen(filename);                         % Open the file

Entire_File=fread(FID, 'uint8=>char')';        % This will fail in HUGE files.
Entire_File=strrep(Entire_File,' xx ','qqqq'); % I use q as a end flag since no Qs appear in the file (oddly enough)

%%ERROR CORRECTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(regexp(Entire_File,'---message from getunr---'))
    fprintf('WARNING: ---message from getunr--- \n');
    Entire_File=regexprep(Entire_File,'---message from getunr--- Warning, negative URR cross sections found, check unresr[\D]+',' ');
end, 
if ~isempty(regexp(Entire_File,'NaN'))
    fprintf('WARNING: Several values are undefined (nan) \n');
    Entire_File=strrep(Entire_File,'NaN+***','nan');
end

% EDIT HERE %%% These lines may be edited to remove errors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT HERE %%% This section allows you to edit for errors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EDIT EDIT %  
if ~isempty(regexp(Entire_File,'search text to indicates a problem')),  %Search for a particular string.   % EDIT EDIT %
    fprintf('WARNING: Please write error message here \n');             %Warning to be display.            % EDIT EDIT %
    % This line replaces the offending string with one that can be parsed.  Use Nan to replace             % EDIT EDIT %
    % errors in numeric values and ' ' for text errors.  Please be careful and see above for examples.     % EDIT EDIT %                                
    Entire_File=regexprep(Entire_File,'reg exp to find','reg exp replacement');  % or use strrep           % EDIT EDIT %
end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EDIT EDIT %
% END EDIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Parameters,Reactions]=Pull_initial_data(Entire_File);

%Is the reaction type a molecule?
type='free '; %0
if ~isempty(regexp(Entire_File,'for mf  6 and mt222')), type='light'; end 
if ~isempty(regexp(Entire_File,'for mf  6 and mt228')), type='heavy'; end

 B.edges        =Parameters.edges;               % Energy Bin edges [ev]
 B.nGroups      =length(Parameters.width);       % number of energy groups
 B.MTs          =sort(unique(Reactions(:,2)));   % List of Reaction types
 B.MT(B.MTs)    =sparse([1:length(B.MTs)]);      % Create MT index (sparse map)
 B.Temps        =sort(unique(Parameters.temps)); % List of Temperatures
 B.Temp(B.Temps)=sparse([1:length(B.Temps)]);    % Create Temp index (sparse map)
 B.s0s          =sort(unique(Parameters.ind_s0));% List of s0 index values
 B.s0(B.s0s)    =sparse([1:length(B.s0s)]);      %This sparse map will fail ...
                         %if s0 is not an integer.  I can deal with this by ...
                         %creating a 12*index, but I'm not sure it is worth it.
                         
% Separate all of the tallies into the array Breakup(i).text
Breakup=regexp(Entire_File,[...
                    'group constants at t=(?<temp>[\d\.\+-eE]+)'...  %Get temps
                    ' (?<what>[degks\.\d\s]+)',...                   %
                    ' for mf\s*(?<mf>\d+) and mt\s*(?<mt>\d+)',...   %Get MF and MT
                    '(?<text>[^q]+)',...                             %Get XS text
                    '            qqqq',...                           %qqqq is used to find end of XS 
                    ],'names');  
clear Entire_File;  %no longer needed (and may be large)

% Extract the data from each of the tallys
for i=1:1:length(Breakup)
    try
    [B]=extract(Breakup,B,i,verbose,sparsed);
    if ~isempty(regexp(Breakup(i).text,'nan')),   fprintf('WARNING: NANs are in your result \n'); end
    if ~isempty(regexp(Breakup(i).text,' \-\d+')),fprintf('WARNING: negative values are in your result \n'); end
    Breakup(i).text=''; %Clear the corresponding text
    catch
        %i,   %dbstop('141')
        fprintf('WARNING: Failed to extract MF=%s MT=%s and Temp=%s (see EDIT EDIT) \n',Breakup(i).mf,Breakup(i).mt,Breakup(i).temp);
    end
end


%% Post Process the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 XSLib=B;
 return
 
 %% NOTHING BELOW USED

% Create tru scatter [THERM, res, fast]
if exist('B.MT(221)'),    
   B=Build_scattering_kernel(B,sparsed,Bin_1ev)
end


% Create 6_4
try q=isempty( B.m(B.MT(4)).kern),  %no 6 4 exists
catch me, q=1; end
 if q && ~sparsed  
   B.m(B.MT(4)).kern=0;
    if B.MT(51)
      for i=51:91
          try
          B.m(B.MT(4)).kern=B.m(B.MT(4)).kern+B.m(B.MT(i)).kern;
          catch me
          %failure means that the MT doesn't exist... so ignore
          end
      end
    end
 elseif q
    SZ=size(B.m(B.MT(51)).p);
    for i=1:1:SZ(1)
     for j=1:1:SZ(2)
        B.m(B.MT(4)).p(i,j).kern=0;
        if B.MT(51)
          for ii=51:91
              try
              B.m(B.MT(4)).kern=B.m(B.MT(4)).p(i,j).kern+B.m(B.MT(ii)).p(i,j).kern;
              catch me
              %failure means that the MT doesn't exist... so ignore
              end
          end
        end
         end
    end
 end 
 
 %compare vectors with Kernel Data...
 % presumes same issues.
if verbose && ~sparsed
     for j=1:1:length(B.MTs)
         i=B.MTs(j);
        if ~isempty(B.m(B.MT(i)).kern)
            SQk=squeeze(sum(B.m(B.MT(i)).kern(:,:,:,1,1),1));
            SQt=B.m(B.MT(i)).xs(:,:,1);
            hold off,
            semilogy(SQk,'-*'); hold on
            semilogy(SQt,'-o'); 
            title(['MT = ' int2str(i)]); 
        end
     end
 end
 
 XSLib=B;
 %create total....
 
 
 
function [Parameters,Reactions]=Pull_initial_data(Entire_File)
% This function reads through Entire_file and pulls a set of parameters.
%
% Parameters.
%            temps  = [  296,  600,  900, ... 1500]    [ K] Temperatures
%            s0     = [10^10, 10^5, 10^3, ...  0.1]    [  ] Self Shielding Values
%            logs0  =     10,    5,    3, ...   -1]    [  ] Log(s0)
%            ind_s0 = [   12,    7,    5, ...    1]    [  ] Log(s0)+2   [to avoid negative indicies]
%            bin_num= [    1,    2,    3, ...  110]    [ #] bin index
%            edge_l = [10e-5,13e-5,16e-5, ... 79e5]    [ev] lhs bin values 
%            edge_r = [13e-5,16e-5,20e-5, ...  1e7]    [ev] rhs bin values
%            edges  = [10e-5,13e-5,16e-5, ...  1e7]    [ev] all bin edge values
%            width  = [26e-6,33e-6,41e-6, ...  2e6]    [ #] the number of bins
%            totbin = 110
%
% Reactions = [3 1 ; 
%              3 2 ;
%              3 4 ;
%              6 2 ;
%              6 18;
%              ... ];
%
    %Pull Parameter Data from groupr 
    Param = regexp(Entire_File,['\W temperatures \(kelvin\) \.+',...
                                          '(?<temps>[\d\.\+-eE ^\x00-\x1F]+)'...
                                 ' sigma zeroes \.+\s+infinity', ...
                                          '(?<S0>[\d\.\+-eE ^\x00-\x1F]+)'...
                                 ' neutron group structure\.+read in'...
                                          '(?<bins>[\d\.\+-eE ^\x00-\x1F]+)'...
                                ' weight function' ...
                                ],'names');
    Parameters.temps = str2num(Param(1).temps);
    Parameters.s0    = [10e10;str2num(Param(1).S0)];
    Parameters.logs0 = log10(Parameters.s0);
    Parameters.ind_s0= Parameters.logs0+2;
                A= regexprep(Param(1).bins,' - ',' ');
                B= str2num(A);
    Parameters.bin_num= B(:,1);
    Parameters.edge_l = B(:,2);
    Parameters.edge_r = B(:,3);
    Parameters.edges  = [B(:,2);B(end,3)];                  clear B;
    Parameters.width  = diff(Parameters.edges);
    Parameters.totbin = length(Parameters.width);

    MFMT_strs = regexp(Entire_File,'for mf\s*(?<mf>\d+) and mt\s*(?<mt>\d+)','names');
    for i=1:1:length(MFMT_strs)
        MFMT_nums(i,:)=[str2num(MFMT_strs(i).mf), str2num(MFMT_strs(i).mt)];
    end
    Reactions=unique(MFMT_nums,'rows');

function [B]=extract(Breakup,B,i,verbose,sparsed)
% This function takes a file segment with cross section data and
% 1. Determines the mf mt temp s0 and legendre 
%    s0 and legendre are flags (0 for only one value) 
% 2. Set the MFMT switch for the right format parsing
% 3. Parse the file
% 4. Render the output
% 5. Add output to the appropriate data structure of B
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    current.mf = str2double(Breakup(i).mf);
    current.mt = str2double(Breakup(i).mt);
    current.tmp= str2num(Breakup(i).temp);
    current.s0 = ~isempty(strfind(Breakup(i).text,' vs sigma zero')); 
    current.leg= ~isempty(strfind(Breakup(i).text,' vs legendre order')); 
    
    index.mt = find(B.MTs==current.mt);
    index.tmp= find(B.Temps==current.tmp);
    
    if verbose==1
    fprintf('\n Currently on MFMT = %d %d and Temp=%d\n',current.mf,current.mt,current.tmp);
    end
    
    %define type  (or maybe better done with a try catch statement)
    if     current.mf == 3 && current.mt == 1, MFMT = 'total'; %3_1
    elseif current.mf == 3 && current.s0 == 0, MFMT = 'vect1'; %3_451
    elseif current.mf == 3 && current.s0 == 1, MFMT = 'vectn'; %3_2
    elseif current.mf == 6 && current.mt ==18, MFMT = 'fissn'; %6_18
    elseif current.mf == 6 && current.mt ==19, MFMT = 'fissn'; %6_18
    %6_18 may also have a multi s0 and multi Leg feature feature.
    elseif current.mf == 6 && current.s0 == 1, MFMT = 'matxn'; %6_2
    elseif current.mf == 6 && current.leg== 1, MFMT = 'matxL'; %6_51 legendre  
    elseif current.mf == 6 && current.s0 == 0, MFMT = 'matx1'; %6_51   
    elseif current.mf == 8 && current.s0 == 1, MFMT = 'matxn'; %8_91 is 6_91
    elseif current.mf == 8 && current.leg== 1, MFMT = 'matxL'; %  
    elseif current.mf == 8 && current.s0 == 0, MFMT = 'matx1'; %   
    else                                       MFMT = 'other';    
    end
    
   %Lines=regexp(Breakup(i).text,'[\r\n\f]+','split'); %\n\f are included 
    Locs=regexp(Breakup(i).text,'[\r\n\f]+');          %for windows output.
    
   %Prep data
     X=Breakup(i).text(Locs(5):end);
     X=regexprep(regexprep(regexprep(X,'+','E+'),'-','E-'),' E-',' -');
     X=regexprep(X,'flux','-1'); 
     X=regexprep(X,'spec','-1');
     X=regexprep(X,'prod','-1');
    
     fprintf('MF=%d  MT=%d \tMFMT=%s \n', current.mf, current.mt, MFMT);
     
switch MFMT %##############################################################
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    case {'total','3_1'} %infinite dilution only
        M=str2num(X);
        
        if ~sparsed
         B.m(index.mt).xs(:,index.tmp,:)= M(1:4:end,3:end);
         B.m(index.mt).L1(:,index.tmp,:)= M(2:4:end,3:end);
         B.m(index.mt).F0(:,index.tmp,:)= M(3:4:end,3:end);
         B.m(index.mt).F1(:,index.tmp,:)= M(4:4:end,3:end);
        else
         for s1=1:length(B.s0s)
          B.m(index.mt).p(index.tmp,s1).xs = M(1:4:end,s1+2);
          B.m(index.mt).p(index.tmp,s1).L1 = M(2:4:end,s1+2);
          B.m(index.mt).p(index.tmp,s1).F0 = M(3:4:end,s1+2);
          B.m(index.mt).p(index.tmp,s1).F1 = M(4:4:end,s1+2);
         end     
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case {'vect1','3_251'} %infinite dilution only
        M=str2num(X);
        if ~sparsed
         B.m(index.mt).xs(1:B.nGroups,index.tmp)= 0;
         B.m(index.mt).xs(M(:,1),index.tmp)= M(:,2);
        else  
          s1=1;  
          B.m(index.mt).p(index.tmp,s1).xs(1:B.nGroups)= 0;
          B.m(index.mt).p(index.tmp,s1).xs(M(:,1)) = M(:,2);
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
    case {'vectn','3_2'} %Given as a function of S0     
        M=str2num(X);         
        if ~sparsed 
            B.m(index.mt).xs(1:B.nGroups,index.tmp,1:length(B.s0s)) = 0;
            B.m(index.mt).xs(M(:,1),index.tmp,:)= M(:,2:end);
        else
         for s1=1:length(B.s0s)
            B.m(index.mt).p(index.tmp,s1).xs(1:B.nGroups)= 0;  
            B.m(index.mt).p(index.tmp,s1).xs(M(:,1))= M(:,s1+1);
         end     
        end
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    case {'fissn','6_18'} %Find the matrix infinite dilution only        
        if current.leg~=0, 'WARNING: Legendre Coeffs', end %I've not seen a version that has Legendre
        Lines=regexp(X,'[\r\n\f]+ ?[\r\n\f]+','split');    %Split into 3 distinct sections for parsing
        if length(Lines{1})<20, q=1; end                   %The first line should be '';
          
   %%%%%%FROM 6_51 - get spec
         M=textscan(Lines{1+q},'%f %f %f %f %f %f %f %f','CollectOutput',1);
         M=M{1};
        
         spec=zeros(B.nGroups,1);  %initialize
        
         [Dy,Dx]=size(M);
          for j1=1:1:Dy 
            for j=0:1:Dx-3
              if ~isnan(M(j1,3+j))   
                 spec(M(j1,2)+j)=M(j1,3+j); %2nd col is location, 3:end is data
              end
            end
          end
        if ~sparsed 
           B.m(index.mt).F0(:,index.tmp)=spec;
        else
           B.m(index.mt).p(index.tmp,1).F0=spec;
        end
        clear M j1 j Dy Dx
        
   %%%%%FROM 6_51 - get prod
        M=textscan(Lines{2+q},'%f %f %f %f %f %f %f %f','CollectOutput',1);
        M=M{1};
        
        prod=zeros(B.nGroups,1); %initialize
        
        [Dy,Dx]=size(M);
          for j1=1:1:Dy
            for j=0:1:Dx-3
              if ~isnan(M(j1,3+j))   
                 prod(M(j1,1)+j)=M(j1,3+j);  %1st col is location, 3:end is data
              end
            end
          end
         
        if ~sparsed 
           B.m(index.mt).F1(:,index.tmp)=prod;
        else
           B.m(index.mt).p(index.tmp,1).F1=prod;
        end
        clear M j1 j Dy Dx
        
  %%%%%%FROM 6_2 - parse matrix
        M=str2num(Lines{3+q});
        Leg=max(M(:,3))>0; %More than one Legendre coeff
          
         if ~sparsed             
              B.m(index.mt).kern(1:B.nGroups,1:B.nGroups,index.tmp,1:length(B.s0s)) = 0;
         else
             for s1=1:1:length(B.s0s)
                 B.m(index.mt).p(index.tmp,s1).kern(1:B.nGroups,1:B.nGroups) = sparse(0);
             end
         end
          
        if Leg
           for j=1:1:length(M) 
             if ~sparsed 
                B.m(index.mt).kern(M(j,2),M(j,1),index.tmp,:,M(j,3)+1)= M(j,4:end);
             else
                for s1=1:length(B.s0s)
%BAD               B.m(index.mt).p(index.tmp,s1,M(j,3)+1).kern(1:B.nGroups,1:B.nGroups) = sparse(0);
                   B.m(index.mt).p(index.tmp,s1,M(j,3)+1).kern(M(j,2),M(j,1))=M(j,s1+3);
                end
             end
           end
        else
           for j=1:1:length(M) 
                if ~sparsed 
                B.m(index.mt).kern(M(j,2),M(j,1),index.tmp,:)= M(j,4:end);
                else
                    for s1=1:length(B.s0s) 
%BAD                   B.m(index.mt).p(index.tmp,s1).kern(1:B.nGroups,1:B.nGroups) = sparse(0);
                       B.m(index.mt).p(index.tmp,s1).kern(M(j,2),M(j,1))=M(j,s1+3);
                    end
                end          %  Output format will look like
           end               %  Z( i -> j);
        end                  %  CrossSection=sum(Z,1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    case {'matx1','6_51'} %Find the matrix infinite dilution only     
        %M=textscan([X '0 0 0 0 0 0 0 0'],'%f %f %f %f %f %f %f %f');
        %if length(M)>8, 'FATAL ERROR', end
        %M=[M{1},M{2},M{3},M{4},M{5},M{6},M{7},M{8}];  
        %M=M(1:end-1,:);
        
        M=textscan(X,'%f %f %f %f %f %f %f %f','CollectOutput',1);
        M=M{1};
       
        Z=zeros(B.nGroups);
        
        [Dy,Dx]=size(M);
          for j1=1:1:Dy
            for j=0:1:Dx-3
              if ~isnan(M(j1,3+j))   
                 Z(M(j1,2)+j,M(j1,1))=M(j1,3+j); 
              end
            end
          end
       
         if ~sparsed 
             B.m(index.mt).kern(:,:,index.tmp)=Z; 
         else 
             B.m(index.mt).p(index.tmp,1).kern=sparse(Z); 
         end   
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
    case {'matxL','6_51leg'} %Find the matrix Legendre is handled but not S0    
        M=str2num(X);
    
        if ~sparsed        
           B.m(index.mt).kern(1:B.nGroups,1:B.nGroups,index.tmp,1,1) = 0;
          for j=1:1:length(M) 
           B.m(index.mt).kern(M(j,2),M(j,1),index.tmp,1,1:length(M(j,3:end)))= M(j,3:end);
          end
        else
          for j=1:1:length(M) 
              B.m(index.mt).p(index.tmp,1,1:length(M(j,3:end))).kern(1:B.nGroups,1:B.nGroups) = sparse(0);
            for s1=1:1:length(M(j,3:end))     
              B.m(index.mt).p(index.tmp,1,s1).kern(M(j,2),M(j,1)) = M(j,s1+2);
            end          %  Output format will look like
          end               %  Z( i -> j);
        end                  %  CrossSection=sum(Z,1);  
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    case {'matxn','6_2'} %Find the matrix S0 and Legendre handled    
        M=str2num(X);
        Leg=max(M(:,3))>0; %More than one Legendre coeff
        
         if ~sparsed             
              B.m(index.mt).kern(1:B.nGroups,1:B.nGroups,index.tmp,1:length(B.s0s)) = 0;
         else
             for s1=1:1:length(B.s0s)
                 B.m(index.mt).p(index.tmp,s1).kern(1:B.nGroups,1:B.nGroups) = sparse(0);
             end
         end
         
        if Leg
           for j=1:1:length(M) 
             if ~sparsed 
                B.m(index.mt).kern(M(j,2),M(j,1),index.tmp,:,M(j,3)+1)= M(j,4:end);
             else
                for s1=1:length(B.s0s)
 %BAD                    B.m(index.mt).p(index.tmp,s1,M(j,3)+1).kern(1:B.nGroups,1:B.nGroups) = sparse(0);
                   B.m(index.mt).p(index.tmp,:,M(j,3)+1).kern(M(j,2),M(j,1))=M(j,s1+3);
                end
             end
           end
        else
           for j=1:1:length(M) 
                if ~sparsed 
                B.m(index.mt).kern(M(j,2),M(j,1),index.tmp,:)= M(j,4:end);
                else
                    for s1=1:length(B.s0s) 
                       B.m(index.mt).p(index.tmp,s1).kern(M(j,2),M(j,1))=M(j,s1+3);
                    end
                end            %  Z( i -> j);
           end                 %  CrossSection=sum(Z,1);   
        end
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    case {'other'}
        fprintf('\nThis code cannot handle MFMT = %d %d\n',current.mf,current.mt);
    otherwise
        fprintf('\nError: This code cannot handle MFMT = %d %d\n',current.mf,current.mt);
end%#######################################################################        
        
%% CREATE KERNELS
function  [B]=Build_scattering_kernel(B,sparsed,Bin_1ev)

 if ~sparsed
   B.kern.no_t=B.m(B.MT(2)).kern;   %save the 2 in a new location
    
   SZ2=size(B.kern.no_t);   
   SZ221 = size(B.m(B.MT(221)).kern);
   
   if length(SZ221)<4, i2(1:SZ2(4))=1; else i2=[1:SZ2]; end %Deal with extra S0 Values
 
   for i=1:1:SZ2(4)   %combine them with the combine function
       for j=1:1:SZ2(3)
            Tempvar(:,:,j,i)=...
                   combineK(...
                            squeeze(B.kern.no_t(:,:,j,i,1)),...
                            squeeze(B.m(B.MT(221)).kern(:,:,j,i2(i))),...
                            Bin_1ev  ); 
       end
   end
   B.m(B.MT(2)).kern = Tempvar;   
 
 else
     B.kern.no_t=B.m(B.MT(2)).p;   %save the 2 in a new location
     SZ2=size(B.kern.no_t);   

     SZ221 = size(B.m(B.MT(221)).p);
   
     if SZ221(2)==1, i2(1:SZ2(2))=1; else i2=[1:SZ2]; end %Deal with extra S0 Values

     
     for j=1:1:SZ2(1) %along temps
      for i=1:1:SZ2(2) % along s0
           Tempvar=combineK(...
                            squeeze(B.kern.no_t(j,i).kern(:,:,1)),...
                            squeeze(B.m(B.MT(221)).p(j,i2(i)).kern(:,:,1)),...
                            Bin_1ev  ); 
           B.m(B.MT(2)).p(j,i).kern = Tempvar; 
      end
     end
     
end    
     
   
%% Combine functions
function [output,error]=combineM(X,Y,Bin_1ev)
         [output,error]=combineK(X,Y,Bin_1ev);  %Matrix and Kernal are the same

function [output,error]=combineK(X,Y,Bin_1ev)
%This function combines the matricies X and Y
%The thermal part of X (below Bin_1ev) is replaced by the region in Y
%  Typically X=6_2 and Y=with 6_221, 6_222, or 6_228.  
%It effectively subs 6_22x in for the thermal part of 6_2.
%
%The full length of Y is maintained and an error structure is available for
%later analysis.

   a=Bin_1ev;
   b=length(Y);
   x=X;
   X(1:b,1:a)=Y(1:b,1:a);
   output=X;
   [A,B]=size(Y);
   error.matrix=abs(x(1:A,1:B)-Y);
                xv=sum(x,1); Yv=sum(Y,1);
   error.vector=abs(xv(1:length(Yv))-Yv); %vector minima is the best bin.
   [error.value,error.index] = min(error.vector); 
                   
function [output,error]=combineV(X,Y,Bin_1ev)
%This function combines vectors X and Y,
%The thermal part of X (below Bin_1ev) is replaced by the region in Y
%Error Checking is done to ensure Bin_1ev is in the 'correct' location.
   a=Bin_1ev;
   x=X;
   X(1:a)=Y(1:a);
   output=X;
            
   error.vector= abs(x(1:length(Y))-Y);  %vector minima is the best bin.
   [error.value,error.index] = min(error.vector); 
   
   if error.index~=Bin_1ev, 
       fprintf('New bin_1ev suggested: %d, old: %d\n',error.index,Bin_1ev); 
   end
