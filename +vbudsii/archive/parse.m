function output1=parse(filename,mat,Bin_1ev,bound)
%This program will parse a NJOY output file.
% filename - is the name of the file
% mat - is the material name (eg U238)
% Bin_1ev is the number of the bin whose LHS is 1 EV
% bound is a random flag to let the program know the data was bound.
%
% Sample: output1=parse('output','U238',40)
%
%
    %filename='output';
    %mat='U238';    
    fid = fopen(filename);
    if nargin==2, 'STOP check 1ev bin', Bin_1ev=4; end
    if nargin==4, bound='yes';, end
    
while 1
         tline = fgetl(fid);
          if ~ischar(tline)
          break 
          end
        
    
    if length(tline)>=37 %30
        
        if strcmp(tline(1:37),' neutron group structure......read in'),  %Get the bin edges
            [bins,edges]=getbins(fid);    
        elseif tline(21:22)=='t='  %Set Temperature
            Temp=str2num(tline(23:31));  TEMP=int2str(Temp);
        elseif tline(2:4)=='for'  %Set MFMT and get data
            MF=tline(10);
            MT=int2str(str2num(tline(18:20)));
            

            %saves data in the format U238_3_1_600
            if bound=='yes' & ( str2num([MF MT])==62 | str2num([MF MT])==31 ),
                eval([mat '_' MF '_' MT '_' TEMP '=MFMT_' MF '_' MT '_bound(fid);']);
            else
                eval([mat '_' MF '_' MT '_' TEMP '=MFMT_' MF '_' MT '(fid);']);
            end
            
            %This line calls COMPARE to correct the thermal region scatter kernals.
            if strcmp(MF,'6') & (strcmp(MT,'221') | strcmp(MT,'222'))
                eval([ mat '_0_2_' TEMP '=combine(' mat '_6_2_' TEMP ',' mat '_6_' MT '_' TEMP ',' 'Bin_1ev' ');']);
            end;
     end;end;
end %while loop
    
%Save the data to a .mat file    
 eval(['save Prop' mat]);
    
 fclose(fid);
 
 output1=mat;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
%% Bin / edge function
function [left,edges]=getbins(fid)
     k=0; tline = fgetl(fid);
     while tline(10)~=' '
       k=k+1;
       left(k)=str2double(tline(10:20)); % str2num(tline(10:20));
       right(k)=str2double(tline(26:36));
       tline = fgetl(fid);
       if length(tline)<=10, break, end
     end
     edges=[left(1) right];
     
%% MFMT parsing functions     
function [output]=MFMT_3_1(fid) %this function is for the not-bound kernels.
          Data=ripdata(fid);
      
          
          k=0;
        for i=1:1:length(Data)
            if Data(i,1)~=-10
                if Data(i,2)==0
                    k=k+1;
                    X(k,:)=Data(i,:);
                else
                    Y(k,:)=Data(i,:);
                end  
            else
              %I don't think this section is needed
                if Data(i,2)==0
                        XX(k,:)=Data(i,:);
                else
                        YY(k,:)=Data(i,:);
                end
              %end of not needed section
            end
        end
        
        output=X(:,3);

function [output]=MFMT_3_1_bound(fid) %this function is for the not-bound kernels.      
          Data=ripdata(fid);
          k=0;
        for i=1:1:length(Data)
            if Data(i,1)~=-10
               k=k+1;
               X(k,:)=Data(i,:);
            end
        end
        
        output=X(:,2);  %this is should be  3 in the unbound version....


function [output]=MFMT_3_2(fid)
          Data=ripdata(fid);
          output=Data(:,2);
          
function [output]=MFMT_6_2(fid)
          Data=ripdata(fid,1);
 
          j=max([Data(:,1);Data(:,2)]);
          Z=zeros(j,j);
          if Data(1,3)~=0, 'ERROR - should use Data(i,3)', end
          for i=1:1:length(Data)
            Z(Data(i,2),Data(i,1))=Data(i,4);  %may be 3 depending on NJOY script
          end
            %  Z( i -> j)
          %CrossSection=sum(Z,1);
          output=Z;
          
function [output]=MFMT_6_2_bound(fid) %for 6.2 for the bound, data is
          output =MFMT_6_221(fid);    %processed like 6.221

          
function [output]=MFMT_3_4(fid)
          Data=ripdata(fid);
          output=[zeros(Data(1,1)-1,1); Data(:,2)];
          %check to be sure you include enough zeros on the end.
          
function [output]=MFMT_3_18(fid)
          Data=ripdata(fid);
          output=Data(:,2);
       
function [output]=MFMT_6_18(fid)
          'skip MFMT 6_18 at this time'
          output='none';          
          
function [output]=MFMT_3_102(fid)
          Data=ripdata(fid);
          output=Data(:,2);
  
function [output]=MFMT_3_221(fid)
          Data=ripdata(fid);
          output=Data(:,2);
    
function [output]=MFMT_6_221(fid)
          Data=ripdata(fid,1);
          
          [Dy,Dx]=size(Data);
          for i=1:1:Dy
            for j=0:1:Dx-3
              if Data(i,3+j)~=-20  
                 Z(Data(i,2)+j,Data(i,1))=Data(i,3+j); 
              end
            end
          end
         output=Z;

          %output=Data(:,3:end)';
          
function [output]=MFMT_3_222(fid)
          Data=ripdata(fid);
          output=Data(:,2);
       
function [output]=MFMT_6_222(fid)
          output =MFMT_6_221(fid);    %processed like 6.221
          %Data=ripdata(fid,1);
          %output=Data(:,3:end)';   %there is a ' here to take care of orientation   
       
function [output]=MFMT_3_251(fid)
          Data=ripdata(fid);
          output=Data(:,2);
  
function [output]=MFMT_3_452(fid)
          Data=ripdata(fid);
          output=Data(:,2);

%% Combine function
function [output]=combine(X,Y,Bin_1ev)
%This function combines the 6_2 (X) with 6_221 or 6_222 (Y).  It effectively subs
%6_22x in for the thermal part of 6_2.
   a=Bin_1ev;
   b=length(Y);
   x=X;   %one needs to be chose carefully
   X(1:b,1:a)=Y(1:b,1:a);
   output=X;
            combine_accuracy=  [sum(X,1); sum(x,1); sum(X,1)-sum(x,1)];
                          %  difference at cutoff /  difference at cutoff+1  /  max difference before cutoff  /  max difference after cutoff  
            combine_values=[combine_accuracy(3,a), sum(Y(:,a+1))-sum(X(:,a+1)),max(abs(combine_accuracy(3,:))) , sum(Y(:,end))-sum(X(:,b)) ] 
   
%% Data ripping function
function [data]=ripdata(fid,varargin)
%This function will skip 4 lines then start collecting data.
%It automatically replaces '+' with 'E+' and '-' with 'E-'.
%It returns the data as an array.  Special care is taken to ensure that
%MF6 MT221 is properly recorded (there are size issues)

       fgetl(fid);fgetl(fid);fgetl(fid);fgetl(fid);
       %ftell(fid)

       k=0; 
       flux=-10; flx=-10;  %NJOY has a term called flux, I'm replacing with -10 for ease of data processing.
       tline = fgetl(fid); 

      while 1
            if length(tline)<11 | tline(11)==' ',  break, end
            k=k+1;
                 if nargin==1 | k==1,
                     
                     
                   eval([    'data(k,:)= [' regexprep(regexprep(regexprep(tline,'+','E+'),'-','E-'),' E-',' -') '];' ])
                   % eval([    'data(k,:)= [' regexprep(regexprep(tline,'+','E+'),'-','E-') '];' ])
                 elseif nargin==2, %check for mismatched data (IE things like MFMT_6_221)
                   eval([    'DD= [' regexprep(regexprep(tline,'+','E+'),'-','E-') '];' ])  
                   mismatch=length(DD)-size(data,2); %this line deals with mf6_221
                     if mismatch > 0
                       data=[data,zeros(k-1,mismatch)];  %this only works if 6_221 is longer
                       data(k,:)=DD;
                     elseif mismatch < 0
                       data(k,:)=[DD, -20*ones(1,-mismatch)];
                     else
                      data(k,:)=DD; 
                     end
                 end
            tline = fgetl(fid); 
      end

