function [done]=findreplace(FID, ...
                            strings2find, ...
                            strings2replace, ...
                            flagreplacements, ...
                            varargin)
%
%
%    [done]=findreplace(FID,{'test'},{'return'},'@')
%
%    String finds may be given in regexp, but replacements may not.
%    eventually we can make this smart.
%
%    At this point, please don't search for @'s
%
%    This program was written for a full file search,  It can do a single
%    string, but it would be SLOW to use that way often.
%
%  MAJOR ISSUES IF % appears why?

if size(varargin) >= 1
    file_out_name = varargin{1};
else
    file_out_name = 'temporaryfile';
end

FID2=fopen(file_out_name,'w');

number = length(strings2find);

%warnings
if length(strings2find)~=length(strings2replace) 
    fprintf('\n ERROR: Number of strings to find must equal the number of strings to replace \n')
end


for i=1:1:number
    for j=2:1:number
        if j>i
           if ~isempty(strfind( strings2find{i} , strings2replace{j})); 
              fprintf('\n WARNING: You will be replacing a replacement if you continue \n')
              fprintf('\n          Please consider changing the order.')
           end
        end
    end
    stringlength(i)=length(strings2find{i});  %not part of the warnings
end

%Flagreplacements
if nargin >= 4 && flagreplacements
    for i=1:1:number
        strings2replace{i}=['@' strings2replace{i} '@'];
    end
end

%Line or FID or File
if isnumeric(FID)  %identifier 
    CURRENT_LOCATION=ftell(FID);
    fseek(FID,0,-1);
    L=fgets(FID);
elseif exist(FID) %file
    FID=fopen(FID); closeit='yes';
    L=fgets(FID);
else %line
    L=FID;
end

    FRs=zeros(1,number);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while L~=-1 %%%%%%%%%%%%%%%%%%% MAIN CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %  L, pause
    %find and replace in string
    L=['@' L '@'];
    for i=1:1:number
        A=strfind(L,strings2find{i})-1;
        B=A+stringlength(i)+1;
        for j=length(B):-1:1
            L=[L(1:A(j)) strings2replace{i} L(B(j):end)];  
        end
        FRs(i)=FRs(i)+length(B);
    end
    L=L(2:end-1);
    
    %write string to file
    fprintf(FID2,L);
    
    %Get next string in file
    if ~isnumeric(FID), break, end
    L=fgets(FID);
    
end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    done=FRs;
    fclose(FID2);

    %Reset the file location (WHO KNOWS IF THIS IS NECESSARY
    if exist('CURRENT_LOCATION')
        fseek(FID,CURRENT_LOCATION,-1);
    elseif isnumeric(FID), fclose(FID)
    else done=L; delete file_out_name, return
    end

    
   
    fprintf('\n\nI have found and replaced following values:\n\n');
    fprintf('\t Number of replacments \t String_found \t \t Replacement \n')
    for i=1:1:number
      fprintf('\t \t %d \t \t %s \t \t \t %s \n',FRs(i),strings2find{i},strings2replace{i})
    end
    
    fprintf('\n\n')
