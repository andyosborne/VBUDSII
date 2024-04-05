function [output]=MCNPX_clean(filename,option)
%
%  This function will open an MCNPX file and go through line by line to
%  make sure there are no disallowed characters in the MCNP input deck.
%
%  Option=0 or blank just returns a warning that there are bad characters
%  and on what lines.
%  Option=1 deletes the offending characters (use with caution!)
%
%  Needed
%     Remove terms and print to file if option is 1
%     Add option 2 where a specific character can be sought.
%     Column will be wrong if a non-printing ascii is there (please fix)


fid=fopen(filename);
R=fread(fid,'uint8=>char')';

S=regexprep(R,'\x9+',' ');       %replace tab with space
S=regexprep(S,'[\x0-\x8]+','');  %remove non-printing 
S=regexprep(S,'[\xB-\x1F]+',''); %remove non-printing 
S=regexprep(S,'[\x80-\x15F]+',''); %remove extended set

fclose(fid);

movefile(filename,[filename '_old']);

fid=fopen(filename,'w');

fprintf(fid,S);

output=1;
return

L=1;

line=0;
ex=0; np=0; st=0;

while L~=-1
    
    L=fgets(fid);  if L==-1, break, end
    line=line+1;
    Ldec=unicode2native(L);
    Lhex=dec2hex(Ldec);

    
    for i=1:1:length(Ldec)
        if Ldec(i)>=128, 
           character=native2unicode(Ldec(i));
           if length(sprintf('T%sT',character))==2, 
               character=' '; 
           end,
           fprintf('An extended ASCII character [%s]\t %03d (%s) was found on line %d column %d. \n',...
                             character ,Ldec(i),Lhex(i,1:2),line,i);
           ex=ex+1;
        elseif Ldec(i)<=31 && Ldec(i)~=13                                
           fprintf('Nonprinting ASCII character [ ]\t %03d (%s) was found on line %d column %d. \n',...
                                            Ldec(i),Lhex(i,:),line,i);
           np=np+1;                             
        elseif ~isempty(find(Ldec(i)==[33 34 35 37 39 42 47 59 61 62 63 64 91 92 93 94 95 96 123 124 125 126 127]))                                
           fprintf('A strange ASCII character   [%s]\t %03d (%s) was found on line %d column %d. \n',...
                             native2unicode(Ldec(i)) ,Ldec(i),Lhex(i,1:2),line,i); 
           st=st+1;              
        end
    end
    
end
       
    if sum(ex+np+st)==0,
         fprintf('This program found no characters that raised flags \n');
    else
         fprintf(['\n \n This program found \n \n'...
                  '  %d nonprinting ASCII characters \n'...
                  '  %d strange characters \n'...
                  '  %d extended ASCII characters \n \n'], np,st,ex);
    end
    
    fclose(fid)
 
return    
    
 % create test string
 % L=native2unicode(1:255)
 
 % create test file
 
  FID=fopen('test_file','w')
  
  for i=1:300;
     
      randstring=native2unicode(floor(256*rand(1,79)));
      fprintf(FID,'%s \n',randstring);
 
  end
  
  fclose(FID)