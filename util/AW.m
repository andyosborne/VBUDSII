function [Atomic_Weight]=AW(ZAID,parseXSDIR)
%
%  [Atomic_Weight]=AW(ZAID,parseXSDIR)
%  
%  Returns the atomic weight of a ZAID or vector of ZAIDs
% 
%  ex:   [Atomic_Weight]=AW([8016,92235,92238],0)
%
% If parseXSDIR==1, then it will parse the file xsdir.txt to obtain updated
% atomic weights.  Otherwise, it will keep an older precompiled table.
%
% ZAID=1  corresponds to 0001 (a neutron)
% ZAID=2  corresponds to      (an electron)
%  
% ZAID=ZZZAAA
%


%ANNOYINGLY xsdir actually gives the mass in units of Neutron masses to fix
%this we multiply the end result by N, the neutron mass in AMU
N=1.0086649156;  %Neutron mass in AMU

if nargin ~= 2, parseXSDIR=0; end;
if parseXSDIR==1, TABLE=parsefile; 
else load Atomic_Weight_Table, end

if iscell(ZAID) || ~isnumeric(ZAID)
    'ZAID must be a vector of numbers'
end

TABLE(2)=0.000543867302;  %electron (M_e/N)
TABLE(1)=1;               %neutron 

Atomic_Weight=N*full(TABLE(ZAID))';









%PARSE XSDIR (This is extra and not to be used often).
function [TABLE]=parsefile
F=fopen('xsdir.txt');
Line_String=fgets(F);
j=0;
TABLE=spalloc(1,118293,3300);  %3310 
while Line_String~=-1  
      Line_Cells=regexp(Line_String,'\d+\s+\d+.\d+','match');
      for i=1:length(Line_Cells)
          Line_Numbers=str2num(Line_Cells{i});
          TABLE(Line_Numbers(1))=Line_Numbers(2);
      end
      j=j+length(Line_Cells);  
      Line_String=fgets(F);
end
fclose(F);
save Atomic_Weight_Table TABLE