function [Tally, TallyProcessed] = TallyPull2(filename)
%
%  Created by: Geoffrey Recktenwald - Sept 2011
%              University of Texas at Austin
%
%   Goal:  This program takes an MCNP(5/X) output file and creates a
%          data structure for the tally output.
%
%   Data:  Tally{tally number}.
%               histories    Number of histories (normalization factor).
%               nps          Number of points
%               type         Tally type / final digit (1,2,4,etc).
%               cell_vol     String array of cell names to match volumes.
%               volumes      Numeric array of cell volumes for tally type 4.
%               surf_area    String array of surface identifiers to match areas.
%               areas        Numeric array of surface areas for tally type 2.
%               surfs        Cell array of the surface identifier strings.
%               cells        Cell array of the cell identifier strings.
%               multi        Cell array of the multiplication factor strings.
%               coses        Cell array of cosine bins strings.
%               energy       Numeric array of the energy bins.
%       ==>     value        Numeric vector of the tally.
%               error        Numeric vector of the tally error.
%               text         Character array of the tally text.
%               forms        Character array of condensed tally information.
%
%   Notes:
%       cell_vol and cells should have the same information in the same order.
%          surf_area and surfs should have the same information in the same order.
%          serface replaces surface because matlab has a surface command
%
%      Call [Tally]=Tallypull('outp')
%
%   perhaps escape characters should just be L(1)='1';
%  111026: Dembia made 2 changes to adapt to MCNP5: detecting MCNP version, the
%  columns where tally numbers are found, and detecting the end of the
%  1materials portion, moved 1analysis break to the top of the loop because it
%  interfered with L = -1, commented out 'save Tallypullsave', started
%  counting linenums, changed the way that the numeric card is managed (turning
%  on a boolean)
%
Tally{1}.void='empty';                 tallynum=[];
fid=fopen(filename);                   n=150;   %L will always be length 150
L=fgetG(fid,n);                        it=0;
linenum = 1;

grabNumericCard = 1;

tallyend2 = ' there ar';
if strfind(L,'MCNP5')
    mcnpversion = '5';
    tallyend1 = ' ========';
else
    mcnpversion = 'X';
    tallyend1 = ' fom = (h';
end

while L~=-1

    if ~isempty(strfind(L, 'fatal error'))
        error(['This output file is for a run that had a fatal ' ...
            'error:\n%s\n'], L);
    end
    
    [num,nps]=findtallystart(L, mcnpversion);
    enter_materials=strcmp(L(1:9),'1material');

    if num~=0
        tallynum=[tallynum num]; %(assumes no repetition)
        it=it+1;
        Tally{it}.num=num;
        Tally{it}.nps=nps;
        Tally{it}.type=mod(num,10);    %if strfind(L,'tally type'), gettallytype, end

        %initialize
        is=1; im=1; ic=1; ie=1; Ie=0; i=0;       %indicies and steps
        cells={}; surfs={};  multi={}; coses={};  %maps for indicies
        surfarea=[]; volume=[]; surf_area=''; cell_vol=''; %voids for loops

        while L~=-1 & ~strcmp(L(1:9),tallyend1) & ~strcmp(L(1:9),tallyend2)
            % fprintf('%s\n',L), pause(0.01)
            Card='';
            if strfind(L,'histories'),   Card='histories';  ,end
            if strcmp(L(2:8),'surface'), Card='serface';    ,end
            if strcmp(L(2:5),'cell'),    Card='cell';       ,end
            if strcmp(L(2:8),'multipl'), Card='multiplier'; ,end
            if strcmp(L(2:7),'cosine'),  Card='cosine';     ,end
            if strfind(L,'cell:'),       Card='volumes';    ,end
            if strfind(L,'surface:'),    Card='serfaces';   ,end
            if strfind(L,'energy'),      Ie=1;              ,end
            if strfind(L,'total'),       Card='total';      ,end
            % if str2num(L),             Card='numeric';    end
                    %this one is dangerous
            N1=isempty(regexp(L,'[^1234567890 Ee+\-\.]'));
            % I contain only numeric characters  ()[]
            N2=isempty(regexp(L,'\d'));
            % I actually have numbers
            if 0 && ~isempty(Card)
                linenum
                Card
            end

            if N1 && ~N2
                Card='numeric';
            end

            switch Card
                % data is stored as   Tally{it}.value{serface,multi}(energy,cos)
                case 'histories'
                    Tally{it}.histories=getnum(L);
                case 'volumes'
                % Presumed that the order here is the same as the order in
                % the actual tallies. To check this, compare
                % Tally.cell_vol with Tally.cells.
                    cell_vol=[cell_vol,L(25:end)];
                    Tally{it}.cell_vol=cell_vol;
                    i=i+1;
                    Tally{it}.Text(i,:)=L;
                    L=fgetG(fid,n);
                    linenum = linenum + 1;
                    volume=[volume,str2num(L)];
                    Tally{it}.volume=volume;
                case 'serfaces'
                    surf_area=[surf_area,L(25:end)];
                    Tally{it}.surf_area=surf_area;
                    i=i+1;
                    Tally{it}.Text(i,:)=L;
                    L=fgetG(fid,n);
                    linenum = linenum + 1;
                    surfarea=[surfarea,str2num(L)];
                    Tally{it}.areas=surfarea;
                case 'serface'
                    serfacename=L(9:end);
                    is=find(strcmp(surfs,serfacename));
                    if isempty(is)
                        is=length(surfs)+1;
                        surfs{is}=serfacename;
                    end
                    Tally{it}.surfs=surfs;
                case 'cell'
                    cellname=L(6:end);
                    is=find(strcmp(cells,cellname));
                    if isempty(is)
                        is=length(cells)+1;
                        cells{is}=cellname;
                    end
                    Tally{it}.cells=cells;
                case 'multiplier'
                    multiplier=L(17:end);
                    im=find(strcmp(multi,multiplier));
                    if isempty(im)
                        im=length(multi)+1;
                        ie = 1; % reset counter
                        multi{im}=multiplier;
                    end
                    Tally{it}.multi=multi;
                case 'cosine'
                    cosbin=L(15:end);
                    ic=find(strcmp(coses,cosbin));
                    if isempty(ic)
                        ic=length(coses)+1;
                        coses{ic}=cosbin;
                    end
                    Tally{it}.coses=coses;
                    %Tally{it}.cosines=
                    %Please remember Cos # to #   multiplier #
                case 'numeric'
                    V=str2num(L);
                    if length(V)==3,
                        Tally{it}.energy{is,im}(ie,ic)=V(1);
                        Tally{it}.value{is,im}(ie,ic) =V(2);
                        Tally{it}.error{is,im}(ie,ic) =V(3);
                    elseif length(V)==2,
                        Tally{it}.value{is,im}(ie,ic)=V(1);
                        Tally{it}.error{is,im}(ie,ic)=V(2);
                    else
                        fprintf('Line %d: V should be length 2 or 3\n',linenum)

                    end
                    ie=ie+Ie; %step forward in energy if needed;
                case 'total'  %Total is stored on the end of the energy
                    V=str2num(L(14:end));
                    Tally{it}.energy{is,im}(ie,ic)=NaN;
                    Tally{it}.value{is,im}(ie,ic) =V(1);
                    Tally{it}.error{is,im}(ie,ic) =V(2);
                otherwise
            end %switch
            i=i+1;

            Tally{it}.Text(i,:)=L;

            % save Tallypullsave
            L=fgetG(fid,n);
            linenum = linenum + 1;

        end %while L ~= -1

        Tally{it}.forms=createFORMs(num,surfs,cells,multi,coses,tallynum);

        num=0;
    end %if num ~= 0

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if enter_materials
        atommass=0;  mo={}; ao={}; matnum=0;
        while L~=-1
            if strfind(L,'component nuclide, atom fraction')
                atommass='atom'; mat_num={}; iso_vec={}; matnum=0;
            elseif strfind(L,'component nuclide, mass fraction')
                atommass='mass'; mat_num={}; iso_vec={}; matnum=0;
            end

            mat_num=regexp(L,'  \d+  ','match');
            if ~isempty(mat_num)
                matnum=str2num(mat_num{1});
            end

            iso_vec=regexp(L,' \d+, [\dEe+\-\.]+ ','match');
            if matnum~=0 && isempty(mat_num)
                isovec=[isovec;str2num(char(iso_vec))];
            else
                isovec=str2num(char(iso_vec));
            end

            if matnum
                switch atommass
                    case 'mass'
                        mo{ matnum } = isovec;
                    case 'atom'
                        ao{ matnum } = isovec;
                end
            end

            L=fgetG(fid,n);
            linenum = linenum + 1;
            if strcmp(L(1:8),'1LAHET p') || strcmp(L(1:6),'1cells')
                             % mcnpx                        mcnp5
                break
            end  %escape
        end
        enter_materials==0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    L=fgetG(fid,n);
    linenum = linenum + 1;
end %while

it=it+1;
for j=1:1:it-1;
    Tally{it}.forms{j}=Tally{j}.forms;
end

Tally{1}=rmfield(Tally{1},'void');
Tally{it}.nums=tallynum;
Tally{it}.it=[1:it-1];
Tally{it}.mo=mo;
Tally{it}.ao=ao;

fclose(fid);

TallyProcessed = postprocess(Tally);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [num,nps]=findtallystart(L, mcnpversion)
% This function looks in string L to find the line '1tally', which
% indicates that the script has reached a tally.  (excusions apply)
num=0; nps=0;
if strcmp(L(1:6),'1tally')                     % Find the tally start
    if strcmp(L(1:10),'1tally flu'), return, end %  ignore fluxuation
    if strfind(L,'print table 30'), return, end  %  ignore
    if mcnpversion == 'X'
        numidxer = 7:10;
    elseif mcnpversion == '5'
        numidxer = 12:15;
    end
    num=str2num(L(numidxer));     % Return the tally number
        %alt: regexp(L, '1tally\s*(?<term\d+>)','names');
    nps=getnum(L,0);           % Return the number of points.
end

end

function [value]=getnum(L,n)
data=regexp(L,'\s\d[\d\.eE+-]+\s','match');
if ~exist('n')
    n=1;
end

if n>=1
    value=str2num(data{n});
elseif n==0
    value=str2num(data{end});
else
    for i=1:1:length(data)
        value(i)=str2num(data{i});
    end
end

end


function [output]=createFORMs(num,surfs,cells,multi,coses,tallynum)
if mod(num,10)==1 || mod(num,10)==2
    R1=strtrim(surfs); R{1}=char([{' Surfaces'},    R1]); [A(1),B(1)]=size(R{1});
elseif mod(num,10)==4
    R1=strtrim(cells); R{1}=char([{' Cells'},       R1]); [A(1),B(1)]=size(R{1});
else
    'Error on num type'
end
R2=strtrim(multi); R{2}=char([{' Multipliers'}, R2]); [A(2),B(2)]=size(R{2});
R3=strtrim(coses); R{3}=char([{' Cosines'},     R3]); [A(3),B(3)]=size(R{3});
R{4}=[ ' Tallies'; reshape(sprintf('   %03d  ',tallynum),8,length(tallynum))'];   %'
[A(4),B(4)]=size(R{4});

Am=max(A); Bm=max(B);

P(1:Am+1,1)='|';
S(1:Am+1,1)=' ';

for i=1:1:4
    R{i}(Am+1,1)=' ';
    R{i}(find(R{i}==native2unicode(0)))=' ';
end

q1='Tally{tallies}.value{surf or cell,multiplier}(energies,angles)';
q2='--------------------------------------------------------------';
q=[q1;q2];

numerMulti = num2str((0:Am)');
formarray=[R{4}, S, P, S R{1}, S, P, S, numerMulti, S, R{2}, S, P, S, R{3}];

RR=length(q);
[A,B]=size(formarray);

if     B>RR, q(1,RR+1:B)=' ';
elseif B<RR, formarray(:,B+1:RR)=' ';
end

output=[q;formarray];

end

function [TallyProcessed]=postprocess(Tally)
% The goal of this function is to postprocess the tally information to
% determine the desired cross sections.  It can format the cross sections
% for input into matlab or maple.
%
% I will presume, for the sake of this code, that cross section tallies
% are created by calling materials with only 1 ZAID ( h20?  Total?)
%

keys0 = cell(1,length(Tally)-1);
values0 = cell(1,length(Tally)-1);
ts = zeros(1,length(Tally)-1);
for tallyidx = 1:(length(Tally) - 1)
    keys0{tallyidx} = Tally{tallyidx}.num;
    ts(tallyidx) = Tally{tallyidx}.num;
%    values0{tallyidx} = zeros(length(Tally{1}.value{1},1);
end

T = struct('fs',cell2mat(keys0),...
           'f',struct('Ms',[],...
                      'Bs',[],...
                      'B2s',[],...
                      'M',[],...
                      'B',[],...
                      'B2',[],...
                      'm',struct('b',struct('value',[],...
                                            'b2',struct('value',[])))));


matlcards2 = [];
B12 = [];
B22 = [];
counter = 0;
    keys = {}; %cell(1,(length(Tally)-1)*nLists);
    values = {}; % cell(1,(length(Tally)-1)*nLists);
for tallyidx = 1:(length(Tally) - 1)

    matlcards = [];
    B1 = [];
    B2 = [];
    nLists = length(Tally{tallyidx}.value);
    for listidx = 1:nLists
        
        counter = counter + 1;

        M = str2num(Tally{tallyidx}.multi{listidx});

            %M1=regexp(M,' [()-\d]+ ','match')
            %M2=regexp(M,'(\d+)','match');
        C = M(1);
        if length(M) == 1
            % we're looking at flux here probably.
            matlcard = 0;
        else
            matlcard = M(2);
        end

        if length(M) >= 3
            B0 = M(3);
        else
            B0 = 0;
        end
        B1 = [B1 B0];

        keys{counter} = [num2str(keys0{tallyidx}) '.' num2str(matlcard) '.' num2str(B0)];
        if length(M) == 4
            B2 = [B2 M(4)];
            keys{counter} = [keys{counter} '.'  num2str(B2(listidx))];
        else
            B2 = [B2 0];
        end

        if isempty(find(matlcards == matlcard))
            matlcards = [matlcards matlcard];
        end

        values{counter} = Tally{tallyidx}.value{listidx};
    end
    T.f(tallyidx).Ms = matlcards;
    T.f(tallyidx).Bs = B1;
    T.f(tallyidx).B2s = B2;

    T.f(tallyidx).M = containers.Map(num2cell(matlcards), num2cell(1:length(matlcards)));
    T.f(tallyidx).B = containers.Map(num2cell(B1), num2cell(1:length(B1)));
    T.f(tallyidx).B2 = containers.Map(num2cell(B2), num2cell(1:length(B2)));
    
    %TallyProcessed(keys{tallyidx}) = containers.Map
    values0{tallyidx} = containers.Map( keys, values);

    matlcards2 = [matlcards2 matlcards];
    B12 = [B12 B1];
    B22 = [B22 B2];

    if 0
    for listidx = 1:nLists
        if B2 == 0
        T.f(tallyidx).m( T.f(tallyidx)).M(matlcards(listidx)).b( T.f(tallyidx)).B(B1(listidx)).value = Tally{tallyidx}.value{listidx};
        else
        T.f(tallyidx).m( T.f(tallyidx)).M(matlcards(listidx)).b( T.f(tallyidx)).B(B1(listidx)).b2(T.f(tallyidx)).B2(B2(listidx)).value = Tally{tallyidx}.value{listidx};
        end
    end
    end

end

TallyProcessed = containers.Map(keys,values);

end

%T = struct('fs',fs,...
%           'f',struct('Ms',[],...
%                      'Bs',[],...
%                      'B2s',[],...
%                      'M',[],...
%                      'B',[],...
%                      'B2',[],...
%                      'm',struct('b',struct('value',[],...
%                                            'b2',struct('value',[])))));
%
%
%T.M(1:length(matlcards2)) = matlcards2;
%T.M = sparse(T.M);
%
%T.B = container.Map( num2cell(1:length(B12)), num2cell(B12) );
%T.B2 = container.Map( num2cell(1:length(B22)), num2cell(B22) );
%


                %switch strtrim(M(36:47))
                %    case '1'
                %    case '16'
                %    case '18'  %fission
                %    case '102' %gamma
                %    case '103'
                %    case '107'
                %    case '-1' %Total
                %    case '-2' %capture
                %    case '-3' %elasic cattering
                %    case '-4' %heating
                %    case '-6' %Total
                %    case '-7' %neutrons/fision
                %    case '-8' %Q (Mev/fission)
                %    otherwise
                %end

                %find cell -> associate with zaids
                %find multiplier -> associate with reaction
                % check for energies or angles
                %fprintf(['sigma[%d,%s]=%3.5f',zaid,reaction,Tally{i1}.value{i2,i3})


%    matl = cell(1,length(matlcards)); 
%    for idx = 1:length(matlcards)
%        matl{idx} = matlcards(idx); 
%    end
%    B = cell(1,length(B1));
%    for idx = 1:length(B1)
%        B{idx} = B1(idx);
%    end
