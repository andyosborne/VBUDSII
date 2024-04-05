function [Meta, L, p] = FarmXS(p)
% TODO:
%   1- Allow for other s0's. 'infinite dilution'

% BOOLEAN TO SHOW ONLY ERRORS OR TO ALSO SHOW ALL NEWS.
% IF IGNstr IS NOT '1', THEN MUST AVOID USING L.GROUPDEF, L.NGRUPS,
% P.P.NGROUPS

disp('It made it here')
disp(p)
if exist([p.metamat '.mat'],'file') == 2
    load(p.metamat);
end

L = 0;
p.fnames = {'tape20', 'script', 'output'};

tic;

% explore: read down the available ENDF/B files from the internet.
if p.explore == 1
    if isfield(p, 'endfversion') && p.endfversion == 7.1
        Meta = ExploreVIIp1(p);
    else
        Meta = Explore(p);
    end
    save(p.metamat,'Meta');
end

% if user wants water (222), we must make sure to request ZAIDs 1001 and 8016.
if sum(p.ZAIDs == 222) > 0
    if sum(p.ZAIDs == 1001) == 0
        p.ZAIDs = [p.ZAIDs 1001];
    end
    if sum(p.ZAIDs == 8016) == 0
        p.ZAIDs = [p.ZAIDs 8016];
    end
end

% mine: pull down ENDF/B text files from the internet.
if p.mine == 1
    Meta = Mine(p, Meta);
    save(p.metamat,'Meta');
end

% plant: plant NJOY scripts.
if p.plant == 1
    Meta = Plant(p, Meta);
    save(p.metamat,'Meta');
end

% fertilize: run NJOY scripts and create NJOY output files.
if p.fertilize == 1
    [Meta, fileStatuses] = Fertilize(p, Meta);
    save(p.metamat,'Meta');
end

% harvest: run MATLAB on NJOY output files to provide cross sections for
% human consumption.
if p.harvest == 1
    [Meta, L, njoysMTs] = Harvest(p, Meta);
    save(p.metamat,'Meta');

    if isfield(p,'combineinelastic') && p.combineinelastic
        [Meta, L] = CombineInelastic(p, Meta, L);
    end

    if sum(p.ZAIDs == 222) > 0
        % water is requested. harvest it!
        [Meta, L] = MakeWater(p, Meta, L);
    end
end

if isfield(p, 'cook') && p.cook

    [Meta, L] = Cook(p, Meta, L, njoysMTs);
end

fprintf('Runtime: %.2f seconds\n',toc);



end



function [Meta, L] = Cook(p, Meta, L, njoysMTs)

if isfield(p.c, 'matperzaid') && p.c.matperzaid && isfield(p.c, 'matpath')
    % L is the full library with all the data, Lib is just the metadata.
    % Eventually this should all be written so that at no one point does all
    % the data need to be loaded into memory.
    matpath = p.c.matpath;
    Lib = initLibraryMeta(p, Meta, njoysMTs);

    if ~strcmp(matpath(end), filesep)
        matpath = [matpath filesep];
    end
    for z = L.ZAIDs
        if z == 1001 && isfield(p.c, 'correctlwtr')
            L.z(L.ZAID(z)).m(L.MT(2)).xs = p.c.correctlwtr;
        end
        ZAIDtables = L.z(L.ZAID(z));
        fname = [matpath 'zaid' num2str(z)];
        Lib.z(L.ZAID(z)).matpath = fname;
        save(fname, 'ZAIDtables');
    end
    save([matpath 'librarymetadata'], 'Lib');
end

end

function Meta = initMeta()
Meta = struct('MAT',[],...
    'isIsomer',[],...
    't2name',[],...
    'numProtons',[],...
    'element',[],...
    'atomNum',[],...
    'ZAID',[],...
    'webPath',[],...
    'nuclideStr',[],...
    'storagePath',[],...
    'writeStatus',[],...
    'isFissionable',[],...
    'maxInelasticMT',[],...
    'MTs',[],...
    'fileStatus',[]);
end

function Meta = ExploreVIIp1(p)
if p.verbose
    disp('[FarmXS/ExploreVIIp1]: Entering.');
end

t2domain = 'http://t2.lanl.gov';
htmlfname = 'endfvii.1-n.html';

if p.e.searchWeb
    urlwrite([t2domain '/data/endf/' htmlfname], htmlfname);
    fid = fopen(htmlfname);
else
    Stack = dbstack('-completenames');
    thisDir = fileparts(Stack(1).file);
    fid = fopen(fullfile(thisDir, htmlfname));
end

if fid == -1
    error(['ERROR: html file current.html is not in the current directory.'...
        ' Set p.e.searchWeb to 1.']);
end

Meta = initMeta();

% Sample entries
% <TR>
% <TD><I><B>67-Ho</B></I></TD>
%   <TD><I><B>-165</B></I></TD>
%     <TD>6725</TD>
%     <TD>-M-4----</TD>
%     <TD><A HREF="/data/data/ENDFB-VII.1-neutron/Ho/165">raw eval</A></TD>
%     <TD>
%     <FORM METHOD="POST" ACTION="/data/ndbrowse.php">
%     <INPUT TYPE=hidden NAME="path" VALUE="data/ENDFB-VII.1-neutron/Ho/165">
%     <INPUT TYPE=hidden NAME="mf" VALUE=0>
%     <INPUT TYPE=hidden NAME="mt" VALUE=0>
%     <INPUT TYPE="submit" VALUE="browse eval">
%     </FORM>
%      <TD><A HREF="/data/endf/endfvii.1-n-pdf/ho165.pdf">view PDF plots</A></TD>
% </TR>
% <TR>
% <TR><TD><BR></TD>
%   <TD><I><B>-166m1</B></I></TD>
%     <TD>6729</TD>
%     <TD>-MU46-g-</TD>
%     <TD><A HREF="/data/data/ENDFB-VII.1-neutron/Ho/166m1">raw eval</A></TD>
%     <TD>
%     <FORM METHOD="POST" ACTION="/data/ndbrowse.php">
%     <INPUT TYPE=hidden NAME="path" VALUE="data/ENDFB-VII.1-neutron/Ho/166m1">
%     <INPUT TYPE=hidden NAME="mf" VALUE=0>
%     <INPUT TYPE=hidden NAME="mt" VALUE=0>
%     <INPUT TYPE="submit" VALUE="browse eval">
%     </FORM>
%      <TD><A HREF="/data/endf/endfvii.1-n-pdf/ho166m1.pdf">view PDF plots</A></TD>
% </TR>

% Index for MATs.
iM = 0;
while ~feof(fid)
    fline = fgetl(fid);
    if length(fline) > 12 && ~isempty(strfind(fline, '<TD><I><B>'))
        % This string occurs at the start of every line that starts the entry
        % of a new isotope. Start parsing for this nuclide.
        iM = iM + 1;

        if ~isempty(strfind(fline, '<B>-'))
            % This is not the first isotope of this element.
            Meta(iM).numProtons = Meta(iM-1).numProtons;
            Meta(iM).element = Meta(iM-1).element;
        else
            % Ex: <TD><I><B>67-Ho</B></I></TD>
            % Get proton number on this line. It's found right after the <B> tag
            % and before the elemental symbol.
            iLeft = strfind(fline, '<B>') + 3;
            iRight = strfind(fline, '-') - 1;
            Meta(iM).numProtons = str2num(fline(iLeft(1):iRight(1)));

            % Element name. Right after the proton number, before any tags.
            iLeft = strfind(fline, '-') + 1;
            iRight = strfind(fline, '</') - 1;
            % Upper-case for consistency with the ENDF-VI parser.
            Meta(iM).element = upper(strtrim(fline(iLeft(1):iRight(1))));
    
            % Move onto the next line.
            % Ex: <TD><I><B>-165</B></I></TD>
            fline = fgetl(fid);
        end

        % Atomic number. Between the tags.
        % May need to deal with an isomeric state.
        % TODO we're assuming there's only one isomeric state for a given
        % proton and atom number.
        iLeft = strfind(fline, '-') + 1;
        iRight = strfind(fline, '</B>') - 1;
        if strcmp(fline(iLeft(1):iRight(1)), 'nat')
            % The naturally occuring version.
            Meta(iM).atomNum = 0;
        else
            if ~isempty(strfind(fline(iLeft(1):iRight(1)), 'm'))
                Meta(iM).isIsomer = true;
                iRight = strfind(fline, 'm') - 1; 
            else
                Meta(iM).isIsomer = false;
            end
            Meta(iM).atomNum = str2num(fline(iLeft(1):iRight(1)));
        end
        if Meta(iM).isIsomer
            % Move to the right of the 'm'.
            idx = iRight + 2;
            stateNo = str2num(fline(idx));
        else
            stateNo = -3;
        end
        % Double zero for nonisomeric states!
        Meta(iM).ZAID = str2num(sprintf('%02i%03i', Meta(iM).numProtons, ...
            Meta(iM).atomNum + Meta(iM).isIsomer * (stateNo + 3) * 100));

        % Used to name files.
        Meta(iM).nuclideStr = [Meta(iM).element '_' num2str(Meta(iM).atomNum)];
        if Meta(iM).isIsomer
            Meta(iM).nuclideStr = [Meta(iM).nuclideStr 'm'];
        end

        % Move onto the next line.
        % Ex:    <TD>6729</TD>
        fline = fgetl(fid);

        % ENDF/B MAT number.
        iLeft = strfind(fline, '>') + 1;
        iRight = strfind(fline, '</') - 1;
        Meta(iM).MAT = str2num(fline(iLeft(1):iRight(1)));

        % Skip a line!
        fline = fgetl(fid);
        fline = fgetl(fid);
   % Ex: <TD><A HREF="/data/data/ENDFB-VII.1-neutron/Ho/165">raw eval</A></TD>

        % link to the tape!
        iLeft = strfind(fline, '"') + 1;
        iRight = strfind(fline, 'raw') - 3;
        Meta(iM).webPath = [t2domain fline(iLeft(1):iRight(1))];
    end
end

fclose(fid);

if p.verbose
    disp('[FarmXS/ExploreVIIp1]: Leaving.');
end

end % function Meta = ExploreVIIp1(p)

function Meta = Explore(p)
%MINE obtains information about ENDF/B files, and downloads ENDF/B files. This
%function does not have anything to do with NJOY.

if p.verbose
    disp('[FarmXS/Explore]: Entering.');
end

if p.e.searchWeb
    urlwrite('http://t2.lanl.gov/cgi-bin/nuclides/endind','current.html');
    fid = fopen('current.html');
else
    Stack = dbstack('-completenames');
    thisDir = fileparts(Stack(1).file);
    % Chris Dembia renamed the local copy of "current.html" to 
    % "endfvi_index.html".
    fid = fopen(fullfile(thisDir, 'endfvi_index.html'));
end

if fid == -1
    error(['ERROR: html file current.html is not in the current directory.'...
        ' Set p.e.searchWeb to 1.']);
end

Meta = initMeta();

% Line number. EDIT.
m = 0;

%% Obtain meta data about all the ENDF/B files listed.
while ~feof(fid)
    fline = fgetl(fid);
    if strcmp(fline(1:2),'<A')
        % Line contains a URL, and so it contains info for an ENFB/F file.

        m = m + 1;

        % Index of the first character of the string for the filename of
        % the webpage where the ENFB/F file is located.
        zidx1 = strfind(fline,'?') + 1;

        % Index of the last character.
        zidx2 = strfind(fline,'">') - 1;

        if m > 1 && ...
                str2num( fline(zidx2 + 9: zidx2 + 14)) == Meta(m-1).MAT &&...
                strcmp( fline( zidx2 + 26), 'M') == Meta(m-1).isIsomer
            % The current file is for the same nuclide; go back one line and
            % use this file instead of the previous one. Files are ordered from
            % oldest to newest, so this ensures that the newest file for each
            % nuclide is used.
            m = m - 1;
        end

        % Extract, for this ENDF/B file, meta data about this file.

        % t2.lanl.gov name for the file.
        Meta(m).t2name = fline(zidx1:zidx2);

        % ENDF/B mat number for the file.
        Lidx = zidx2 + 9; Ridx = Lidx + 5;
        Meta(m).MAT = str2num( fline(Lidx:Ridx));

        % Proton number for the file.
        Lidx = zidx2 + 16; Ridx = Lidx + 2;
        Meta(m).numProtons = str2num( fline(Lidx:Ridx) );

        % Element name for the file.
        Lidx = zidx2 + 20; Ridx = Lidx + 1;
        Meta(m).element = strtrim( fline(Lidx:Ridx) );

        % Atomic number for the file.
        Lidx = zidx2 + 23; Ridx = Lidx + 2;
        Meta(m).atomNum = str2num( fline(Lidx:Ridx) );

        % Isomer index for the file.
        Lidx = zidx2 + 26;
        Meta(m).isIsomer = strcmp( fline(Lidx), 'M' );

        % ZAID for the file.
        % Get the correct number of zeros.
        prezaid = num2str(10^(3 - length(num2str(Meta(m).atomNum))));
        Meta(m).ZAID = str2num([num2str(Meta(m).numProtons) prezaid(2:end)...
            num2str(Meta(m).atomNum)]);

        % IF THIS ZAID MATCHES A REQUESTED ZAID, PUT IT INTO A SEPARATE LIST.

        % Obtain the web addresses of the ENDF/B files.
        if Meta(m).atomNum == 0
            % This file is for the naturally-occuring nuclide, not a specific
            % isotope.
            % Index of the last character of the URL
            urlidx = strfind( Meta(m).t2name, 'nat') - 1;
        else
            urlidx = strfind( Meta(m).t2name, int2str( Meta(m).atomNum )) - 1;
        end
        % The webpath contains a directory name, then a file name. The directory
        % name is just the element abbreviation, but with the first letter
        % capitalized and the subsequent letters as lower case. After the slash
        % is the atomic number (or nata / natf for files for natural nuclides).
        Meta(m).webPath = ['http://t2.lanl.gov/data/data/ENDFB-VI-neutron/'...
            Meta(m).element(1) Meta(m).t2name(2:urlidx) ...
            '/' ...
            Meta(m).t2name(urlidx + 1:end)];

        Meta(m).nuclideStr = [Meta(m).element '_' num2str(Meta(m).atomNum)];
        if Meta(m).isIsomer
            Meta(m).nuclideStr = [Meta(m).nuclideStr 'm'];
        end
        % DEAL WITH ISOMERS. EDIT.
        % IS ISOTOPE ON THE FISSILE LIST? EDIT.
    end
end

% bound materials EDIT.

% Close current.html.
fclose(fid);

if p.verbose
    disp('[FarmXS/Explore]: Leaving.');
end

end

function Meta = Mine(p, Meta)
% Download ENDF/B files, for the ZAIDs requested, in the p.address folder.

if p.verbose
    disp('[FarmXS/Mine]: Entering.');
end

mainidxs = MapZAIDs(Meta, p.ZAIDs);

for m = mainidxs % EDIT.

    Meta(m).storagePath = fullfile(p.address, Meta(m).nuclideStr);

    if p.m.downloadTape20s

        if ~isdir(Meta(m).storagePath)
            mkdir(Meta(m).storagePath);
            disp(['[FarmXS/Mine]: Created a new directory ' ...
                p.address filesep Meta(m).nuclideStr '.']);
        end

        if Meta(m).ZAID == 26056
            % endf7 path TODO
            Meta(m).webPath = ...
                'http://t2.lanl.gov/data/data/ENDFB-VII-neutron/Fe/56';
        end
        [F, status] = urlwrite(Meta(m).webPath, ...
            fullfile(Meta(m).storagePath, ...
            'tape20'));

        % If there was no error downloading the file, status = 0.
        Meta(m).writeStatus = status;

        if status ~= 1
            disp(['WARNING [FarmXS/Mine]: ' Meta(m).nuclideStr ', ZAID '...
                num2str(Meta(m).ZAID) ': URLWRITE error code ' ...
                num2str(status)]);
        end
    else

        % Value is normally 2 is the file exists, so we take half of that.
        Meta(m).writeStatus = .5*exist(fullfile(Meta(m).storagePath, 'tape20'));

    end
end

% TODO List the files that failed to download.
if p.verbose
    disp('[FarmXS/Mine]: Leaving.');
end

end


function Meta = Plant(p, Meta)

if p.verbose
    disp('[FarmXS/Plant]: Entering.');
end

nTs = length(p.p.Ts);
nS0s = length(p.p.S0s);

% Error check input
if isfield(p.p, 'legendreOrder')
    legOrder = p.p.legendreOrder;
else
    legOrder = '0';
end

p.p.IGNstr = num2str(p.p.IGNstr);
p.p.IWTstr = num2str(p.p.IWTstr);

maxNJOYGroups = 500;

minInelasticMT = 51;

if ~strcmp(p.p.IGNstr, '1') && (length(p.p.groupDef)-1) > maxNJOYGroups
    % Eventually automate the creation of multiple scripts.
    error(['ERROR [FarmXS/Plant]: current number of groups is ' ...
        num2str(length(p.p.groupDef)-1) ', maximum number of groups '...
        'per script is ' num2str(maxNJOYGroups) '. ' ...
        'Break run into multiple scripts']);
end

mainidxs = MapZAIDs(Meta, p.ZAIDs); % Call this every time?


for m = mainidxs % EDIT.

    [Meta(m).isFissionable, Meta(m).maxInelasticMT, Meta(m).MTs] = ...
        GetMFMTs( Meta(m).storagePath, Meta(m).writeStatus);

    % Inspect comment section of ENDF/B file to update Fissile and
    % maxinelastic. EDIT.
    % THIS FUNCTION SHOULD REALLY BE INCORPORATED INTO EXPLORE.

    % DEAL WITH THE FACT THAT WATER IS AVAILABLE AT ONLY CERTAIN TEMPERATURES.

    % Select MFMTs.

    MFMTs = p.p.MFMTs;
    if isfield(p.p,'MFs');
        MFs = p.p.MFs;
    end

    % Write the NJOY script!
    if Meta(m).writeStatus ~= 1
        disp(['WARNING [FarmXS/Plant]: Not generating a script for ' ...
            Meta(m).nuclideStr ' ZAID ' num2str(Meta(m).ZAID) ...
            ' because of a write error.']);
    else
        fid = fopen( fullfile(p.address, Meta(m).nuclideStr,...
            'script'),'wt');

        MATstr = int2str(Meta(m).MAT);
        Tstr = strtrim( sprintf('%d ', p.p.Ts));
        nTstr = int2str(nTs);
        % THIS LINE HERE IS PROBABLY NOT OKAY FOR NJOY NUMBER INPUTS. EDIT.
        S0str = strtrim( sprintf('%1.1e ', p.p.S0s));
        nS0str = int2str(nS0s);

        fprintf(fid,[...
            ' moder \n' ... % MODER - converts tape20 to binary
            ' 20 -21/ \n' ... % input, output tape number
            ' reconr \n' ... % RECONR
            ' -21 -22/ \n' ... % input, output tape number
            ' ''pendf file''/ \n' ... % label
            ' ' MATstr ' 1/ \n' ... % ENDF/B mat, # cards
            ' .05/ \n' ... % tolerance
            ' ''<pendf>-file''/ \n' ... % label
            ' 0/ \n' ... % end RECONR
            ' broadr \n' ... % BROADR
            ' -21 -22 -23/ \n' ... %input, input, output
            ' ' MATstr ' ' nTstr '/ \n' ... % ENDF/B mat, # temps
            ' .005/ \n' ... % tolerance (thinning)
            ' ' Tstr '/ \n' ... % temperatures (Kelvin)
            ' 0/ \n' ... % terminate BROADR
            ' unresr \n' ... % UNRESR
            ' -21 -23 -24/ \n' ... % input, input, output
            ' ' MATstr ' ' nTstr ' ' nS0str ' 1/ \n' ... % # temps, # S0s
            ' ' Tstr '/ \n' ... % temperatures (Kelvin)
            ' ' S0str '/ \n' ... % S0s (cm^{-1})
            ' 0/ \n' ... % terminate UNRESR
            ' thermr \n' ... % THERMR
            ' 0 -24 -25/ \n' ... % ENDF input, PENDF input, output
            ...   % (ENDF not needed)
            ' 0 ' MATstr ' 16 ' nTstr ' 1 0 1 221 1/ \n' ...
            ...   % ENDF mat, PENDF mat, # angles, # temps, inelastic,
            ...   % elastic, nAtoms, inelasticMT, verbose
            ' ' Tstr '/ \n' ... % temperatures (Kelvin)
            ' .001 10./ \n' ... % tolerance, max energy for thermal treatment
            ' groupr \n' ... % GROUPR
            ' 20 -25 0 -26/ \n' ... % ENDF input, PENDF input, NGOUT1, NGOUT2
            ' ' MATstr ' ' p.p.IGNstr ' 0 ' p.p.IWTstr ' ' legOrder ' ' nTstr ' ' nS0str '/ \n' ...
            ...   % ENDF mat, IGN (groupdef), IGG (gamma), IWT (weight),
            ...   % Legendre order, # temps, # S0s
            ' ''Group constants''/ \n' ... % label
            ' ' Tstr '/ \n' ... % temperatures (Kelvin)
            ' ' S0str '/ \n']); % S0s (cm^{-1})
        if strcmp(p.p.IGNstr,'1')
            % Using our own group structure. p.p.groupDef must be defined.
            fprintf(fid,[' ' int2str(length(p.p.groupDef)-1) '/ \n']); % nGroups
            fprintf(fid, '%1.8f ', p.p.groupDef); % print groupDef
            fprintf(fid, '/ \n'); % End groupDef line.
        end

        % Exception handling.
        if ~Meta(m).isFissionable
            % Do not find MTs related to only fissionable isotopes.
            for mt = [18 19 20 21 452] % There are even more! 453? TODO
                if find(MFMTs == mt)
                    [rowidx colidx] = find(MFMTs == mt);
                    MFMTs = MFMTs([1:rowidx-1 rowidx+1:end], :);
                end
            end
        end

        if (sum(Meta(m).ZAID == 1001) >= 1) && (sum(MFMTs(:,2) == 4) == 1)
            % Water does not have inelastic scattering.
            disp('WARNING [FarmXS/Plant]: ZAID 1001 does not have MT 4');
            [rowidx colidx] = find(MFMTs == 4);
            MFMTs = MFMTs([1:rowidx-1 rowidx+1:end], :);
        end

        % Put in the MFMTs and MFs
        for tempidx = 1:nTs

            % Get all MF = 3 cross sections.
            % fprintf(fid,' 3 / \n');
            fprintf(fid,' %i %i / \n',MFMTs');
            if exist('MFs','var')
                fprintf(fid,' %i / \n',MFs');
            end

            if isfield(p.p,'farmInelastic') && ...
                    p.p.farmInelastic && ...
                    sum(Meta(m).MTs == minInelasticMT) > 0 && ...
                    Meta(m).maxInelasticMT >= minInelasticMT
                    
                iInelastic = find(Meta(m).MTs >= 51 && Meta(m).MTs <= 91);
                inelasticMTs = Meta(m).MTs(iInelastic);

                if Meta(m).ZAID == 92235
                    % issue with MF 8 being for MT 6
                    for mt = 1:length(inelasticMTs)
                        if mt == 91
                            fprintf(fid, ' 8 %i / \n', Meta(m).MTs(mt));
                        else
                            fprintf(fid, ' 6 %i / \n', Meta(m).MTs(mt));
                        end
                    end
                else
                    % we are asked to get inelastic cross sections and
                    % the ENDF tape has the first inelastic cross section.
                    fprintf(fid, ' 6 %i / \n', ...
                        iInelastic);
                    %minInelasticMT:Meta(m).maxInelasticMT;
                end
            end
            %fprintf(fid,' 3 251 / \n');
            %if Meta(m).numProtons >= 92 % shouldn't this be the other way? EDIT
            %    fprintf(fid,' 3 452 / \n');
            %end
            % Get all MF = 6 cross sections.
            % fprintf(fid, ' 6 / \n');
            fprintf(fid,' 0/ \n');
        end

        fprintf(fid,[...
            ' 0/ \n' ... % terminate THERMR
            ' moder \n' ... % MODER converts back from binary
            ' -26 27/ \n' ... % input, output
            ' stop \n' ... % terminate NJOY
            'EOF \n']); % end of file marker
        fclose(fid);

        if Meta(m).isFissionable
            fissString = '';
        else
            fissString = 'non-';
        end

        if p.verbose
            fprintf(['[FarmXS/Plant]: Script written for ' ...
                '%sfissile nuclide %s ZAID %d\n'],...
                fissString, Meta(m).nuclideStr, Meta(m).ZAID);
        end

    end
end

% bound materials EDIT.

if p.verbose
    disp('[FarmXS/Plant]: Leaving.');
end

end

function [Meta, fileStatuses]= Fertilize(p, Meta)

if p.verbose
    disp('[FarmXS/Fertilize]: Entering.');
end

% change the size of this?
nFilesPerNuclide = 4;
fileStatuses = zeros(length(Meta), nFilesPerNuclide);

disp('[FarmXS/Fertilize]: Running NJOY.');

fnames = p.fnames;

mainidxs = MapZAIDs(Meta, p.ZAIDs); % Call this every time?

for m = mainidxs

    if p.verbose
        disp(['[FarmXS/Fertilize]: Running NJOY on ' Meta(m).nuclideStr '.']);
    end

    for i = 1:length(fnames)
        Meta(m).fileStatus(i) = exist([Meta(m).storagePath filesep fnames{i}]);
    end

    if isempty(Meta(m).writeStatus)
        error(['ERROR [FarmXS/Fertilize]: Meta(m).writestatus is empty. '...
            'Must run FarmXS/Mine.']);
    end

    if Meta(m).writeStatus ~= 1 || min(Meta(m).fileStatus(1:2)) == 0
        % Files don't exist to run NJOY
        Meta(m).fileStatus(4) = -1;

        %disp(['WARNING [FarmXS/Fertilize]: ' Meta(m).nuclideStr ...
        %    ' missing the file ' ...
        %    fnames{find(Meta(m).fileStatus(1:2) == 0)} '.']);
        disp(['WARNING [FarmXS/Fertilize]: ' Meta(m).nuclideStr ...
            ' is being skipped because of a file error.']);

    elseif Meta(m).fileStatus(3) == 2 && p.f.overwriteOutput == 0

        if p.verbose
            disp(['[FarmXS/Fertilize]: ' Meta(m).nuclideStr ...
                ' output exists. Skipping file.']);
        end
    else
        % Run NJOY.
        Meta(m).fileStatus(4) = 0;

        % Prepare.
        for i = 1:2 % tape20 and script
            copyfile([Meta(m).storagePath filesep fnames{i}],fnames{i});
        end

        if p.f.logResults
            logFileStream = '>> njoy_log';
        else
            logFileStream = ' ';
        end

        % Call NJOY from the system.
        if isunix
            system(['xnjoy<script' logFileStream]);
        elseif ispc
            system(['njoy<script' logFileStream]);
        end

        % Clean up.
        delete tape*
        copyfile(fnames{3}, [Meta(m).storagePath filesep fnames{3}]);

        delete(fnames{3});

        if p.f.logResults == 1
            logFilename = logFileStream(4:end);
            copyfile(logFilename, [Meta(m).storagePath filesep ...
                logFilename]);
            delete(logFilename);
        end

    end

    fileStatuses(m,:) = Meta(m).fileStatus;

end

% "getmats" bound materials EDIT.

if p.verbose
    disp('[FarmXS/Fertilize]: Leaving.');
end

end

function [Meta, L, njoysMTs] = Harvest(p, Meta)

if p.verbose
    disp('[FarmXS/Harvest]: Entering.');
end

minInelasticMT = 51;

% See assignment of p.fnames at the top of this file.
fnames = p.fnames;

mainidxs = MapZAIDs(Meta, p.ZAIDs); % Call this every time?

% We know some files have issues. EDIT.
%Errors([300 312 313]) = [1 1 1];

% Get all the MTs used in all the output files we want to parse.
if p.h.makeLibrary
    njoysMTs = [];
    NJOYGroupDef = [];
    for m = mainidxs

        % EDIT: inefficient to call this twice, but necessary if only Harvest
        % is being called.
        [Meta(m).isFissionable, Meta(m).maxInelasticMT, Meta(m).MTs] = ...
            GetMFMTs(Meta(m).storagePath, Meta(m).writeStatus);
        fdirname = [Meta(m).storagePath filesep fnames{3}];
        fid = fopen(fdirname);
        if fid == -1
            error(['ERROR [FarmXS/Harvest]: Failed to open file ' fdirname '.']);
        end

        while ~feof(fid)
            fline = fgetl(fid);

            % Read in MFMT info and get cross section data
            if length(fline) >= 4 && strcmp(fline(2:4), 'for')
                MF = str2num(fline(10));
                MT = str2num(fline(18:20));
                if ProcessMFMT(MF, MT)
                    njoysMTs = unique([njoysMTs MT]);
                end
            elseif strcmp(p.p.IGNstr, '1') && ...
                    length(fline) >= 25 && ...
                    strcmp(fline(1:25), ' neutron group structure.')
                prospectiveGroupDef = GetNJOYGroupDef(fid);
                if isempty(NJOYGroupDef)
                    NJOYGroupDef = prospectiveGroupDef;
                elseif ~isequal(NJOYGroupDef, prospectiveGroupDef)
                    disp(['ERROR [FarmXS/Harvest]: NJOY outputs have ' ...
                        'different group definitions (conflict at ' ...
                        Meta(m).storagePath]);
                else
                end
            end
        end
        fclose(fid);
    end
    p.p.groupDef = NJOYGroupDef';
    L = initLibrary(p, Meta, njoysMTs);

end

for m = mainidxs

    fdirname = [Meta(m).storagePath filesep fnames{3}];
    fid = fopen(fdirname);
    if Meta(m).writeStatus ~= 1
        disp('[FarmXS/Harvest]: FLAG parse error 1');
        % Remove bad numbers and replace with NaN. This function is Geoff's. I
        % should incorporate it more somehow.
        %done = findreplace(fid,{' -Inf+***'},{'NaN'});
        done = findreplace(fid,{' NaN+***'},{'NaN'});
        % Close all files.
        fclose all;
        fid = fopen('temporaryfile'); % EDIT.
    end

    if fid == -1
        error(['ERROR [FarmXS/Harvest]: Failed to open file ' fdirname '.']);
    end

    %{
    % Grab the list of MTs in the NJOY output file, use this in the Library.
    if p.h.makeLibrary && isFirstFile
        njoysMTs = [];
        while ~feof(fid)
            fline = fgetl(fid);
            % Read in MFMT info and get cross section data
            if length(fline) >= 4 && strcmp(fline(2:4), 'for')
                njoysMTs = [njoysMTs str2num(fline(18:20))];
            end
        end
        njoysMTs = unique(njoysMTs);
        frewind(fid);
    end
    %}

    while ~feof(fid)
        fline = fgetl(fid);

        if ~isempty( regexp(fline, '.***error in findf***.', 'once'))
            % Crete a parse error log file.
            ferrorid = fopen('PARSE_ERRORS', 'a'); % rename EDIT.
            fprintf(ferrorid, 'NJOY output file %s is corrupted.\n',...
                Meta(m).nuclideStr);
            fprintf(ferrorid, '  %s\n', fline);
            fprintf(ferrorid, '  Data$ %d %s %s\n\n',...
                m, Meta(m).nuclideStr, Meta(m).storagePath);
            fclose(ferrorid);
            error('Errror in findf');
            break; % EDIT.
        end

        % maybe we loop through the file first to gather the list of MFMTs
        if length(fline) >= 37 % 30? What is the ground for this number? EDIT.
            % Only lines that have over a length of 37 are useful to us? What
            % if we only ask for 1 temperature? EDIT!

            % Pull data from the file.

            if strcmp(fline(21:22), 't=')
                % We can read in the current temperature.
                temp = str2num(fline(23:31));
                tempStr = int2str(temp);

            elseif strcmp(fline(2:4), 'for')

                % Read in MFMT info and get cross section data
                MFstr = fline(10);
                MF = str2num(MFstr);
                MT = str2num(fline(18:20));
                MTstr = int2str(MT);
                MFMTstr = [MFstr ' ' MTstr];

                % LIBRARY
                if p.h.makeLibrary
                    ZAID = Meta(m).ZAID;
                    L.z(L.ZAID(ZAID)).isFissionable = Meta(m).isFissionable;
                    if sum(Meta(m).MTs == minInelasticMT) == 1
                        iInelastic = ...
                            find(Meta(m).MTs >= 51 & Meta(m).MTs <= 91);
                        L.z(L.ZAID(ZAID)).inelasticMTs = ...
                            Meta(m).MTs(iInelastic);
                    else
                        L.z(L.ZAID(ZAID)).inelasticMTs = [];
                    end
                end

                % ARRAY
                arrayStr = [Meta(m).element int2str(Meta(m).atomNum) ...
                    '_' MFstr '_' MTstr '_' tempStr];

                if p.verbose
                    fprintf('[FarmXS/Harvest]: Harvesting %s\n',...
                        arrayStr);
                end

                if ProcessMFMT(MF, MT)
                    [XS, S0s] = GetNJOYXS(fid);

                    logS0s = log10(S0s);
                    nLogS0s = length(logS0s);

                    logS0strs = cell(1,nLogS0s);

                    for iS0 = 1:nLogS0s
                        if logS0s(iS0) < 0
                            logS0strs{iS0} = ['_m' num2str(-logS0s(iS0))];
                        else
                            logS0strs{iS0} = ['_' num2str(logS0s(iS0))];
                        end
                    end
                    if strcmp(MFMTstr, '3 1')
                        L = Harvest31(p, Meta, m, L, ZAID, MT, temp, logS0s,...
                            nLogS0s, XS);
                    elseif strcmp(MFstr, '3')
                        L = HarvestMF3(p, Meta, m, L, ZAID, MT, temp, ...
                            logS0s, nLogS0s, XS);
                    elseif strcmp(MFstr, '6') || strcmp(MFMTstr, '8 91')
                        L = HarvestMF6(p, Meta, m, L, ZAID, MT, temp, ...
                            logS0s, nLogS0s, XS);
                    else
                        disp(['ERROR [FarmXS/Harvest]: ' ...
                            Meta(m).nuclideStr ...
                            ' MF ' MFstr ' MT ' MTstr ...
                            ': is not being harvested. There is no case '...
                            'for it.']);
                    end
                else
                    disp('[FarmXS/Harvest]: Skipping MFMT, don''t want it.');
                end % if ProcessMFMT

                % Call "datatypes"
                %
                % load CurrentMFMT
                %
                % COMPARE to correct for thermal region (old code)
            end


        end % if strcmp(fline(1:37)...
        % end % if length(fline)...
    end % while ~feof(fid)

    fclose all;

end % for m = mainidxs

if isfield(p.h,'arrayMatName');
    save(p.h.arrayMatName);
end

if p.verbose
    disp('[FarmXS/Harvest]: Leaving.');
end

end % function Harvest

function [Meta, L] = CombineInelastic(p, Meta, L)
if p.MFs ~= 6
    warning('MF 6 needs to be requested to combine inelastic.');
end
for z = L.ZAIDs
    %if z ~= 222
        for m = L.z(L.ZAID(z)).inelasticMTs
          %  if m ~= 91
                L.z(L.ZAID(z)).m(L.MT(6)).xs = ...
                    L.z(L.ZAID(z)).m(L.MT(6)).xs + ...
                    L.z(L.ZAID(z)).m(L.MT(m)).xs;
          %  end
        end
%        L.z(L.ZAID(z)).m(L.MT(2)).xs = L.z(L.ZAID(z)).m(L.MT(2)).xs + ...
%            L.z(L.ZAID(z)).m(L.MT(6)).xs;
    %end
end
end % function CombineScattering

function [Meta, L] = MakeWater(p, Meta, L)
%{

i think i need to create the inelastic kernel here.
for m = L.MTs % NOT SCATTERING. also, only want to do this for specified ZAIDs,
    namely 102 and 251
    xshbound = L.z(L.MT(1001)).m(L.MT(m)).xs;
    xshbound(1:nr,1:nc) = L.z(L.MT(1001)).m(L.MT(??)).xs(1:nr,1:nc);
    L.z(L.ZAID(222)).m(L.MT(m)).xs = 2*xshbound + ...
    L.z(L.ZAID(8016)).m(L.MT(m)).xs;
end


%}
end % function HarvestWater

function truefalse = ProcessMFMT(MF, MT)
% Determines if this is an MFMT to be processed

if MF == 8 && MT ~= 91
    truefalse = false;
    return;
elseif (MF == 6 && (MT ~= 2 && MT < 51))
    truefalse = false;
    return;
elseif (MT > 9 && MT < 16) || (MT > 26 && MT < 38) || (MT > 38 && MT < 51)
    truefalse = false;
    return;
else
    truefalse = true;
end
end % function ProcessMFMT

% Harvest31 (MFMT 3 1)
function L = Harvest31(p, Meta, m, L, ZAID, MT, temp, logS0s, nLogS0s, XS);
% Lengedre order 0
XSLORD0 = XS(1:4:end,:);
% Legendre order 1
XSLORD1 = XS(3:4:end,:);
% Flux Legendre order 0
FluxLORD0 = XS(2:4:end,:);
% Flux Legendre order 1
FluxLORD1 = XS(4:4:end,:);

for iS0 = 1:nLogS0s
    if p.h.makeArrays
        eval([arrayStr logS0strs{iS0} ...
            ' = XSLORD0(:, iS0 + 2);']);
        eval([arrayStr logS0strs{iS0} '_L1' ...
            ' = XSLORD1(:, iS0 + 2);']);
        eval([arrayStr logS0strs{iS0} '_F0' ...
            ' = FluxLORD0(:, iS0 + 2);']);
        eval([arrayStr logS0strs{iS0} '_F1' ...
            ' = FluxLORD1(:, iS0 + 2);']);
    end

    if p.h.makeLibrary
        %dbstop in FarmXS.m at 708t
        L.z(L.ZAID(ZAID)).m(L.MT(MT)).xs(:,L.T(temp),L.S0(logS0s(iS0))) = ...
            XSLORD0(:, iS0 + 2);
    end
end
clear XSLORD0 XSLORD1 FluxLORD0 FluxLORD1;
end % function Harvest31

% HarvestMF3
function L = HarvestMF3(p, Meta, m, L, ZAID, MT, temp, logS0s, nLogS0s, XS);
% EDIT this has not been generalized.
if nLogS0s ~= length(p.p.S0s)
    disp(['WARNING [FarmXS/Harvest]: ' ...
        'S0 mismatch. ' ...
        'Copying output for all S0.']);
   % Using first set of NJOY '...
   %     'data for requested S0 outputs.']);
        %'number of S0s for which this MFMT has'...
        %' cross sections is not the same as '...
        %'the number of S0s requested. '...
    for iS0 = 1:length(p.p.S0s)
        if p.h.makeArrays
            eval([arrayStr logS0strs{iS0} ...
                ' = XS(:, 2);']);
        end
        if  p.h.makeLibrary
            L.z(L.ZAID(ZAID)).m(L.MT(MT)).xs(XS(:,1),L.T(temp),iS0) = XS(:, 2);
        end
    end
else
    for iS0 = 1:nLogS0s
        if p.h.makeArrays
            eval([arrayStr logS0strs{iS0} ...
                ' = XS(:, iS0 + 1);']);
        end
        if  p.h.makeLibrary
            L.z(L.ZAID(ZAID)).m(L.MT(MT)).xs(XS(:,1),L.T(temp),L.S0(logS0s(iS0))) = ...
                XS(:, iS0 + 1);
        end
    end
end
end % function HarvestMF3

%% HarvestMF6
function L = HarvestMF6(p, Meta, m, L, ZAID, MT, temp, logS0s, nLogS0s, XS)
% Must manage Legendre order!

if XS(1,3) ~= 0
    disp(['ERROR [FarmXS/Harvest]: ' ...
        Meta(m).nuclideStr ...
        ' MF 6 MT 2: should use XS(i, j + 2)']);
end

% Initialize for the upcoming loop.
[lengthXS junkvariable] = size(XS);

kernel = zeros(L.nGroups);
for iS0 = 1:nLogS0s
    for iXS = 1:lengthXS
        % kernel(i -> j)
        % XS's first and second columns are bin indices.
        % Third column is Legendre order.
        kernel(XS(iXS,2), XS(iXS,1)) = XS(iXS, iS0 + 3);
    end

    if p.h.makeArrays
        eval([arrayStr logS0strs{iS0} ...
            ' = kernel;']);
    end
    if p.h.makeLibrary
        L.z(L.ZAID(ZAID)).m(L.MT(MT)).xs(:,:,L.T(temp),L.S0(logS0s(iS0))) = ...
            kernel;
    end
end
clear kernel;
end % function HarvestMF6

function [NJOYGroupDef] = GetNJOYGroupDef(fid)

binidx = 0;
fline = fgetl(fid);
while ~isempty(fline)
    binidx = binidx + 1;
    leftGroupBound(binidx) = str2double(fline(10:20));
    rightGroupBound(binidx) = str2double(fline(26:36));

    fline = fgetl(fid);
end

NJOYGroupDef = [leftGroupBound(1) rightGroupBound]';

end


function [XS, S0s] = GetNJOYXS(fid)
% cannot yet handle the kernel format for if there is only one s0/infinite
% dilution.

XS = [];
% Move down enough lines to get to the header.
%while 1
%    fline = fgetl(fid)
%    ~isempty( regexp(fline, '.infinit.', 'once'))
%    if ~isempty( regexp(fline, '.infinit.', 'once')) || ...
%            ~isempty( regexp(fline, '.group  group.', 'once'))
%        % EDIT. I do not understand these conditions.
%        % disp('[FarmXS/GetNJOYXS]: Skipping a cross section. Investigate.');
%        break;
%    end
%end
fline = fgetl(fid);
while isempty( regexp(fline, '.infinit.', 'once')) && ...
        isempty( regexp(fline, '.group  group.', 'once'))
    fline = fgetl(fid);
end


% Get header of the XS table for this MFMT. Put in E's into the header so the
% header can be converted to numbers.
fline = regexprep(regexprep(regexprep(fline, '+', 'E+'), '-', 'E-'),' E-',' -');
% Index of where S0 numbers begin in the header.
firstindex = regexp(fline, '[\dE.+-]+', 'once');
% This next line assumes we are always pulling infinite dilution.
S0s = [1e10, str2num( fline(firstindex:end))];
% Move to the next (empty) line.
fgetl(fid);

% We do not know beforehand the number of rows for XS, unless we preparsed or
% something. EDIT.

% This line is needed to evaluate the regular expression eval below. The
% variable flux appears amongst the total cross section data. I am not sure why
% we don't replace this with NaN. EDIT.
flux = NaN;
flx = NaN;

% Read the first line of data.
fline = fgetl(fid);
lineidx = 0;
while isempty(strfind(fline, 'xx ans')) ...
        || ((length(strfind(fline, 'xx ans')) == 1) ...
        && (strfind(fline, 'xx ans') ~= 59))
    %length(fline) >= 11 % && ~strcmp(fline(11), ' ')
    % I don't know where the 11 condition comes from.

    %fline = regexprep(fline,'NaN\+\*\*\*','NaN');
    fline = regexprep(fline,'NaN\+\*\*\*','0');

    if (length(fline) < 7) ...
            || (length(strfind(fline, '---message from get')) == 1 ...
            && (strfind(fline, '---message from get') == 2))
        % takes care of (ignores) warning messages from getunr and getsed
    else
        lineidx = lineidx + 1;
        eval(['XS(lineidx,:) = [' regexprep( regexprep( regexprep( fline, ...
            '+', 'E+'), '-', 'E-'), ' E-', ' -') '];']);
    end

    fline = fgetl(fid);

end

%    XS(lineidx,:) = str2num(regexprep( regexprep( regexprep( fline, '+', 'E+'), '-', 'E-'), ' E-', ' -');
%{
    eval(['tempp = [' regexprep( regexprep( regexprep( fline, ...
            '+', 'E+'), '-', 'E-'), ' E-', ' -') '];']);
    fline
    XS(lineidx,:)
    tempp
%}

% There is a bunch of extra code here to deal with mismatched vectors, for
% stuff like MFMT = '6 221'. EDIT.

end % GetNJOYXS


function [isFissionable, maxInelasticMT, MTs] = GetMFMTs(storagePath, writeStatus)
% Geoff's scripts have some valuable comments for here.
% EDIT this script is selective, will censor the available MTs.

if writeStatus == 0
    isFissionable = 0;
    maxInelasticMT = 0;
    disp(['WARNING [FarmXS/GetMFMTs]: No tape20 at ' storagePath ...
        '. Check input options.']);

else

    [fid message] = fopen(fullfile(storagePath,'tape20'));
    
    if fid == -1
        error([message ': ' storagePath]);
    end

    % Move to the second line of the file.
    fgetl(fid);
    fline = fgetl(fid);

    % Obtain ZAID from second line of file. EDIT should compare to the ZAID
    % grabbed in Explore.
    % ZAID = str2num([fline(2:9) 'E' fline(10:11)]);
    % atomMass = str2num([fline(13:20) 'E' fline(21:22)]);

    % Store all the MTs that are found.
    isFissionable = 0;
    MTs = [];
    has91 = 0;
    runLoop = 1;
    while ~feof(fid) && runLoop
        fline = fgetl(fid);
        try 
            if ~strcmp( fline(73:75), '451') 
                % 451 appears to the right side of the ENDFB tape for all of the
                % front matter of the file.
    
                while runLoop
                    fline = fgetl(fid);
    
                    MF = str2num(fline(72));
                    MT = str2num(fline(73:75));
    
                    if MT == 18
                        isFissionable = 1;
                    end
    
                    % Only consider scattering.
                    if MF == 3
                        if MT >= 91
                            if MT == 91
                                has91 = 1;
                            end
                            % Looked far enough.
                            runLoop = 0;
                        elseif ~sum(MT == MTs)
                            % If we have this MT already, move on. If not,
                            % collect it.
                            MTs = [MTs, MT];
                        end
                    elseif length(MTs) > 2
                        runLoop = 0;
                        disp(['[FarmXS/GetMFMTs]: Not many MTs found,'...
                            ' break from MF ~= 3 for ' storagePath]);
                    end
                end
            end
        catch
            disp('          IS YOUR INTERNET WORKING?');
            rethrow;
        end
    end

    fclose(fid);

    maxInelasticMT = max(MTs);
    %disp('[FarmXS/GetMFMTs]: Not including MT 91');
    if has91
        MTs = [MTs 91];
    end

end

end

function matidxs = MapZAIDs(Meta, ZAIDs)


nZAIDs = length(ZAIDs);
AllZAIDs = zeros(nZAIDs,1);

for m = 1:length(Meta)
    AllZAIDs(m) = Meta(m).ZAID;
    % check for uniqueness of ZAID inputs?
end

ZAIDsENDFDoesNotHave = [];

% EDIT manage AllZAIDs having duplicates.
for z = 1:nZAIDs
    matidx = find(AllZAIDs == ZAIDs(z));
    if isempty(matidx)
        ZAIDsENDFDoesNotHave = [ZAIDsENDFDoesNotHave ZAIDs(z)];
    else
        % EDIT not okay for isomerics!
        matidxs(z) = matidx(1);
    end
end

if ~isempty(ZAIDsENDFDoesNotHave)
    error(['[FarmXS/MapZAIDs]: ENDF does not have a tape20 for the ' ...
        'requested ZAID(s): %i \n'], ZAIDsENDFDoesNotHave);
end

end

function L = initLibraryMeta(p, Meta, njoysMTs)

%[height width] = size(p.p.MFMTs);
%for i = 1:height
%    if p.p.MFMTs(i,2) < 0 % is negative
%    end
%
%end

ZAIDs = sort([p.ZAIDs]); %[1001 8016 11023 92235 92238 222];
%MTs = [p.p.MFMTs(:,2)' 6 7 8 9]; %[1 2 4 6 7 8 9 18 102 221 251 452]; T120103
MTs = [6 7 8 9 sort(unique([njoysMTs p.p.MFMTs(:,2)']))]; %[1 2 4 6 7 8 9 18 102 221 251 452];
%MTs = [1 2 4 6:9 18 51:91 102 221 251 452];
mainMTs = [1 2 6 7 8 9 18 102 251 452];
Ts = p.p.Ts; %[300 600 900 1200 1500];
S0s = p.p.S0s;

L = struct('groupDef',p.p.groupDef,...
    'nGroups',length(p.p.groupDef)-1,...
    'ZAIDs',ZAIDs,...
    'MTs',MTs,...
    'mainMTs',mainMTs,...
    'Ts',Ts,...
    'S0s',log10(S0s),...
    'ZAID',[],...
    'MT',[],...
    'mainMT',[],...
    'T',[],...
    'S0',[],...
    'MTmap',[],...
    'z',struct('isFissionable',[],...
    'inelasticMTs',[],...
    'm',struct('hasResonances',[],...
    'xs',[])));
%           'MTnames', cell(1,length(MTs)),...

L.ZAID(ZAIDs) = 1:length(ZAIDs);
L.MT(MTs) = 1:length(MTs);
L.mainMT(MTs) = 1:length(MTs);
L.T(Ts) = 1:length(Ts);
keys = num2cell(log10(S0s));
values = num2cell(1:length(keys));
L.S0 = containers.Map(keys, values);

L.ZAID = sparse(L.ZAID);
L.MT = sparse(L.MT);
L.mainMT = sparse(L.mainMT);
L.T = sparse(L.T);
nTs = length(Ts);
nS0s = length(S0s);

MTcellmap1 = {1;
    2;
    4;
    6;
    7;
    8;
    9;
    18;
    102;
    251;
    452};
MTcellmap2 = {'total';
    'elastic';
    'inelastic';
    'totalinelastic';
    'vbudsiitotal';
    'vbudsiitransport';
    'vbudsiinufission';
    'fission';
    'capture';
    'mubar';
    'nubar'};

L.MTmap = containers.Map(MTcellmap1, MTcellmap2);
L.MTmapback = containers.Map(MTcellmap2, MTcellmap1);
end

function L = initLibrary(p, Meta, njoysMTs)

L = initLibraryMeta(p, Meta, njoysMTs);
nTs = length(L.Ts);
nS0s = length(L.S0s);

for i = L.ZAIDs
    L.z(L.ZAID(i)).isFissionable = 0;
    for j = L.MTs
        L.z(L.ZAID(i)).m(L.MT(j)).hasResonances = 1;
        if j == 2  || j == 6 || j == 221 || (j >= 51 && j <= 91) % VERY DANGEROUS. EDIT.
            L.z(L.ZAID(i)).m(L.MT(j)).xs = ...
                zeros(L.nGroups,L.nGroups,nTs,nS0s);
        %elseif j == 452
        %    % THIS IS TEMPORARY, SINCE ANDY DOES NOT STORE 452
        %    L.z(L.ZAID(i)).m(L.MT(j)).xs = ones(L.nGroups,nTs,nS0s);
        else
            L.z(L.ZAID(i)).m(L.MT(j)).xs = zeros(L.nGroups,nTs,nS0s);
        end
    end
end

%{
for mtidx = 1:length(MTs)
    if isKey(MTmap, MTs(mtidx))
        L.MTnames{mtidx} = MTmap(MTs{mtidx});
    else
        L.MTnames{mtidx} = 'MT name unknown';
    end
end
%}

end


% TRASH
%% MFMT for all files.
%        MFMTall = [3 1;
%                   6 2;
%                   3 102;
%                   3 251];
%                 % 3 2;
%
%        % MFMT at high temp. EDIT should this be max(nGroups) or max(groupDef)?
%        if max(p.p.groupDef) < 1000
%            MFMThightemp = [];
%        else
%            MFMThightemp = [3 4];
%        end
%
%        % MFMT at low temp only. EDIT same as previous
%        if min(p.p.groupDef) > 10
%            MFMTlowtemp = [];
%        else
%            MFMTlowtemp = [3 221;
%                           6 221;
%                           3 222;
%                           6 222];
%        end
%
%        % MFMT related to fission.
%        if Meta(m).isFissionable
%            MFMTfiss = [3 18;
%                        3 452];
%                      % 6 18;
%        else
%            MFMTfiss = [];
%        end
%
%        % Combine everything
%        MFMTs = [MFMTall;
%                 MFMThightemp;
%                 MFMTlowtemp;
%                 MFMTfiss];




%{
                    case {'1000'} % {'6 51', '6 221'} % EDIT NOT WORKING 221.
                        [nrows, ncols] = size(XS);

                        for rowidx = 1:nrows
                            for colidx = 0:ncols - 3
                                if XS(rowidx, colidx + 3) ~= -20
                                kernel( XS(rowidx,2)+colidx, XS(rowidx,1)) ...
                                        = XS(rowidx, colidx + 3);
                                end
                            end
                        end

                        for iS0 = 1:nLogS0s
                            if p.h.makeArrays
                                eval([arrayStr logS0strs{iS0} ...
                                    ' = kernel;']);
                            end
                            if p.h.makeLibrary
    L.z(L.ZAID(ZAID)).m(L.MT(MT)).xs(:,:,L.T(temp),L.S0(logS0s(iS0))) = ...
        kernel;
                            end

                        end
                        clear kernel;
%}
% I think something is incorrect with how this loop
% works.
% EDIT.
%                                if ...
%    size(L.z(L.ZAID(ZAID)).    m(L.MT(MT)).xs(:,L.T(temp),L.S0(logS0s(iS0)))) ...
%        ~= size(XS(:, iS0     + 1))
%                disp(['ERRO    R [FarmXS/Harvest]: ' ...
%                    Meta(m)    .nuclideStr ...
%                    ' MF '     MFstr ' MT ' MTstr ': size error, not harvested.']);
%                                    else
