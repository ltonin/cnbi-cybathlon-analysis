function findRehearsalStyleRaces(SubID, GDFSesPath, TargetPath)

if nargin == 2
    TargetPath = [GDFSesPath '/' SubID '_RaceMat/'];
end
TargetPath = regexprep(TargetPath, '//+', '/');

% Add this
addpath(genpath('/home/sperdikis/Git/cnbi-cybathlon-analysis/'));

% Add biosig
addpath(genpath('~/Git/cnbi-smrtrain/toolboxes/biosig/'));

% Find existing session folders
SesFolders = dir([GDFSesPath '/' SubID '_*']);

for ses=1:length(SesFolders)
    SessionRaceCounter = 0;
    RunDate = {};
    RunTime = [];
    SessionDate = [];
    Runs = [];
    
    % Find rehearsal-style race runs (*.cybathlon.*)
    Runs = dir([GDFSesPath '/' SesFolders(ses).name '/' SubID '.cybathlon.*.gdf']);
    if(isempty(Runs))
        continue;
    end
    
    % Organize in temporal order
    for run=1:length(Runs)
        RunDots = strfind(Runs(run).name,'.');
        RunDate{run} = Runs(run).name(RunDots(2)+1:RunDots(3)-1);
        RunTime(run) = datenum(datetime(Runs(run).name(RunDots(3)+1:RunDots(4)-1),'InputFormat','HHmmss'));
    end
    if(length(unique(RunDate))>1)
        disp(['Runs in session folder have different dates, exiting...']);
        return;
    else
        SessionDate = unique(RunDate);
    end
    % Reorganize ascending
    [mpla sortedInd] = sort(RunTime,'ascend');
    Runs = Runs(sortedInd);
    
    for run=1:length(Runs)
        try
            [data header] = sload([GDFSesPath '/' SesFolders(ses).name '/' Runs(run).name]);
          
        catch
            disp(['Skipping run ' Runs(run).name ' , GDF file is corrupted!']);
            continue;
        end
        
        % Check if there are triggers
        TriggersInd = find(ismember(header.EVENT.TYP,[783 769 770 771 772 773 774 775]));
        
        % Remove triggers that are essentially the same one (group triggers)
        remPOS = [];
        for i=2:length(TriggersInd)
            if ((header.EVENT.POS(TriggersInd(i)) - header.EVENT.POS(TriggersInd(i-1)))/512 < 0.1)
                remPOS = [remPOS header.EVENT.POS(TriggersInd(i))];
            end
        end
        remInd = find(ismember(header.EVENT.POS,remPOS));
        header.EVENT.TYP(remInd) = [];
        header.EVENT.POS(remInd) = [];
        header.EVENT.DUR(remInd) = [];
        
        % Remove weird triggers (for online modification of thresholds etc)
        remInd = find(ismember(header.EVENT.TYP,[0 101 102 103 104]));
        header.EVENT.TYP(remInd) = [];
        header.EVENT.POS(remInd) = [];
        header.EVENT.DUR(remInd) = [];
        
        
        % Recalculate existing trigger positions
        TriggersInd = find(ismember(header.EVENT.TYP,[783 769 770 771 772 773 774 775]));
        NTriggers = length(TriggersInd);
        
        if(NTriggers < 18)
            % Skip run
            disp(['Skipping run ' Runs(run).name ' , there are not enough triggers, incomplete race!']);
            continue;
        end
        
        % Find the first and last 783 to be sure we have a full race in
        % between
        Ind783 = find(header.EVENT.TYP(TriggersInd)==783);
        First783 = Ind783(1);
        Last783 = Ind783(end);
        NPadsAvailable = Last783-First783-1;
        NFakeRaces = floor(NPadsAvailable/16);
        
        FirstPadTYP = header.EVENT.TYP(TriggersInd(First783):TriggersInd(First783+1));
        FirstPadPOS = header.EVENT.POS(TriggersInd(First783):TriggersInd(First783+1)) - header.EVENT.POS(TriggersInd(First783)) + 1;
        FirstOffset = FirstPadPOS(end)-1;
        FirstPadTYP = FirstPadTYP(1:end-1);
        FirstPadPOS = FirstPadPOS(1:end-1);
        FirstPadData = data(header.EVENT.POS(TriggersInd(First783)):header.EVENT.POS(TriggersInd(First783))+FirstOffset,:);
        
        OneSecondBefore = data(header.EVENT.POS(TriggersInd(First783))-512:header.EVENT.POS(TriggersInd(First783))-1,:);
        
        LastPadTYP = header.EVENT.TYP(TriggersInd(Last783):end);
        LastPadPOS = header.EVENT.POS(TriggersInd(Last783):end) - header.EVENT.POS(TriggersInd(Last783)) + 1;
        LastPadData = data(header.EVENT.POS(TriggersInd(Last783)):end,:);
        for r=1:NFakeRaces
            MiddlePadTYP = header.EVENT.TYP(TriggersInd(First783 + (r-1)*16 + 1) : TriggersInd(First783 + (r-1)*16 + 17));
            MiddlePadPOS = header.EVENT.POS(TriggersInd(First783 + (r-1)*16 + 1) : TriggersInd(First783 + (r-1)*16 + 17)) ...
                - header.EVENT.POS(TriggersInd(First783 + (r-1)*16 + 1)) ...
                + FirstOffset+1;
            MiddleOffset = MiddlePadPOS(end)-1;
            MiddlePadTYP = MiddlePadTYP(1:end-1);
            MiddlePadPOS = MiddlePadPOS(1:end-1);
            MiddlePadData = data(header.EVENT.POS(TriggersInd(First783 + (r-1)*16 + 1)) : header.EVENT.POS(TriggersInd(First783 + (r-1)*16 + 17))-1,:);
            
            newTYP = [FirstPadTYP ; MiddlePadTYP ; LastPadTYP];
            newPOS = [FirstPadPOS ; MiddlePadPOS ; LastPadPOS+MiddleOffset];
            newdata = [FirstPadData ; MiddlePadData ; LastPadData];
            
            % Find race end
            Ind783Last = find(newTYP == 783);
            Ind783Last = Ind783Last(end);
            CommandsAfterPad = (newPOS(Ind783Last+1:end) - newPOS(Ind783Last))/512;
            TimeOnLastPad = decideRaceEnd(CommandsAfterPad)
            POS666 = newPOS(Ind783Last) + round(TimeOnLastPad*512);
            newPOS = [newPOS ; POS666];
            newTYP = [newTYP ; 666];
            [newPOS sortInd] = sort(newPOS,'ascend');
            newTYP = newTYP(sortInd);
            Ind666 = find(newTYP==666);
            newPOS(Ind666+1:end) = [];
            newTYP(Ind666+1:end) = [];
            newdata = newdata(1:newPOS(end)+512,:); % Keep also 1 sec after 666
            
            % Now, add also 1 sec before the first trigger and transpose
            % the POS as needed
            newdata = [OneSecondBefore ; newdata];
            newPOS = newPOS + 512;
            
            % Check if there are enough commands
            CommandsInd = find(ismember(newTYP,[10 20 30]));
            if(length(CommandsInd) < 5)
                % Skip race
                disp(['Skipping race ' Runs(run).name ' , there are not enough commands!']);
            end
            
            % Prepare Race structure
            Race.data = newdata;
            Race.EVENT.POS = newPOS;
            Race.EVENT.TYP = newTYP;
            Race.EVENT.DUR = zeros(length(Race.EVENT.TYP),1);
            
            if(sum(ismember(Race.EVENT.TYP,770)))
                Race.EVENT.TYP(Race.EVENT.TYP==10) = hex2dec('6000') + 770; % Hands = 10 = 25346 (0x0302 + 0x6000)
            elseif(sum(ismember(Race.EVENT.TYP,773)))
                Race.EVENT.TYP(Race.EVENT.TYP==10) = hex2dec('6000') + 773; % Both Hands = 10 = 25349 (0x0305 + 0x6000)
            else
                disp(['Not sure which command goes right!']);
                continue;
            end
            
            Race.EVENT.TYP(Race.EVENT.TYP==20) = hex2dec('6000') + 771; % Both Feet = 20 = 25347 (0x0303 + 0x6000)
            Race.EVENT.TYP(Race.EVENT.TYP==30) = hex2dec('6000') + 769; % Left Hand = 30 = 25345 (0x0301 + 0x6000)
            
            Race.protocol = 'rehearsalstyle';
            Race.AlignmentOffset = NaN;
            Race.RaceTime = Race.EVENT.POS(end)/512;
            Race.gamecommands.type = NaN;
            Race.gamecommands.pos = NaN;
            Race.Date = SessionDate;
            Race.StartTime = NaN;
            Race.resfile = NaN;
            Race.GDFFile = Runs(run).name;
            Race.PlayerIndex = 1;
            Race.GameIndexAll = NaN;
            Race.GameIndexFound = NaN;
            RaceTriggersInd = find(ismember(Race.EVENT.TYP,[783 769 770 771 772 773 774 775]));
            lvlInt = Race.EVENT.TYP(RaceTriggersInd);
            Race.Level = '';
            for l=1:length(lvlInt)
                switch(lvlInt(l))
                    case 783
                        Race.Level = [Race.Level '0'];
                    case 769
                        Race.Level = [Race.Level '3'];
                    case 770
                        Race.Level = [Race.Level '1'];
                    case 771
                        Race.Level = [Race.Level '2'];
                    case 773
                        Race.Level = [Race.Level '1'];
                    otherwise
                        disp('Weird level, skipping');
                        continue;
                end
            end
            Race.GameStartInd = NaN;
            Race.GameEndInd = NaN;
            
            if(isequal(ismember([769 770 771 773], lvlInt),[0 1 1 0]))
                Taskset = 'rhbf';
            elseif(isequal(ismember([769 770 771 773], lvlInt),[1 1 0 0]))
                Taskset = 'rhlh';
            elseif(isequal(ismember([769 770 771 773], lvlInt),[1 1 1 0]))
                Taskset = 'rhlhbf';
            elseif(isequal(ismember([769 770 771 773], lvlInt),[0 0 1 1]))
                Taskset = 'bhbf';
            else
                disp('Weird Taskset, skipping!');
                continue;
            end
            
            % Add extra GTYP for compliancy
            Race.EVENT.GTYP = Race.EVENT.TYP;
            
            % Just check that the fake race does not have any huge gap, if
            % so, di
            if(max(diff(Race.EVENT.POS(ismember(Race.EVENT.TYP,[768 769 770 771 773 783 666]))/512)) > 20.0)
                % Skipping this race 
                disp('This fake race has pads with larger duration than it should, skipping!');
                continue;
            end
            
            % Save the Race output
            SessionRaceCounter = SessionRaceCounter + 1;

            % Creating folder if does not exist
            % 20161220 - ltonin
            [~, savepath] = cnbiutil_mkdir(TargetPath);
            tmpName = Race.GDFFile(1:end-4);
            IndOnline = strfind(tmpName,'cybathlon');
            tmpName(IndOnline:IndOnline+9)=[];
            tmpName = tmpName(1:22);
            tmpName = [tmpName '.race.mi.mi_' Taskset];
            save([savepath '/' tmpName '.race' num2str(SessionRaceCounter)  '.mat'  ],'Race');
            clear Race newPOS newTYP newdata
            fclose all;            
        end
        
    end
end