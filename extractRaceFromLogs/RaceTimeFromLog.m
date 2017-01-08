function RaceTimeFromLogOnly(SubID, DesiredSessionDate, player, LogFilePath, TargetPath)

if nargin == 3
    TargetPath = [GDFSesPath '/' SubID '_RaceMat/'];
end
TargetPath = regexprep(TargetPath, '//+', '/');

% Load results csv files
ResFiles = dir([LogFilePath  '*.csv']);
ResFiles(find(strcmp({ResFiles.name},'gameevents.csv')))=[]; % Get rid of game events
[ResDate ResTime] = cellfun(@splitdatetime,{ResFiles.name},'UniformOutput',false);
ResDateTime = strcat(ResDate,ResTime); 

% Load the huge log csv file
fid = fopen([LogFilePath '/gameevents.csv']);
LogData = textscan(fid,'%s%s%s%s%s%s%s%s%s%s%s%s',Inf,'Delimiter',',');
fclose(fid);

% Compute useful fields
[RaceDate RaceTime] = cellfun(@splitdatetime,LogData{1},'UniformOutput',false);
RaceCommandTime = cell2mat(cellfun(@gettimes,LogData{2},'UniformOutput',false));
Commands = cellfun(@trimquotes,LogData{6},'UniformOutput',false);
PlayerTrack = cellfun(@trimquotes,LogData{3},'UniformOutput',false);
PlayerName = cellfun(@trimquotes,LogData{4},'UniformOutput',false);
Trigger = cellfun(@trimquotes,LogData{5},'UniformOutput',false);
Pad = cellfun(@extractpad,LogData{11},'UniformOutput',false);
End = cellfun(@trimquotes,LogData{7},'UniformOutput',false);
LevelInd = find(strcmp(Commands,'loadlevel'));
Levels = End(LevelInd);
Levels = cellfun(@processlevels,Levels,'UniformOutput',false);


% Find possible games

% This is if we consider everything between app/startgame and any next app (game closed, reset, etc)
GameStartInd = find(strcmp(Commands,'startgame'));
%GameEndInd = setdiff(find(strcmp(Trigger,'app')),GameStartInd);
% Actually, only consider the End (app, state, End)
GameEndInd = find(strcmp(End,'End'));
StartEndInd = union(GameStartInd,GameEndInd);
StartEndIndVal = zeros(1,length(StartEndInd));
StartEndIndVal(find(ismember(StartEndInd,GameStartInd)))=1;
StartEndIndVal(find(ismember(StartEndInd,GameEndInd)))=2;
GameStartInd  = StartEndInd(strfind(StartEndIndVal,[1 2]));
GameEndInd  = StartEndInd(strfind(StartEndIndVal,[1 2])+1);

GamesFound = 0;
PrevSessionDate = '';
for g=1:length(GameStartInd)

    disp(['Checking game ' num2str(g) '/' num2str(length(GameStartInd))]);
    
    % Find levels for this game
    thisGLevel = Levels(max(find(LevelInd-GameStartInd(g) < 0)));
    
    % Find players in this game
    GPlayers = unique(PlayerTrack(GameStartInd(g):GameEndInd(g)));
    % The first one is always the empty string, remove it
    GPlayers(1)=[];
    if(isempty(GPlayers))
        disp(['This game has no players! Skipping'])
        continue; % This game has no players, never started!
    end
    
    thisGDate = RaceDate(GameStartInd(g));
    thisGTime = RaceTime(GameStartInd(g));
    
    if(strcmp(thisGDate{1},DesiredSessionDate)==0)
        disp(['This race was not ran on the desired date! Skipping'])
        continue; % This game has no players, never started!        
    end
    
    % Reset session race counter if this is a new session
    if(strcmp(thisGDate{1},PrevSessionDate)==0)
        SessionRaceCounter = 0;
    end
    PrevSessionDate = thisGDate{1};

       % Find race end from game result csv
    thisResFile = find(datenum(datetime(ResDateTime,'InputFormat','yyyyMMddHHmmss')) -...
    datenum(datetime(strcat(thisGDate{1},RaceTime(GameEndInd(g))),'InputFormat','yyyyMMddHHmmss'))<0);

    if(isempty(thisResFile))
        thisResFile = 1;
    else
        thisResFile = thisResFile(end) + 1;    
    end
    thisResFile = ['result' ResDate{thisResFile} ResTime{thisResFile} '.csv'];
    fid2 = fopen([LogFilePath '/' thisResFile]);
    ResData = textscan(fid2,'%s%s%s%s%s%s%s%s',Inf,'Delimiter',',');

    % Find the player row
    PlInd = find(strcmp(ResData{2},['"player' num2str(player) '"']));
    % Race time
    RT = ResData{5}(PlInd);RT = str2num(RT{1}(2:end-1))

    Race.RaceTime = RT;
    Race.data = NaN;
    Race.EVENT.POS = NaN;
    Race.EVENT.TYP = NaN;
    Race.EVENT.DUR = NaN;
    Race.protocol = 'racetime';
    Race.AlignmentOffset = NaN;
    Race.gamecommands.type = NaN;
    Race.gamecommands.pos = NaN;
    Race.Date = thisGDate{1};
    Race.StartTime = NaN;
    Race.resfile = NaN;
    Race.GDFFile = NaN;
    Race.PlayerIndex = player;
    Race.GameIndexAll = NaN;
    Race.GameIndexFound = NaN;
    Race.Level = thisGLevel{1};
    Race.GameStartInd = NaN;
    Race.GameEndInd = NaN;

    % Save the Race output
    SessionRaceCounter = SessionRaceCounter + 1;

    % Creating folder if does not exist
    % 20161220 - ltonin
    [~, savepath] = cnbiutil_mkdir(TargetPath);

    save([savepath '/' SubID '.' thisGDate{1} '.' thisGTime{1} '.online.mi.mi_bhbf.racetime' num2str(SessionRaceCounter)  '.mat'  ],'Race');
    clear Race
    fclose all;         
end
