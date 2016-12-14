function [] = syncGDFlog(SubID, GDFSesPath, LogFilePath)

% Add biosig
addpath(genpath('~/Git/cnbi-smrtrain/toolboxes/biosig/'));

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

% For all players
for p=1:4
    IndTrigP{p} = intersect(find(strcmp(PlayerTrack,['P' num2str(p)])),find(strcmp(Trigger,'trigger')));
    tmpCommSpeed = intersect(find(strcmp(PlayerTrack,['P' num2str(p)])),find(strcmp(Commands,'Speed')));
    tmpCommJump = intersect(find(strcmp(PlayerTrack,['P' num2str(p)])),find(strcmp(Commands,'Jump')));
    tmpCommKick = intersect(find(strcmp(PlayerTrack,['P' num2str(p)])),find(strcmp(Commands,'Kick')));
    IndActiveCommandP{p} = union(tmpCommSpeed,tmpCommJump);  
    IndActiveCommandP{p} = union(IndActiveCommandP{p},tmpCommKick);
    IndP{p} = find(strcmp(PlayerTrack,['P' num2str(p)]));
end


% Find possible games

% This is if we consider everything between app/startgame and any next app (game closed, reset, etc)
GameStartInd = find(strcmp(Commands,'startgame'));
GameEndInd = setdiff(find(strcmp(Trigger,'app')),GameStartInd);
StartEndInd = union(GameStartInd,GameEndInd);
StartEndIndVal = zeros(1,length(StartEndInd));
StartEndIndVal(find(ismember(StartEndInd,GameStartInd)))=1;
StartEndIndVal(find(ismember(StartEndInd,GameEndInd)))=2;
GameStartInd  = StartEndInd(strfind(StartEndIndVal,[1 2]));
GameEndInd  = StartEndInd(strfind(StartEndIndVal,[1 2])+1);

% Games found
GamesFound = 0;
% Search for each game
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
    
    % First, find the file where this game should be in
    % Check if there is such session folder:
    thisGDate = RaceDate(GameStartInd(g));
    thisGTime = RaceTime(GameStartInd(g));
    thisSesFold = [GDFSesPath '/' SubID '_' thisGDate{1}];
    if(exist(thisSesFold,'dir')==7)
        % Candidate runs
        thisCandRuns = dir([thisSesFold '/' SubID '.' thisGDate{1} '.*.online.mi.mi_bhbf.gdf']);
        if(isempty(thisCandRuns))
            disp('There are no game runs in this session. Skipping game');
            continue;
        end
        
        % Check if there are one or more suitable run times, BEFORE the game start time
        % For this, first find the run times!
        thisRunTimes = [];
        thisGTimeNum = datenum(datetime(thisGTime{1},'InputFormat','HHmmss'));
        for r=1:length(thisCandRuns)
            dots = strfind(thisCandRuns(r).name,'.');
            thisRunTimes(r) = datenum(datetime(thisCandRuns(r).name(dots(2)+1:dots(3)-1),'InputFormat','HHmmss'));
        end
        % Check if there are negative values (otherwise all runs were recorded AFTER this game)
        %if(isempty(find(thisRunTimes-thisGTimeNum < 0)))
        %    disp('There are no game runs started before this game in this session. Skipping game');
        %    continue;
        %end
        % Find smallest negative difference
        [~, minInd] = min(1./(thisRunTimes-thisGTimeNum));
        % OK! So, this should be the run I should look into!
        thisChosenRun = [GDFSesPath '/' SubID '_' thisGDate{1} '/' thisCandRuns(minInd).name]; 
        % Do not use this because the laptop times were not aligned. Check
        % every run in the same date
    else
        disp(['There is no session ' thisGDate{1} ' in the database. Skipping game']);
        continue;
    end
    
    % Get indices (tracks) of existent players
    GPlayerInd = [];
    for p=1:length(GPlayers)
        GPlayerInd(p) = str2num(GPlayers{p}(2));
    end    
    % Search for each player
    for p=GPlayerInd
        disp(['Checking game ' num2str(g) ' and Player ' num2str(p)]);
        
        % Find all trigger indices for this player and game
        thisPGtrigInd = intersect(IndTrigP{p}(find(IndTrigP{p}>=GameStartInd(g))),...
            IndTrigP{p}(find(IndTrigP{p}<=GameEndInd(g))));
        thisPGtrigTime = RaceCommandTime(thisPGtrigInd);
        
        % Group the triggers
        try
            [thisPGtrigInd, thisPGtrigTime] = groupTriggers(thisPGtrigInd,thisPGtrigTime);
        catch
            disp('Skipping this player and game, problem grouping triggers!!');
            continue;
        end
        
%         % Group timings of triggers, since they are fucked up in the log file and there are extra ones that should not be there
%         % Since it seems the extra ones appear immediately, they could be
%         % considred as one. Note that the pad of field Pads cannot be
%         % trusted, that's why I use the loadlevel entry for each game race
%         length(thisPGtrig) - length(find(diff(thisPGtrigTime) <0.2 )) % Below 200 msec it cannot possibly be a new trigger!
        % Even after grouping there are misalignments. I think I can only
        % trust the loadlevel, and then simply assume that timings are
        % correct (as many of them as can be found)
        
        % Now find commands of this player to be used for the alignment
        thisPGActiveCommandsInd = intersect(IndActiveCommandP{p}(find(IndActiveCommandP{p}>=GameStartInd(g))),...
            IndActiveCommandP{p}(find(IndActiveCommandP{p}<=GameEndInd(g))));
        thisPGActiveCommands = Commands(thisPGActiveCommandsInd);
        thisPGActiveCommandsTime = RaceCommandTime(thisPGActiveCommandsInd);
        if(length(thisPGActiveCommands) < 5)
            disp('Skipping this game/player due to too few commands found!');
            continue;
        end
        
        for gd=1:length(thisCandRuns)
            disp(['Checking game ' num2str(g) ' and Player ' num2str(p) ' and run: ' thisCandRuns(gd).name]);
            %% Now load the header of the chosen run
            try
                thisHeader = sopen([GDFSesPath '/' SubID '_' thisGDate{1} '/' thisCandRuns(gd).name]); 
            catch
                disp(['This run is corrupted. Skipping']);
                continue;
            end
            % Process the loaded GDF data
            thisGDFTimings = thisHeader.EVENT.POS(find(ismember(thisHeader.EVENT.TYP,[5 6 7 25347 25349])))/512;
            thisGDFTypes = thisHeader.EVENT.TYP(find(ismember(thisHeader.EVENT.TYP,[5 6 7 25347 25349])));

            % Now attempt alignment
            DT = alignGameGDFSeries(thisPGActiveCommandsTime,thisPGActiveCommands,...
            thisGDFTimings,thisGDFTypes);
            % Translate the trigger times according to alignment
            thisPGtrigTimeAligned = thisPGtrigTime - DT;
            thisPGActiveCommandsTimeAligned = thisPGActiveCommandsTime - DT;
            thisPGStartGameAligned = RaceCommandTime(GameStartInd(g))-DT;
            thisPGEndGameAligned = RaceCommandTime(GameEndInd(g))-DT;
            % Convert triggers, commands and start into GDF POS indices
            thisTrigPOS = ceil(thisPGtrigTimeAligned*512);
            thisCommPOS = ceil(thisPGActiveCommandsTimeAligned*512);
            thisStartPOS = ceil(thisPGStartGameAligned*512);
            thisEndPOS = ceil(thisPGEndGameAligned*512);

            if(~isnan(DT))  
            
                GamesFound = GamesFound +1;
                ShowDT(GamesFound) = DT;
                % Check if trigger missing (they should be 17 = 16 + ending pad (no trigger for starting pad, this can be retrived form startgame ). Often we get 16 or less)
                if(length(thisTrigPOS) == 17)
                    disp(['All triggers found!']);
                elseif(length(thisTrigPOS) > 17)
                    disp(['Found more tirggers than expected! That ''' 's a new one!!! Skipping...']);
                    GamesFound = GamesFound -1;
                    continue;
                elseif(length(thisTrigPOS) < 17) % We are missing triggers for sure sespite syncing is possible

                    disp(['Missing ' num2str(17 - length(thisTrigPOS)) ' triggers!!']);

                    % Check if there is a levelElementEnd(Clone) as a sign
                    % that this player really finished. Otherwise, it seems
                    % that simply this is an incomplete game, we would
                    % rather not take it into account
                    HasEnd = length(intersect(find(intersect(IndP{p},find(strcmp(Pad,'End'))) >...
                        GameStartInd(g)),find(intersect(IndP{p},find(strcmp(Pad,'End'))) < GameEndInd(g))));
                    if(HasEnd == 0)
                        disp(['Missing finish line event, this must be an incomplete race, skipping...']);
                        GamesFound = GamesFound -1;
                        continue;
                    else
                        disp(['There is a finish line event for this player, this must be a complete race! Trying to reinstate missing triggers']);
                    end
                    
                    % Find Kick command indices in this game
                    thisPGKickCommInd = thisPGActiveCommandsInd(find(strcmp(thisPGActiveCommands,'Kick')));
                    thisPGPadsInd = intersect(IndActiveCommandP{p}(find(IndActiveCommandP{p}>=GameStartInd(g))),...
                        IndActiveCommandP{p}(find(IndActiveCommandP{p}<=GameEndInd(g))));
                    % Find Kick pads in this game
                    IndKickP = intersect(IndP{p},find(strcmp(Pad,'Kick')));
                    thisPGKickPadInd = intersect(IndKickP(find(IndKickP >= GameStartInd(g))),IndKickP(find(IndKickP <= GameEndInd(g))));
                    MissingTrigCandidates = intersect(thisPGKickCommInd,thisPGKickPadInd);
                    % Remove double entries
                    NewMissingTrigCandidates = [];
                    NewMissingTrigCandidates = MissingTrigCandidates(1);
                    for m=2:length(MissingTrigCandidates)
                        if( (RaceCommandTime(MissingTrigCandidates(m)) - RaceCommandTime(MissingTrigCandidates(m-1))) > 0.2)

                            NewMissingTrigCandidates = [ NewMissingTrigCandidates MissingTrigCandidates(m)];
                        end
                    end
                    MissingTrigCandidates = NewMissingTrigCandidates;clear NewMissingTrigCandidates;
                    
                    % There must be a "Run" very close after "Kick", say
                    % 2.0 sec and there must be no trigger soon after, say, 0.1
                    % sec after the Run
                    thisGInd = [GameStartInd(g):GameEndInd(g)]';
                    thisPGRunInd = intersect(intersect(IndP{p},thisGInd),find(strcmp(Commands,'Run')));
                    
                    FinalCandidateInd = [];
                    for mm=1:length(MissingTrigCandidates)
                        tmpDiffTime = RaceCommandTime(thisPGRunInd)'-RaceCommandTime(MissingTrigCandidates(mm));
                        tmpDiffTime(tmpDiffTime<=0) = Inf;
                        [minValRun minIndRun] = min(tmpDiffTime);
                        if(minValRun <= 2.0)
                            % There is Run very close afterwards, look now for trigger
                            % very close
                            if(isempty(find(min(abs(RaceCommandTime(thisPGtrigInd)-RaceCommandTime(thisPGRunInd(minIndRun)))) < 0.05)))
                                % There is no trigger very close, so, this
                                % must be a missing trigger position! Save
                                % it
                                FinalCandidateInd = [FinalCandidateInd ; [MissingTrigCandidates(mm) thisPGRunInd(minIndRun)]];
                            end
                        end                       
                    end
                    if(size(FinalCandidateInd,1) == (17 - length(thisTrigPOS)))
                        disp('There are equal number of missing triggers and candidate positions! Reinsating triggers');
                        thisTrigPOS = [thisTrigPOS ceil((RaceCommandTime(FinalCandidateInd(:,2))-DT)*512)'];
                        thisTrigPOS = sort(thisTrigPOS,'ascend');
                    else
                        if(strcmp([thisGDate{1} thisGTime{1}],'20160831151542'))
                            FinalCandidateInd = [57257 57262]; % This is a race of Eric where there is a second bad candidate AFTER the end of the race...
                        else
                            disp('There are NOT equal numbers of missing triggers and candidate positions! Skipping...');
                            GamesFound = GamesFound -1;
                            continue;
                        end
                    end    
                end

                % Add in the beginning the StartGame trigger
                thisTrigPOS = [ceil((RaceCommandTime(GameStartInd(g))-DT)*512) thisTrigPOS];

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
                PlInd = find(strcmp(ResData{2},['"player' num2str(p) '"']));
                % Race time
                RT = ResData{5}(PlInd);RT = str2num(RT{1}(2:end-1));

                % Race end POS trigger
                thisTrigPOS = [thisTrigPOS thisTrigPOS(1)+ceil(RT*512)];
                %sclose(thisHeader);

                % Make final output
                FinalRaceTimes(GamesFound) = RT;

                % load the GDF file since we want to also save the data
                [data header] = sload([GDFSesPath '/' SubID '_' thisGDate{1} '/' thisCandRuns(gd).name]);
                %Crop data between starting and ending index
                Race.eeg = data(thisTrigPOS(1):thisTrigPOS(end),:);
                IndUseful = intersect(find(header.EVENT.POS >= thisTrigPOS(1)) , find(header.EVENT.POS <= thisTrigPOS(end)));
                Race.trigger.pos = header.EVENT.POS(IndUseful) - thisTrigPOS(1);
                Race.trigger.type = header.EVENT.TYP(IndUseful);

                % Remove triggers from UDP messages in case they exist
                % They are useful for debugging but the gameevetns are
                % more accurate wrt timing
                RemInd = ismember(Race.trigger.type,[783 768 771 773]);
                if(~isempty(RemInd))
                    Race.trigger.pos(RemInd) = [];
                    Race.trigger.type(RemInd) = [];
                end

                % Convert pad triggers
                thisTrigPOS = thisTrigPOS - thisTrigPOS(1);
                thisGLevelInt = nan(1,16);
                thisGLevelInt(find(thisGLevel{1}=='0')) = 783;
                thisGLevelInt(find(thisGLevel{1}=='1')) = 773;
                thisGLevelInt(find(thisGLevel{1}=='2')) = 771;
                thisGLevelInt(find(thisGLevel{1}=='3')) = 768;
                thisTrigTYP = [100 thisGLevelInt 200 300];
                tmpAllPOS = [thisTrigPOS Race.trigger.pos'];
                tmpAllTYP = [thisTrigTYP Race.trigger.type'];
                [sorted sortInd] = sort(tmpAllPOS,'ascend');
                Race.trigger.pos = tmpAllPOS(sortInd);
                Race.trigger.type = tmpAllTYP(sortInd);
                Race.AlignmentOffset = DT;
                Race.RaceTime = RT;
                % Add also the game commands, can be useful to know
                % which of the bh and bf led to "Slide"
                Race.gamecommands.type = thisPGActiveCommands;
                Race.gamecommands.pos = ceil((thisPGActiveCommandsTimeAligned - ...
                    (RaceCommandTime(GameStartInd(g))-DT))*512);
                Race.Date = thisGDate{1};
                Race.StartTime = thisGTime{1};
                Race.resfile = thisResFile;
                Race.GDFFile = thisCandRuns(gd).name;
                Race.PlayerIndex = p;
                Race.GameIndexAll = g;
                Race.GameIndexFound = GamesFound;
                Race.Level = thisGLevel{1};
                Race.GameStartInd = GameStartInd;
                Race.GameEndInd = GameEndInd;
                % Save the Race output
                save([GDFSesPath '/RaceMat/Race_' num2str(Race.GameIndexFound) '_' Race.GDFFile(1:end-4) '.mat'  ],'Race');
                clear Race data header thisGLevelInt tmpAllPOS tmpAllTYP thisTrigPOS thisTrigTYP IndUseful
            end        
        end
    end
end

disp(['Found ' num2str(GamesFound) '/' num2str(length(GameStartInd))]);