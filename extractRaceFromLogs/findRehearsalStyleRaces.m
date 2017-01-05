function findRehearsalStyleRaces(SubID, GDFSesPath)

% Add biosig
addpath(genpath('~/Git/cnbi-smrtrain/toolboxes/biosig/'));

% Find existing session folders
SesFolders = dir([GDFSesPath '/' SubID '_*']);

for ses=1:length(SesFolders)
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
        disp('a');        
        
    end


end