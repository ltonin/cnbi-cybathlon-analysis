function DT = alignGameGDFSeries(GameTime, GameType, GDFTime, GDFType)

% Transform types to 1/2/3
% GDF
NewGDFType = [];
if(length(find(ismember(GDFType,[5 6 7])))>0)
    NewGDFType = GDFType-4; % 5 = Speed, 6=Jump, 7=Kick
elseif(length(find(ismember(GDFType,[25347 25349])))>0)
    % TODO
    NewGDFType(find(GDFType==25347)) = 2; % 25347 = 771 = Both Feet = Jump (or Kick...)
    NewGDFType(find(GDFType==25349)) = 1; % 25349 = 773 = Both Hands = Speed (or Kick...)
else
    disp('UNKNOWN PROTOCOL!');
    DT = NaN;
    return;
end

% Game
NewGameType = zeros(1,length(GameType));
NewGameType(find(not(cellfun('isempty', strfind(GameType,'Speed'))))) = 1;
NewGameType(find(not(cellfun('isempty', strfind(GameType,'Jump'))))) = 2;
NewGameType(find(not(cellfun('isempty', strfind(GameType,'Kick'))))) = 3;

% Slide GDF onto Game
pDiff = {};
MatchIndex = [];
TD = [];
i=0;
while(true)
    i=i+1;
    % Find time difference
    TD(i) = GameTime(1) - GDFTime(i);
    tmpGameTime = GameTime - TD(i);
    if(tmpGameTime(end) > max(GDFTime))
        break;
    end
    
    % Check all game commands, whether they have a GDF point really close
    % in time. I check all GDF points and not simply in order (as I had
    % done originally), because fucking game log has double entries (like
    % it happens for the triggers). By checking only for Matched GAME points, 
    % there is no problem if I have GDF events that led to no command (like
    % for instance when the artifact rejection was on!)
    pDiff{i} = [];
    for pGame=1:length(tmpGameTime)
        pDiff{i} = [pDiff{i} min(abs(GDFTime - tmpGameTime(pGame))) ];
    end
    %MatchIndex(i) = mean(pDiff{i});
    MatchIndex(i) = max(pDiff{i}); % Make it even more strict!
end

MatchPoint =  find(MatchIndex < 0.1);
if(length(MatchPoint)>1)
    disp(['Found more than one matching points! This is degenerate situations that means we had the same commands all the time. Skipping this.']);
    DT = NaN;
    return
end
if(~isempty(MatchPoint))
    DT = TD(MatchPoint);
else
    DT = NaN;
end


