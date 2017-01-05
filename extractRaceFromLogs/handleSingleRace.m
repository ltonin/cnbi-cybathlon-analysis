function Race = handleSingleRace(data, TYP, POS)

TriggersInd = find(ismember(TYP,[783 769 770 771 772 773 774 775]));

% Check that there are at least 5 commands during the race
% (excluding last pad)
CommandsInd = find(ismember(TYP,[10 20 30 40]));
RaceCommInd = intersect(CommandsInd,find(POS >= POS(TriggersInd(1))));
RaceCommInd = intersect(RaceCommInd,find(POS <= POS(TriggersInd(end))));

if(length(RaceCommInd) < 5)
    disp(['Skipping run ' Runs(run).name ' , too few commands found!']);
    Race = NaN;
    return;
end

% Remove potential weird 0 triggers
Ind0 = find(TYP==0);
TYP(Ind0) = [];
POS(Ind0) = [];
DUR(Ind0) = [];

% Find the race end knowing that at full speed (no false positive)
% the finish line is reached in 3 seconds, while at low speed (full false positives)
% the finish line is reached in 10 seconds
TimeOnLastPad = ...
decideRaceEnd(POS(intersect(find(POS >...
    POS(TriggersInd(end))),CommandsInd))/512-...
    POS(TriggersInd(end))/512);
disp(['TimeOnLastPad = ' num2str(TimeOnLastPad)]);