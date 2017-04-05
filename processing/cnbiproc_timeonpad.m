function [TOPPad TOPTask SesInd topall] = cnbiproc_timeonpad(padevents, PadTypeId, PadTypeLb, Dk)

if(nargin < 4)
    Dk = nan(1,padevents.POS(end)+padevents.DUR(end)+100); % Just a long useless vector with NaNs
end

SPadT = {'Speed','Jump','Slide','Rest'};

% The start and end pads are always a multiple of 18 (or we are in serious trouble)
if(mod(length(padevents.TYP),18) ~= 0)
    disp(['These events are not a full list of races!']);
    return;
end

% Find starting and ending pads
startInd = [];
for r=1:length(padevents.TYP)/18
    startInd = [startInd ; (r-1)*18+1];
end
endInd = [];
for r=1:length(padevents.TYP)/18
    endInd = [endInd ; r*18];
end
startendInd = union(startInd, endInd);


for p=1:3
    TOPPad{p} = padevents.DUR(ismember(padevents.TYP,PadTypeId(strcmp(PadTypeLb,SPadT{p}))))*0.0625;
    SesInd{p} = Dk(padevents.POS(ismember(padevents.TYP,PadTypeId(strcmp(PadTypeLb,SPadT{p})))));
end
TOPPad{4} = padevents.DUR(setdiff(find(padevents.TYP==783),startendInd))*0.0625;
SesInd{4} = Dk(padevents.POS(setdiff(find(padevents.TYP==783),startendInd)));
TOPPad{5} = padevents.DUR(startInd)*0.0625;
SesInd{5} = Dk(padevents.POS(startInd));
TOPPad{6} = padevents.DUR(endInd)*0.0625;
SesInd{6} = Dk(padevents.POS(endInd));

% Normally I should separate start/end from rest also here, but....screw it
for t=1:length(PadTypeId)
    TOPTask{p} = padevents.DUR(padevents.TYP==PadTypeId(t))*0.0625;
end

% Return them also in order
TOPPadAll = padevents.DUR*0.0625;
TOPPadAllLb = padevents.TYP;

TOPPadAllLb(find(ismember(TOPPadAllLb,PadTypeId(strcmp('Speed',PadTypeLb)))))=1;
TOPPadAllLb(find(ismember(TOPPadAllLb,PadTypeId(strcmp('Jump',PadTypeLb)))))=2;
TOPPadAllLb(find(ismember(TOPPadAllLb,PadTypeId(strcmp('Slide',PadTypeLb)))))=3;
TOPPadAllLb(find(ismember(TOPPadAllLb,PadTypeId(strcmp('Rest',PadTypeLb)))))=4;
TOPPadAllLb(startInd)=5;
TOPPadAllLb(endInd)=6;
topall.values = TOPPadAll;
topall.labels = TOPPadAllLb;