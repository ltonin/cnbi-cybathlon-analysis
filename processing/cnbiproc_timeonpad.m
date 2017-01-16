function [TOPPad TOPTask] = cnbiproc_timeonpad(padevents, PadTypeId, PadTypeLb)

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
end
TOPPad{4} = padevents.DUR(setdiff(find(padevents.TYP==783),startendInd))*0.0625;
TOPPad{5} = padevents.DUR(startInd)*0.0625;
TOPPad{6} = padevents.DUR(endInd)*0.0625;

% Normally I should separate start/end from rest also here, but....screw it
for t=1:length(PadTypeId)
    TOPTask{p} = padevents.DUR(padevents.TYP==PadTypeId(t))*0.0625;
end