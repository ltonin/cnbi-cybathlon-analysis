function [TPFPPad, TPFPTask, SpeedPad, SpeedTask, FPPerSec] = cnbiproc_commacc(CommLb, TrialEvents, PadTypeId, PadTypeLb, PadTypeInd)

TPFPPad = zeros(4,2);
TPFPTask = zeros(length(PadTypeId),2);
SpeedPad  = cell(4,1);
SpeedTask  = cell(length(PadTypeId),1);
TotalINCDur = 0;
TotalFP = 0;
for tr=1:length(TrialEvents.TYP)
    TaskTyp = TrialEvents.TYP(tr);
    TrialInd = PadTypeInd(find(PadTypeId==TaskTyp));
    TrialCommLb = CommLb(TrialEvents.POS(tr):TrialEvents.POS(tr)+TrialEvents.DUR(tr));
    if(isempty(find(TrialCommLb~=0)))
        TrialCommLb(end) = 783 + hex2dec('6000');
    end
    
    if(TaskTyp == 783)
        TotalINCDur = TotalINCDur + TrialEvents.DUR(tr)*0.0625;
        TotalFP = TotalFP + sum(ismember(TrialCommLb,[25344 25347 25349]));
    end
    
    if(sum(TrialCommLb == (TaskTyp + hex2dec('6000')))>0)
        TPFPPad(TrialInd,1) = TPFPPad(TrialInd,1)+1;
        TPFPTask(find(PadTypeId==TaskTyp),1) = TPFPTask(find(PadTypeId==TaskTyp),1)+1;
        SpeedPad{TrialInd} = [SpeedPad{TrialInd} ; 0.0625*TrialEvents.DUR(tr)];
        SpeedTask{find(PadTypeId==TaskTyp)} = [SpeedTask{find(PadTypeId==TaskTyp)} ; 0.0625*TrialEvents.DUR(tr)];
    else
        TPFPPad(TrialInd,2) = TPFPPad(TrialInd,2)+1;
        TPFPTask(find(PadTypeId==TaskTyp),2) = TPFPTask(find(PadTypeId==TaskTyp),2)+1;
    end
end

TPFPPad = 100*TPFPPad./repmat(sum(TPFPPad,2),1,2);
TPFPTask = 100*TPFPTask./repmat(sum(TPFPTask,2),1,2);
FPPerSec = TotalFP/TotalINCDur;