function [NewtrigInd, NewtrigTime, NewtrigPad] = groupTriggers(trigInd, trigTime)

NewtrigInd(1) = trigInd(1);
NewtrigTime(1) = trigTime(1);

for i=2:length(trigInd)
    if( (trigTime(i) - trigTime(i-1)) > 0.2)
        NewtrigInd = [ NewtrigInd trigInd(i)];
        NewtrigTime = [ NewtrigTime trigTime(i)];
    end
end