function GrouppedTimes = groupTriggerTimes(Times)

% Sort times ascending
Times = sort(Times, 'ascend');
GrouppedTimes = Times(1);
for i=2:length(Times)
    if(Times(i)-Times(i-1) > 0.5)
        GrouppedTimes = [GrouppedTimes; Times(i)];
    end
end

