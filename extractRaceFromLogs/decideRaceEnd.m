function TimeOnLastPad = decideRaceEnd(CommandsAfterPad)

if(isempty(CommandsAfterPad))
    TimeOnLastPad = 3;
    return;
end

if(CommandsAfterPad(1) > 3)
    TimeOnLastPad = 3;
    return;
end

Start(1) = 0;
End(1) = CommandsAfterPad(1);
Slope(1) = 1;

i=2;
while(i <= length(CommandsAfterPad))
    
    if(CommandsAfterPad(i)-CommandsAfterPad(i-1) >= 4)
        % From the previous command and for 4 sec it goes slow
        Start = [Start ; End(end)];
        End = [End ; CommandsAfterPad(i-1)+4];
        Slope = [Slope ; 0.3];
        
        % From the previous point until current command, it goes fast
        Start = [Start ; End(end)];
        End = [End ; CommandsAfterPad(i)];
        Slope = [Slope ; 1];
    end
    i=i+1;
end
Start = [Start ; End(end)];
End = [End ; CommandsAfterPad(end)+4];
Slope = [Slope ; 0.3];
Start = [Start ; End(end)];
End = [End ; NaN];
Slope = [Slope ; 1];

% Build piecewise function
syms v t x
strexpr = 'v = 0';
for i=1:length(Start)-1
    strexpr = [ strexpr ' + (heaviside(t-' num2str(Start(i)) ') - heaviside(t-' num2str(End(i)) '))* ' num2str(Slope(i))];
end
strexpr = [ strexpr ' + (heaviside(t-' num2str(Start(end)) '))* ' num2str(Slope(end)) ';'];
eval(strexpr);
x = int(v,t)-3;

TimeOnLastPad = double(solve(x,t));