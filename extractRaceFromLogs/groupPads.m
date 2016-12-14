function [FPads FStrPads] =  groupPads(Pads, RCT)

StrPads = pads2str(Pads);

FStrPads{1} = StrPads(1);
FPads{1} = Pads{1};
f=1;
for i=2:length(Pads)
    if( (RCT(i) - RCT(i-1)) < 0.2)
        FPads{f} = [ FPads{f} '/' Pads{i}];
        FStrPads{f} = [ FStrPads{f} '/' StrPads(i)];
    else
        f=f+1;
        FPads{f}  = Pads{i};
        FStrPads{f}  = StrPads(i);
    end
end