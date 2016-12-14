function s = pads2str(P)

s='';
for i=1:length(P)
    switch P{i}
        case 'Walk'
            s(i) = '0';
        case 'Speed'
            s(i) = '1';
        case 'Jump'
            s(i) = '2';
        case 'Kick'
            s(i) = '3';
        otherwise
            s(i) = 'X';
    end
end