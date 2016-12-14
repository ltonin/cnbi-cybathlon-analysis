function s = processlevels(s)

s(strfind(s,';'))=[];
if(length(s)==17)
    s=s(2:end); % Get rid of initial 0 that is given, but not always there
end
