function s = trimquotes(s)

if strcmp(s,'""')
    s='';
else
    s= s(2:end-1);    
end
