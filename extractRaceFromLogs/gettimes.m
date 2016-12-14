function s = gettimes(s)


if strcmp(s,'""')
    s=[];
else
    s= str2num(s(2:end-1));    
end
