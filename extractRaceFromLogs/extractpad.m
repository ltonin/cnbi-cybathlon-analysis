function pad = extractpad(pad)

if strcmp(pad,'""')
    pad='';
else
    ParenthesisL = strfind(pad,'(');
    pad = pad(14:ParenthesisL-1);
end
