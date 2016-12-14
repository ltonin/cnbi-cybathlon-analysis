function [Date Time] = splitdatetime(datetime)

if(~isempty(strfind(datetime,'result')))
    datetime = datetime(6:end-3);
end

if strcmp(datetime,'""')
    Date=[];
    Time = [];
else
    Date = datetime(2:9);
    Time = datetime(10:end-1);
end


