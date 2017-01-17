function log = cnbiutil_read_logfile(logpath, targetname) 


    fid = fopen(logpath);
    
    tline = fgetl(fid);
    
    while ischar(tline)
        if(isempty(strfind(tline, targetname)) == false)
            classifier  = regexp(tline, '(?<=classifier=)\w*.mat', 'match');
            rejection   = regexp(tline, '(?<=rejection=)[-]?+\d+\.\d+', 'match');
            integration = regexp(tline, '(?<=integration=)[-]?+\d+\.\d+', 'match');
            thresholds  = regexp(tline, '(?<=thresholds=)\(\s(?<first>[-]?+\d+\.\d+)\s(?<second>[-]?+\d+\.\d+)\s\)', 'names');
            % TO DO: thresholds
            break
        end
        tline = fgetl(fid);
    end
    
    fclose(fid);
    
    log.classifier  = char(classifier);
    log.rejection   = str2double(rejection);
    log.integration = str2double(integration);
    log.thresholds = [str2double(thresholds.first) str2double(thresholds.second)];

end