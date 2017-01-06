function info = cnbiutil_getfile_info(filename)
% info = cnbiutil_getfile_info(filename)
%
% Given a standard filename format, the function return a structure with
% fields retrieved from the filename. Standard filename format is defined
% as follows:   SUBJECT.DATE.TIME.MODALITY.TASK.TASKSET.EXTENSION
%
% The returned structure has the following fields:
%
%   info
%       .subject
%       .date
%       .time
%       .modality
%       .task
%       .taskset
%       .extra
%       .filepath
%       .extension
%
% Example:
%
% filename = /mnt/data/Research/cvsa/online/20120704_b4//b4.20120704.124619.offline.va.va_brbl.bdf;
% info = util_getfile_info(filename);
%
% info = 
% 
%       subject: 'b4'
%          date: '20120704'
%          time: '124619'
%      modality: 'offline'
%          task: 'va'
%       taskset: 'va_brbl'
%       extra:   '';
%      filepath: '/mnt/data/Research/cvsa/online//20120704_b4/'
%     extension: '.bdf'
%

    [path, name, ext] = fileparts(filename);
    
    fields = regexp(name, '\.', 'split');
    
    info.subject    = fields{1};
    info.date       = fields{2};
    info.time       = fields{3};
    info.modality   = fields{4};
    info.task       = fields{5};
    info.taskset    = fields{6};
    info.extra      = '';
    if length(fields) == 7
        info.extra      = fields(7);
    elseif length(fields) > 7
        info.extra = fields(7:length(fields));
    end
    info.filepath   = path;
    info.extension  = ext;
end