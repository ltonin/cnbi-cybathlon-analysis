SubList = {'AN14VE', 'MA25VE'};
NumSubjects = length(SubList);

datapath  = [pwd '/analysis/'];
figuredir = './figures/';

for sId = 1:NumSubjects
    
    csubject = SubList{sId};

    cracetimefilepath = [datapath '/' csubject '.time.race.mat'];
    cpadaccfilepath = [datapath '/' csubject '.padaccuracy.race.mat'];
    cdiscrfilepath = [datapath '/' csubject '.discriminancy.maps.mat'];
    
    cdatatime = load(cracetimefilepath);
    cdatapadacc = load(cpadaccfilepath);
    cdatadiscr = load(cdiscrfilepath);
    
    crt = cdatatime.time.run.values;
    crtDk = cdatatime.time.run.label.Dk;
    crtDl = mat2cell(cdatatime.time.run.label.Dl,ones(1,size(cdatatime.time.run.label.Dl,1)),8);
    
    cpa = cdatapadacc.pad.accuracy.run.values;
    cpaDk = cdatapadacc.pad.label.run.Dk;
    cpaDl = mat2cell(cdatapadacc.pad.label.session.Dl,ones(1,size(cdatapadacc.pad.label.session.Dl,1)),8);
    
    cfs = cdatadiscr.discriminancy.run.fisherscore;
    cfsDk = cdatadiscr.discriminancy.run.label.Dk;
    cfsDl = mat2cell(cdatadiscr.discriminancy.run.label.Dl,ones(1,size(cdatadiscr.discriminancy.run.label.Dl,1)),8);
    
    % Create the common attributes
    AllSessions = union(union(crtDl,cpaDl),cfsDl);
    
    ss = 0;
    for s=1:length(AllSessions)
        % Check if this session exists forall metrics
        ExistsIn = sum(strcmp(AllSessions{s},crtDl))+sum(strcmp(AllSessions{s},cpaDl))...
            +sum(strcmp(AllSessions{s},cfsDl));
        if(ExistsIn == 3)
            ss=ss+1;
            crtF{ss} = crt(find(crtDk==find(strcmp(AllSessions{s},crtDl))));
            cpaF{ss} = nanmean(cpa(1:3,find(cpaDk==find(strcmp(AllSessions{s},cpaDl)))));
            cfsF{ss} = cfs(:,find(cfsDk==find(strcmp(AllSessions{s},cfsDl))),4);
            disp('a');
        end
    end
    disp('a');
end
