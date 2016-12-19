function race2csv(filename,tableMat)

load(tableMat);
cellArray = table2cell(raceLog);

datei = fopen(filename,'w');
for z=1:size(cellArray,1)
    for s=1:size(cellArray,2)

        var = eval(['cellArray{z,s}']);

        if size(var,1) == 0
            var = '';
        end

        if isnumeric(var) == 1
            var = num2str(var);
        end
        if (s == size(cellArray,2))
            fprintf(datei,'"%s"\n',char(var));
        elseif(s==3)
            fprintf(datei,'"P%s",',char(var));
        else
            fprintf(datei,'"%s",',char(var));
        end
    end
end
fclose(datei);
