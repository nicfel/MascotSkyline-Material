%% 
clear

g = fopen('nodesupport.tsv', 'w');
fprintf(g, 'method\trun\tsupport\n');


files=dir('./mcc/*mascot.tree');


for i = 1: length(files)
    f = fopen(['./mcc/' files(i).name]);
    tmp2 = strsplit(files(i).name, '.');
    values = zeros(0,0);

    while ~feof(f)
        line = fgets(f);
        if contains(line, 'tree TREE1')
            tmp = strsplit(line, '[&');
            for j=2:length(tmp)
                if contains(tmp{j}, 'Sampling')
                    continue;
                end
                true_val = regexp(tmp{j}, 'location=(\d*)', 'match');
                true_val = strsplit(true_val{1}, '=');
                support = regexp(tmp{j}, 'state(\d)=(\d*).(\d*)', 'match');
                support = strsplit(support{str2double(true_val{2})+1}, '=');
                if str2double(support{2})>1
                    values(end+1) = 0.0;
                else
                    values(end+1) = str2double(support{2});
                end
                fprintf(g, 'mascot\t%s\t%f\n', strrep(tmp2{1}, 'S',''), values(end));


            end 
        end
    end
    fclose(f);

    
    values = zeros(0,0);

    f = fopen(['./mcc/' strrep(files(i).name, 'mascot','dta') ]);
    while ~feof(f)
        line = fgets(f);
        if contains(line, 'tree TREE1')
            tmp = strsplit(line, '[&');
            for j=2:length(tmp)
                if contains(tmp{j}, 'Sampling')
                    continue;
                end
                
                true_val = regexp(tmp{j}, 'location=(\d*)', 'match');
                true_val = strsplit(true_val{1}, '=');
                
                tmp3 = strsplit(tmp{j}, '{');
                tmp4 = strsplit(tmp3{2}, '}');
                tmp5 = strsplit(tmp3{3}, '}');
                
                states = strsplit(tmp4{1}, ',');
                probs = strsplit(tmp5{1}, ',');
                
                if ismember(['state' true_val{2}], states)
                    values(end+1) = str2double(probs{ismember(states,['state' true_val{2}])});
                else
                    values(end+1) = 0;
                end
                
                fprintf(g, 'dta\t%s\t%f\n', strrep(tmp2{1}, 'S',''), values(end));



            end 
        end
    end
    fclose(f);
%     fprintf(g, 'dta\t%s\t%f\t%f\n', strrep(tmp2{1}, 'S',''), mean(values), values(end));


end

fclose(g);