%% build mcc target trees
clear
% create lisco files
tree_files = dir('master/*.tree');

system('rm -r mcc');
system('mkdir mcc');
for i = 1 : length(tree_files)
    if contains(tree_files(i).name, 'single')
        continue
    end
    
    f = fopen(sprintf('../out/%s', strrep(tree_files(i).name, 'master.tree','1mascot.H3N2.trees')));
    g = fopen(['mcc/' tree_files(i).name], 'w');
    while ~feof(f)
        line = fgets(f);
        if contains(line, 'tree STATE_0')
            break;
        end
        fprintf(g, line);
    end
    fclose(f);
    f = fopen(['master/' tree_files(i).name]);
    while ~feof(f)
        line = fgets(f);
        if contains(line, 'tree TREE')
            fprintf(g, line);
            fprintf(g, 'End;');            
        end        
    end
    fclose('all');
    
    tmp = strsplit(tree_files(i).name, '_');

    system(sprintf('/Applications/BEAST\\ 2.7.1/bin/treeannotator -target mcc/SIR_%s_master.tree ../out/SIR_%s_1mascot.H3N2.trees mcc/%s.mascot.tree',...
    tmp{2}, tmp{2}, tmp{2}));

    system(sprintf('/Applications/BEAST\\ 2.7.1/bin/treeannotator -target mcc/SIR_%s_master.tree ../out/DTASIR_%s_1dta.trees mcc/%s.dta.tree',...
    tmp{2}, tmp{2}, tmp{2}));
    
end