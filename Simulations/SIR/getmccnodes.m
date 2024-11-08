% create lisco files
tree_files = dir('master/*.tree');

system('rm -r mcc');
system('mkdir mcc');
for i = 1 : length(tree_files)
    if contains(tree_files(i).name, 'single') || contains(tree_files(i).name, 'original') 
        continue
    end
    disp(tree_files(i).name)
    % read tree files
    g = fopen(['master/' tree_files(i).name],'r');
    t = textscan(g,'%s'); fclose(g);

    % get the number of the tree files that is S%d
    tmps = strsplit(tree_files(i).name,'_');
    tmps = strrep(tmps{2}, 'S','');

    
    if (str2double(tmps)<801)
    
        % coalescing
        tree_tmp3 = regexprep(t{1}(end-1),'&type="I",location="(\d*)",reaction="Sampling",time=(\d*).(\d*)','');
        tree_tmp2 = regexprep(t{1}(end-1),'&type="I",location="(\d*)",reaction="Infection",time=(\d*).(\d*)','');


        %migrating
        tree_tmp1 = regexprep(tree_tmp2,'&type="I",location="(\d*)",reaction="Migration",time=(\d*).(\d*)','');

        
        % make the MASTER tree compatible with BEAST2
        % sampling
        tree_tmp1 = regexprep(tree_tmp1,'E[-](\d)]',']');
        tip_locs = regexp(tree_tmp1,'[&type="I",location="(\d*)",reaction="Sampling",time=(\d*)\.(\d*)\]','match');
        for j = 1 : length(tip_locs{1})
            loc = regexprep(tip_locs{1}{j},'[&type="I",location="(\d*)",reaction="Sampling",time=(\d*)\.(\d*)\]','$1');
            time = regexprep(tip_locs{1}{j},'[&type="I",location="(\d*)",reaction="Sampling",time=(\d*)\.(\d*)\]','$2.$3');
            tree_tmp1 = strrep(tree_tmp1,tip_locs{1}{j},['loc_' loc '_time_' time 'kickout']);
            tree_tmp1 = strrep(tree_tmp1,'kickout','');
        end

        tree_tmp = regexprep(tree_tmp1,'(\d*)loc_','inv$1loc_');
        
        tree = strrep(tree_tmp{1},'[]','');
        if ~isempty(strfind(tree,']'))
            b = strfind(tree,']');
            c = tree((b-50):(b+50));
            disp(tree_files(i).name)
        end
    else
        continue;
    end


    % get the leafnames
    ptree = phytreeread(tree);
    leafnames = get(ptree,'leafnames');
    % for all leafnames, remove the _time_.... part
    for j = 1 : length(leafnames)
        leafnames{j} = regexprep(leafnames{j},'_time_\d*.\d*','');
    end

    
    print_tree = tree;

    % open a new file that is named target
    g = fopen(['mcc/' strrep(tree_files(i).name, 'master','target')], 'w');
    fprintf(g, '#NEXUS\n\n');

    fprintf(g, 'Begin taxa;\n');
    fprintf(g, '\tDimensions ntax=%d;\n', length(leafnames));
    fprintf(g, '\t\tTaxlabels\n');
    for j = 1 : length(leafnames)
        fprintf(g, '\t\t\t%s\n', leafnames{j});
    end
    fprintf(g, ';\n');
    fprintf(g, 'End;\n');
    fprintf(g, 'Begin trees;\n');
    fprintf(g, '\tTranslate\n');
    for j = 1 : length(leafnames)-1
        fprintf(g, '\t\t%d %s,\n', j, leafnames{j});
    end
    fprintf(g, '\t\t%d %s\n', length(leafnames), leafnames{length(leafnames)});
    fprintf(g, ';\n');
    fprintf(g, 'tree STATE_0 =  %s\n', t{1}{end-1});
    fprintf(g, 'End;\n');
end