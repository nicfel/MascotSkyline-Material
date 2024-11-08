%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% creates structured coalescent xmls from the master trees. Always creates
% 3 xmls per tree with different initial values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
% create lisco files
tree_files = dir('master/*.tree');

system('rm -r xmls');
system('mkdir xmls');

states = 2;

% read in the file rates.txt that is a tsv table, with the first columb being the name of the run, 
% the second the avg migration rate and the third the sample number, all other columns are not needed.
rates = readtable('rates.txt','Delimiter','\t', 'HeaderLines',1);


for i = 1 : length(tree_files)
    if contains(tree_files(i).name, 'single') || contains(tree_files(i).name, 'original') || contains(tree_files(i).name, 'background')
        continue
    end

    disp(tree_files(i).name)
    % read tree files
    g = fopen(['master/' tree_files(i).name],'r');
    t = textscan(g,'%s'); fclose(g);

    % get the number of the tree files that is S%d
    tmps = strsplit(tree_files(i).name,'_');
    tmps = strrep(tmps{2}, 'S','');

    % from rates, get the column with the last instance of the first variable = tmps
    rate = rates.Var2(find(rates.Var1==str2double(tmps)));   

    if (str2double(tmps)<801)
        % coalescing
        tree_tmp3 = regexprep(t{1}(end-1),'&type="I",location="(\d*)",reaction="Sampling",time=(\d*).(\d*)','');
        tree_tmp2 = regexprep(t{1}(end-1),'&type="I",location="(\d*)",reaction="Infection",time=(\d*).(\d*)','');    
        % migrating
        tree_tmp1 = regexprep(tree_tmp2,'&type="I",location="(\d*)",reaction="Migration",time=(\d*).(\d*)','');        
        % make the MASTER tree compatible with BEAST2
        % sampling
        tree_tmp1 = regexprep(tree_tmp1,'E[-](\d)]',']');
        tip_locs = regexp(tree_tmp1,'[&type="I",location="(\d*)",reaction="Sampling",time=(\d*)\.(\d*)\]','match');
        for j = 1 : length(tip_locs{1})
            loc = regexprep(tip_locs{1}{j},'[&type="I",location="(\d*)",reaction="Sampling",time=(\d*)\.(\d*)\]','$1');
            tree_tmp1 = strrep(tree_tmp1,tip_locs{1}{j},['loc_' loc 'kickout']);
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
        tree = t{1}{1};
        if str2double(tmps)<901
            rate = 5;
        else
            rate = 25;
        end
        % remove everything of the form _time_"number" or "number.number"
        tree = regexprep(tree,'_time_(\d*).(\d*)','');
    end


    % get the leafnames
    ptree = phytreeread(tree);
    leafnames = get(ptree,'leafnames');
    
    print_tree = tree;
    
    dist = pdist(ptree,'nodes','all','squareform',true);
    
    % get the tree height
    root_height = max(dist(end,1:length(leafnames)));

    
    
    % get the covariates
    
    % make tripletts of all runs with different random initial values
    for tr = 1 : 1    
        g = fopen('Template.xml');
        
        flog = strrep(tree_files(i).name,'master.tree',sprintf('%dmascot',tr));
        fname = sprintf('xmls/%s.xml',flog);
        f = fopen(fname,'w');
        while ~feof(g)
            line = fgets(g);
            if contains(line, 'insert_sequences')
                for j = 1 : length(leafnames)
                    fprintf(f,'\t\t<sequence id="%s" taxon="%s" totalcount="4" value="??"/>\n',leafnames{j},leafnames{j});
                end
            elseif contains(line, 'insert_init')        
                fprintf(f,'\t\t<init spec="beast.base.evolution.tree.TreeParser" id="NewickTree.t:Species" adjustTipHeights="false"\n');
                fprintf(f,'\t\t\tinitial="@Tree.t:H3N2" taxa="@H3N2"\n');
                fprintf(f,'\t\t\tIsLabelledNewick="true"\n');
                fprintf(f,'\t\t\tnewick="%s"/>\n',print_tree);
            elseif contains(line, 'insert_states')       
                tstates = 'rem';
                for j = 1 : length(leafnames)
                    tmp1 = strsplit(leafnames{j},'_');
                    tstates = [tstates sprintf(',%s=state%s',leafnames{j},tmp1{end})];
                end
                fprintf(f, strrep(line, 'insert_states', strrep(tstates, 'rem,','')));
            elseif contains(line, 'insert_mean')
                tmps = strsplit(flog,'_');
                tmps = strrep(tmps{2}, 'S','');
                fprintf(f, strrep(line, 'insert_mean', num2str(rate(end))));
            elseif contains(line ,'insert_rate_shifts')
                fprintf(f, strrep(line, 'insert_rate_shifts', sprintf('%f ', [1:25]/25)));
            else 
                fprintf(f, line);
            end
        end
        fclose(f);fclose(g);    
    end
end
