%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% creates structured coalescent xmls from the master trees. Always creates
% 3 xmls per tree with different initial values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
% create lisco files
tree_files = dir('master/*.tree');

system('rm -r dta');
system('mkdir dta');

states = 3;

for i = 1 : length(tree_files)
    if contains(tree_files(i).name, 'single')
        continue
    end
    disp(tree_files(i).name)
    % read tree files
    g = fopen(['master/' tree_files(i).name],'r');
    t = textscan(g,'%s'); fclose(g);
    
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
        g = fopen('dta_template.xml');
        flog = ['DTA' strrep(tree_files(i).name,'master.tree',sprintf('%ddta',tr))];
        fname = sprintf('dta/%s.xml',flog);
        f = fopen(fname,'w');
        while ~feof(g)
            line = fgets(g);
            if contains(line, 'insert_seqs')
                for j = 1 : length(leafnames)
                    fprintf(f,'\t\t\t<sequence>\n');
                    fprintf(f,'\t\t\t\t<taxon idref="%s"/>\n', leafnames{j});
                    fprintf(f,'\t\t\t\tNN\n');
                    fprintf(f,'\t\t\t</sequence>\n');
                    
%                     fprintf(f,'\t\t<sequence id="%s" taxon="%s" totalcount="4" value="??"/>\n',leafnames{j},leafnames{j});
                end
            elseif contains(line, 'insert_starting_tree')    
                fprintf(f,'\t\t\t%s\n', print_tree);
%                 fprintf(f,'\t\t<init spec="beast.base.evolution.tree.TreeParser" id="NewickTree.t:Species" adjustTipHeights="false"\n');
%                 fprintf(f,'\t\t\tinitial="@Tree.t:H3N2" taxa="@H3N2"\n');
%                 fprintf(f,'\t\t\tIsLabelledNewick="true"\n');
%                 fprintf(f,'\t\t\tnewick="%s"/>\n',print_tree);
            elseif contains(line, 'insert_dat')                
                for j = 1 : length(leafnames)
                    tmp1 = strsplit(leafnames{j},'_');
                    fprintf(f,'\t\t<taxon id="%s">\n', leafnames{j});
                    fprintf(f,'\t\t\t<date value="%s" direction="forwards" units="years"/>', tmp1{end});
                    fprintf(f,'\t\t\t<attr name="division">\n');
                    fprintf(f,'\t\t\t\tstate%s\n', tmp1{2});
                    fprintf(f,'\t\t\t</attr>\n');
                    fprintf(f,'\t\t</taxon>\n');
                end
            elseif contains(line, 'insert_mean')
                tmps = strsplit(flog,'_');
                tmps = strrep(tmps{2}, 'S','');
                if str2double(tmps)>50
                    fprintf(f, strrep(line, 'insert_mean', '5.0'));
                else
                    fprintf(f, strrep(line, 'insert_mean', '1.0'));
                end

            elseif contains(line ,'$(filename)')
                fprintf(f, strrep(line, '$(filename)', flog));
            else 
                fprintf(f, line);
            end
        end
        fclose(f);fclose(g);    
    end
end
