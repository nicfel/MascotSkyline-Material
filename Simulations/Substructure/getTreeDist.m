clear
% script to compute the distance of the reconstructed tree to the true tree
% using the simulated trees in /master/

% first, get all the mascot trees that end in events.trees in ../out
mascot_trees = dir('../out/SIR*events.trees');
% loop over all trees
for i = 1 : length(mascot_trees)
    % get the true tree corresponding to this run by getting the S%d number from the name
    S = strsplit(mascot_trees(i).name,'_');
    S = S{2};
    % read in the true tree line by line
    true_tree = fopen(sprintf('./master/SIR_%s_master.tree',S));
    true_tree = textscan(true_tree,'%s','Delimiter','\n');
    % get the line that starts with tree TREE
    true_tree = true_tree{1}{cellfun(@(x) ~isempty(strfind(x,'tree TREE')),true_tree{1})};
    % get all the locations="%d", by getting the %d value
    true_locations = regexp(true_tree,'(?<=location=")\d+(?=")','match');
    % convert to numbers
    true_locations = cellfun(@str2num,true_locations);
    % remive tree TREE = from the beginning
    true_tree = regexprep(true_tree,'tree TREE = ','');
    % label all nodes based on the number of )
    true_tree = regexprep(true_tree,'\)',')n_');
    % replace everything between [ and ] with the corresponding true_locations value
    true_tree = sprintf(regexprep(true_tree,'\[.*?\]','[%d]'),true_locations);
    % replace all instances of %d[%d], with the first number only
    true_tree = regexprep(true_tree,'(?<=\d)\[\d+\]','');
    % replace all instances of [ or ] with nothing
    true_tree = regexprep(true_tree,'[\[\]]','');    
    % replace all n with n%d, use %d as the number of n's
    true_tree = sprintf(regexprep(true_tree,'n','n%d'),1:length(strfind(true_tree,'n')));
    % read in the true_tree and define all the locations based on the leaf nodes below that node
    true_tree = phytreeread(true_tree);
    % get the  nodes
    nodenames = get(true_tree,'NodeNames');
    % for each nodenames, get all the leaf nodes below that node and save the ordered list of locations
    true_corresponding_leafs = cell(floor(length(nodenames)/2),1);
    true_location = zeros(floor(length(nodenames)/2),1);
    for j = ceil(length(nodenames)/2)+1 : length(nodenames)
        % get the leaf nodes below this node
        st = subtree(true_tree,j);
        %get all leafs in the subtree
        leafs = get(st,'LeafNames');
        % get the number of the node, by getting the first number after n in n%d_%d in nodenames
        node_number = regexp(nodenames{j},'(?<=n)\d+','match');
        % get the second number to get the actual location
        node_loc = regexp(nodenames{j},'\d+','match');
        % convert to number
        node_number = str2num(node_number{1});  
        if node_number>500
            dsa
        end
        % join the names of all leafs seperated by , after ordering them
        true_corresponding_leafs{node_number} = strjoin(sort(leafs),',');
        true_location(node_number) = str2num(node_loc{length(node_loc)});
    end

    % read in the mascot tree line by line
    mascot_tree = fopen(sprintf('../out/%s',mascot_trees(i).name));
    mascot_tree = textscan(mascot_tree,'%s','Delimiter','\n');
    % keeps track of if this is the first tree analysed
    first_tree = 1;

    diff_tree = zeros(0,1);

    % loop over all lines in the mascot tree
    for j = 1 : length(mascot_tree{1})
        % if the line starts with tree TREE
        if ~isempty(strfind(mascot_tree{1}{j},'tree STATE_'))
            % get all blocks of )[...] and get the state%d value within each instance of )[...]
            mascot_locations = regexp(mascot_tree{1}{j},'(?<=\)\[).*?(?=\])','match');
            % keep the last char in each block and convert it to a number
            mascot_locations = cellfun(@(x) str2num(x(end)),mascot_locations);
            if first_tree
                % do the same as for the true tree, but here for the mascot tree
                mt = mascot_tree{1}{j};
                mt = regexprep(mt,'tree STATE_(\d*) = ','');
                mt = regexprep(mt,'\)',')n_');
                % get all the "state%d," by getting the %d value
                ml = regexp(mascot_tree{1}{j},'(?<=state)\d+(?=)','match');
                % convert to numbers
                ml = cellfun(@str2num,ml);
                mt = sprintf(regexprep(mt,'\[.*?\]','[%d]'),ml);
                mt = regexprep(mt,'(?<=\d)\[\d+\]','');
                mt = regexprep(mt,'[\[\]]','');
                mt = sprintf(regexprep(mt,'n','n%d'),1:length(strfind(mt,'n')));
                mt = phytreeread(mt);
                nodenames = get(mt,'NodeNames');
                corresponding_leafs = cell(floor(length(nodenames)/2),1);
                for k = ceil(length(nodenames)/2)+1 : length(nodenames)
                    st = subtree(mt,k);
                    leafs = get(st,'LeafNames');
                    node_number = regexp(nodenames{k},'(?<=n)\d+','match');
                    node_number = str2num(node_number{1});
                    corresponding_leafs{node_number} = strjoin(sort(leafs),',');
                    if length(leafs)==0
                        dsa
                    end
                end

                map = zeros(floor(length(nodenames)/2),1);
                for k = 1 : length(map)
                    map(k) = find(strcmp(corresponding_leafs{k},true_corresponding_leafs));
                end               
            end
            % remove everything between [ and ] from the mascot tree
            mascot_tree{1}{j} = regexprep(mascot_tree{1}{j},'\[.*?\]','');
            %remove all numbers after : includeing E- and E+
            mascot_tree{1}{j} = regexprep(mascot_tree{1}{j},':\d+.\d+E[+-]\d+','');
            disp(mascot_tree{1}{j})

            
            % compare how often the true and mascot node states diffe
            diff_tree(end+1) = sum(true_location(map) ~= mascot_locations');
        end
    end
end
