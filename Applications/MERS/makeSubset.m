f = fopen('xmls/MERS_skyline.xml');
% g = fopen('zikv_red.xml', 'w');
seq = cell(0,0);
loc = cell(0,0);
while ~feof(f)
    line = fgets(f);
    if contains(line,'<typeTrait')
        disp(line)
        tmp = strsplit(line, '"');
        tmp2 = strsplit(tmp{8}, ',');
        for i = 1: length(tmp2)
            tmp3 = strsplit(tmp2{i}, '=');
            seq{i,1} = tmp3{1};
            loc{i,1} = tmp3{2};
        end
    end
end
fclose('all')
% get unique locations
uni_locs = unique(loc);
keep_locs = uni_locs;
system('rm -r redxml');
system('mkdir redxml');
rng(45568);

in_ne = seq(ismember(loc, 'human'));
date_ne = zeros(0,0);
for i = 1 : length(in_ne)
    tmp = strsplit(in_ne{i}, '|');
    date_ne(i) = str2double(tmp{end});    
end

% randomize order
for repetition = 10 : 30 : 190    
    % open template
    f = fopen('xmls/MERS_skyline.xml');
    
    if nr_seq <= 100
        per(1) = nr_seq;
        per(2) = 100;
    else
        per(1) = 100;
        per(2) = abs(nr_seq-200);        
    end

    r_seqs = zeros(0,0);
    for i = 1 : length(keep_locs)
        in_loc = ismember(loc, keep_locs{i});
        nr_seqs=floor(sum(in_loc)*(100-per(i))/100);   
        r_seqs=[r_seqs; randsample(find(in_loc), nr_seqs)];
        
    end
    
    remove_seqs = seq(sort(r_seqs));
                
    g = fopen(sprintf('redxml/mascot_mers_red.%d.xml',nr_seq), 'w');
    first = true;

    while ~feof(f)
        line = fgets(f);
        if contains(line,'sequence')
            tmp = strsplit(line, '"');            
            if ~ismember(tmp{6}, remove_seqs)
                fprintf(g, line);
            end
        elseif contains(line, 'leaf.prior"/>')
            tmp = strsplit(line, '"');
            tmp2 = strrep(tmp{2}, '.leaf.prior', '');
            if ~ismember(tmp2, remove_seqs)
                fprintf(g, line);
            end           
        elseif contains(line,'leaf.prior') 
            tmp = strsplit(line, '"');
            tmp2 = strrep(tmp{2}, '.leaf.prior', '');
            if ismember(tmp2, remove_seqs)
                fgets(f);fgets(f);fgets(f);fgets(f);fgets(f);
            else
                fprintf(g, line);
            end
        elseif contains(line,'TipDatesRandomWalker') 
            tmp = strsplit(line, '"');
            tmp2 = strrep(tmp{2}, 'TipDatesRandomWalker.', '');
            tmp2 = strrep(tmp2, '.leaf', '');
            if ~ismember(tmp2, remove_seqs)
                fprintf(g, line);
            end
        elseif contains(line,'|camel|') || contains(line,'|human|')      

            for i = 1 : length(remove_seqs)
                line = strrep(line, remove_seqs{i}, '');
            end
            if first
                for i = 1 : length(remove_seqs)
                    line =regexprep(line, ',=(\d*)-(\d*)-(\d*)', '');
                    line =regexprep(line, '"=(\d*)-(\d*)-(\d*)', '"');
                end
            end

            for i = 1 : length(remove_seqs)
                line =regexprep(line, ',=(\w*)', '');
                line =regexprep(line, '"=(\w*)', '"');
            end
            
            for i = 1 : length(remove_seqs)
                line =regexprep(line, '",', '"');
            end

            fprintf(g, line);
            
            first = false;
        else
            fprintf(g, line);
        end
    end
    fclose('all');
end


% randomize order
for nr_seq = 10 : 30 : 190    

    f = fopen('xmls/MERS_skyline_nolo.xml');
    
    if nr_seq <= 100
        per(1) = nr_seq;
        per(2)=100;
    else
        per(1)=100;
        per(2)=abs(nr_seq-200);        
    end
    r_seqs = zeros(0,0);
    for i = 1 : length(keep_locs)
        in_loc = ismember(loc, keep_locs{i});
        nr_seqs=floor(sum(in_loc)*(100-per(i))/100);   
        r_seqs=[r_seqs; randsample(find(in_loc), nr_seqs)];
        
    end
    
    remove_seqs = seq(sort(r_seqs));
                
    g = fopen(sprintf('redxml/mascotnolo_mers_red.%d.xml',nr_seq), 'w');
    first = true;

    while ~feof(f)
        line = fgets(f);
        if contains(line,'sequence')
            tmp = strsplit(line, '"');            
            if ~ismember(tmp{6}, remove_seqs)
                fprintf(g, line);
            end
        elseif contains(line, 'leaf.prior"/>')
            tmp = strsplit(line, '"');
            tmp2 = strrep(tmp{2}, '.leaf.prior', '');
            if ~ismember(tmp2, remove_seqs)
                fprintf(g, line);
            end           
        elseif contains(line,'leaf.prior') 
            tmp = strsplit(line, '"');
            tmp2 = strrep(tmp{2}, '.leaf.prior', '');
            if ismember(tmp2, remove_seqs)
                fgets(f);fgets(f);fgets(f);fgets(f);fgets(f);
            else
                fprintf(g, line);
            end
        elseif contains(line,'TipDatesRandomWalker') 
            tmp = strsplit(line, '"');
            tmp2 = strrep(tmp{2}, 'TipDatesRandomWalker.', '');
            tmp2 = strrep(tmp2, '.leaf', '');
            if ~ismember(tmp2, remove_seqs)
                fprintf(g, line);
            end
        elseif contains(line,'|camel|') || contains(line,'|human|')      

            for i = 1 : length(remove_seqs)
                line = strrep(line, remove_seqs{i}, '');
            end
            if first
                for i = 1 : length(remove_seqs)
                    line =regexprep(line, ',=(\d*)-(\d*)-(\d*)', '');
                    line =regexprep(line, '"=(\d*)-(\d*)-(\d*)', '"');
                end
            end

            for i = 1 : length(remove_seqs)
                line =regexprep(line, ',=(\w*)', '');
                line =regexprep(line, '"=(\w*)', '"');
            end
            
            for i = 1 : length(remove_seqs)
                line =regexprep(line, '",', '"');
            end

            fprintf(g, line);
            
            first = false;
        else
            fprintf(g, line);
        end
    end
    fclose('all');
end



%% convert xmls to dta
xmls = dir('redxml/mascot*.xml');
for r = 1 : length(xmls)
    f = fopen(['redxml/' xmls(r).name]);
    id = cell(0,0);
    sequence = cell(0,0);
    state = cell(0,0);
    date = cell(0,0);
    while ~feof(f)
        line = fgets(f);
        if contains(line,'sequence')
            tmp = strsplit(line, '"');
            id{end+1,1} = tmp{6};
            sequence{end+1,1} = tmp{10};                     
        end
    end
    fclose('all');
    
    
    f = fopen('dta_template.xml');
    g = fopen(['redxml/' strrep(xmls(r).name, 'mascot','dta')], 'w');
    while ~feof(f)
        line = fgets(f);
        if contains(line, 'insert_seqs')
            for i = 1 : length(id)
                fprintf(g,'\t\t<sequence>\n');
                fprintf(g,'\t\t\t<taxon idref="%s"/>\n', id{i});
                fprintf(g,'\t\t\t%s\n', sequence{i});
                fprintf(g,'\t\t</sequence>\n');
            end

            
        elseif contains(line, 'insert_dat')
            for i = 1 : length(id)
                tmp = strsplit(id{i}, '|', 'CollapseDelimiters',false);
                
                tmp2 = strsplit(tmp{4}, '-');
                
                if length(tmp2)==2
                    tmp{4} = [tmp{4} '-15'];                    
                end
                tmp2 = strsplit(tmp{4}, '-');
                
                
                
                fprintf(g,'\t\t<taxon id="%s">\n', id{i});
                fprintf(g,'\t\t\t<date value="%s/%s/%s" direction="forwards" units="years"/>\n', tmp2{3}, tmp2{2}, tmp2{1});
                fprintf(g,'\t\t\t<attr name="division">\n');
                fprintf(g,'\t\t\t\t%s\n',tmp{3});
                fprintf(g,'\t\t\t</attr>\n');
                fprintf(g,'\t\t</taxon>\n');
            end
        
        elseif contains(line, 'insert_filename')
            fname = strrep(xmls(r).name, 'mascot','dta');
            fprintf(g, strrep(line, 'insert_filename', strrep(fname, '.xml','')));
        else
            fprintf(g, line);
        end
    end
        
end


% %% convert xmls to add "other" deme
% xmls = dir('redxml/mascot*.xml');
% for r = 1 : length(xmls)
%     f = fopen(['redxml/' xmls(r).name]);
%     g = fopen(['redxml/' strrep(xmls(r).name, 'mascot','mascotother')], 'w');
%     while ~feof(f)
%         line = fgets(f);
%         if contains(line, 'insert_seqs')
%             for i = 1 : length(id)
%                 fprintf(g,'\t\t<sequence>\n');
%                 fprintf(g,'\t\t\t<taxon idref="%s"/>\n', id{i});
%                 fprintf(g,'\t\t\t%s\n', sequence{i});
%                 fprintf(g,'\t\t</sequence>\n');
%             end
%         elseif contains(line, '<log idref="SkylineNe.human.t:MERS_CoV_274"/>')
%             fprintf(g, line);
%             fprintf(g, '\t\t<log idref="SkylineNe.other.t:MERS_CoV_274"/>\n');        
%         elseif contains(line, 'insert_filename')
%             fname = strrep(xmls(r).name, 'mascot','dta');
%             fprintf(g, strrep(line, 'insert_filename', strrep(fname, '.xml','')));
%         else
%             fprintf(g, line);
%         end
%     end
        
% end

