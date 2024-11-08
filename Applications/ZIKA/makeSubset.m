f = fopen('xmls/ZIKV_skyline_rep0.xml');
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
uni_locs = unique(loc);

keep_locs = uni_locs([2,3,4,6,7]);

system('rm -r redxml');
system('mkdir redxml');
rng(45568);
in_ne = seq(ismember(loc, 'Brazil_Northeast'));
date_ne = zeros(0,0);
for i = 1 : length(in_ne)
    tmp = strsplit(in_ne{i}, '|');
    date_ne(i) = str2double(tmp{end});
    
end


[a,b] = sort(date_ne);

% [~,b] = sort(date_ne);
% in_ne = in_ne(b);
% randomize order
for nr_seq = 1 : 1 : 10
    

    f = fopen('xmls/ZIKV_skyline_rep0.xml');
    remove_seqs = seq(~ismember(loc, keep_locs));
    for i=1:nr_seq
        remove_seqs = [remove_seqs; in_ne(b==i)];
    end
    g = fopen(sprintf('redxml/mascot_zikv_red.%d.xml',nr_seq), 'w');
    first = true;

    while ~feof(f)
        line = fgets(f);
        if contains(line,'sequence')
            tmp = strsplit(line, '"');            
            if ~ismember(tmp{6}, remove_seqs)
                fprintf(g, line);
            end
        elseif contains(line,'ZIKV|') || contains(line,'ZBR')      

            for i = 1 : length(remove_seqs)
                line =strrep(line, remove_seqs{i}, '');
            end
            if first

                for i = 1 : length(remove_seqs)
                    line =regexprep(line, ',=(\d*)\.(\d*)', '');
                    line =regexprep(line, '"=(\d*)\.(\d*)', '"');
                end
                for i = 1 : length(remove_seqs)
                    line =regexprep(line, ',=(\d*)', '');
                    line =regexprep(line, '"=(\d*)', '"');
               end
            end

            for i = 1 : length(remove_seqs)
                line =regexprep(line, ',=(\w*)_(\w*)', '');
                line =regexprep(line, '"=(\w*)_(\w*)', '"');
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
            sequence{end+1,1} = tmp{10}(1:length(tmp{10}/10));
        elseif contains(line,'ZIKV|') || contains(line,'ZBR')
            tmp = strsplit(line, '"');
            
            tmp2 = strsplit(tmp{8}, ',');
            for i = 1: length(tmp2)
                tmp3 = strsplit(tmp2{i}, '=');
                ind = find(ismember(id,tmp3{1}));              
                
                if contains(tmp3{2}, '.')
                    date{ind,1} = tmp3{2};                    
                else
                    state{ind,1} = tmp3{2};
                end
            end
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
                fprintf(g,'\t\t<taxon id="%s">\n', id{i});
                fprintf(g,'\t\t\t<date value="%s" direction="forwards" units="years"/>\n', date{i});
                fprintf(g,'\t\t\t<attr name="division">\n');
                fprintf(g,'\t\t\t\t%s\n',state{i});
                fprintf(g,'\t\t\t</attr>\n');
                fprintf(g,'\t\t</taxon>\n');
            end
        else
            fprintf(g, line);
        end
    end
        
end

