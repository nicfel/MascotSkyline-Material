clear

f = fopen('template/MERS_skyline.xml');

ids = cell(0,0);
seq = cell(0,0);

g = fopen('xmls/MERS_skyline.xml', 'w');

while ~feof(f)
    line = fgets(f);
    if contains(line, '<sequence')
        tmp = strsplit(line, '"');
        ids{end+1,1} = tmp{6};
        seq{end+1,1} = tmp{10};
        
    elseif contains(line, ' id="dateTrait')
        for i = 1 : length(ids)
            tmp1 = strsplit(ids{i}, '|');
            tmp2 = strsplit(tmp1{end}, '-');
            if length(tmp2) == 2
                line = strrep(line, [ids{i} '=' tmp1{end}], [ids{i} '=' tmp1{end} '-15']);                
            end
        end
    elseif contains(line, 'insert_tip_prior')
        for i = 1 : length(ids)
            tmp1 = strsplit(ids{i}, '|');
            tmp2 = strsplit(tmp1{end}, '-');
            if length(tmp2) == 2
                fprintf(g, '\t\t\t\t<distribution id="%s.leaf.prior" spec="beast.base.evolution.tree.MRCAPrior" monophyletic="true" tipsonly="true" tree="@Tree.t:MERS_CoV_274">\n',ids{i});
                fprintf(g, '\t\t\t\t\t<taxonset id="%s.leaf" spec="TaxonSet">\n',ids{i});
                fprintf(g, '\t\t\t\t\t\t<taxon id="%s" spec="Taxon"/>\n',ids{i});
                fprintf(g, '\t\t\t\t\t</taxonset>\n');
                from = datenum([tmp1{end} '-01']);
                to = addtodate(from, 1, 'month');
                fromy = datenum([tmp2{1} '-01-01']);
                toy = addtodate(fromy, 1, 'year');
                diffy = toy-fromy;
                
                from_dec = str2double(tmp2{1}) + (from-fromy)/diffy;
                to_dec = str2double(tmp2{1}) + (to-fromy)/diffy;

                fprintf(g, '\t\t\t\t\t<Uniform id="Uniform.%s" lower="%f" name="distr" upper="%f"/>\n',ids{i},from_dec, to_dec);
                fprintf(g, '\t\t\t\t</distribution>\n');
            end
        end
        line = fgets(f);
    elseif contains(line, 'insert_tip_ops')
        for i = 1 : length(ids)
            tmp1 = strsplit(ids{i}, '|');
            tmp2 = strsplit(tmp1{end}, '-');
            if length(tmp2) == 2
                fprintf(g, '\t\t<operator id="TipDatesRandomWalker.%s.leaf" spec="TipDatesRandomWalker" taxonset="@%s.leaf" tree="@Tree.t:MERS_CoV_274" weight="0.01" windowSize="0.01"/>\n', ids{i}, ids{i});
            end
        end
        line = fgets(f);
    elseif contains(line, 'insert_tip_log')
        for i = 1 : length(ids)
            tmp1 = strsplit(ids{i}, '|');
            tmp2 = strsplit(tmp1{end}, '-');
            if length(tmp2) == 2
                fprintf(g, '\t\t<log idref="%s.leaf.prior"/>\n', ids{i});
            end
        end
        line = fgets(f);

    end
    fprintf(g, line);
end

fclose('all');


%% check which sequences are from the same local camel outbreak



comb_name = cell(0,0);
for i = 1 : length(ids)
    if contains(ids{i}, 'camel')
        tmp = strsplit(ids{i}, '/');
        tmp2 = strsplit(ids{i}, '|');
        comb_name{i,1} = [tmp{1} '|' tmp2{end} '|' seq{i}];
    else
        comb_name{i,1} = 'NA';
    end
end


[u, ind]  = getUniqueStrings(comb_name);

%% sample which camel samples to keep

remove_seqs = cell(0,0);
for i=2 : max(ind)
    j=find(ind==i);
    if length(j)>1
       remove_seqs = [remove_seqs; randsample(ids(j), length(j)-1)];
    end       
end
%%



f = fopen('xmls/MERS_skyline.xml');

                
g = fopen('xmls/MERS_skyline_nolo.xml', 'w');
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

