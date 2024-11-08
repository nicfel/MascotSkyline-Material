% removes some percentage of samples randomly
clear

system('rm -r xmls');
system('mkdir xmls');

fasta = fastaread('data/GenomicFastaResults.fasta');
clear all_sequences
c = 1;
location = cell(0,0);
year = zeros(0,0);
for i  = 1 : length(fasta)
    tmp = strsplit(fasta(i).Header, '|');
    if length(tmp)==4
        tmp2=strsplit(tmp{3},'/');
        if length(tmp2)==3
            all_sequences(c) = fasta(i);
            location{c} = tmp{4};
            year(c) = str2double(tmp2{3});
            c=c+1;
        end
    end
end

rng(546546543);

is_nz = randsample(find(ismember(location, {'New_Zealand', 'Australia'})),200);
is_notnz = find(~ismember(location, {'New_Zealand', 'Australia'}));


% sample 100 sequences randomly and weighted by first year, then location
use_notnzsample = zeros(0,0);
for i = 1 : 50
    target_year = randi([2000 2005],1);
    disp(target_year)
    indices = is_notnz(find(year(is_notnz)==target_year));
    uni_locs = unique(location(indices));
    target = indices(ismember(location(indices), randsample(uni_locs,1)));
    if length(target)==1
        use_notnzsample(end+1,1) = target;
    else
        use_notnzsample(end+1,1) = randsample(target,1);
    end
    is_notnz(is_notnz==use_notnzsample(end)) = [];    
end

use_notnzsample = randsample(use_notnzsample, length(use_notnzsample));

for not_nz_sample = 0:50:50

    if not_nz_sample>0
        use = [is_nz use_notnzsample(1:not_nz_sample)'];
    else
        use = is_nz;
    end

    sequences = all_sequences(use);

    f = fopen('template.xml');
    g = fopen(sprintf('xmls/h3n2_structured_%d_rep0.xml',not_nz_sample), 'w');

    while ~feof(f)
        line = fgets(f);
        if contains(line, 'insert_sequence')
            for i = 1 : length(sequences)
                fprintf(g, '\t\t<sequence id="seq_%s" spec="Sequence" taxon="%s" totalcount="4" value="%s"/>\n', sequences(i).Header, sequences(i).Header, sequences(i).Sequence);
            end
        elseif contains(line, 'insert_date')
            date = 'rem';
            for i = 1 : length(sequences)
                tmp = strsplit(sequences(i).Header, '|');
                date = [date ',' sequences(i).Header  '=' tmp{3}];
            end

            fprintf(g, strrep(line, 'insert_date', strrep(date, 'rem,','')));
        elseif contains(line, 'insert_type')
            date = 'rem';
            for i = 1 : length(sequences)
                tmp = strsplit(sequences(i).Header, '|');
                if contains(sequences(i).Header, 'New_Zealand') || contains(sequences(i).Header, 'Australia')
                    date = [date ',' sequences(i).Header  '=Oceania'];
                else
                    date = [date ',' sequences(i).Header  '=Other'];
                end
            end

            fprintf(g, strrep(line, 'insert_type', strrep(date, 'rem,','')));
        else
            fprintf(g, line);
        end
    end
    fclose(g);
    fclose(f);
    
    
    f = fopen('dta_template.xml');
    g = fopen(sprintf('xmls/h3n2_dta_%d_rep0.xml',not_nz_sample), 'w');
    while ~feof(f)
        line = fgets(f);
        if contains(line, 'insert_seqs')
            for i = 1 : length(sequences)
                fprintf(g,'\t\t<sequence>\n');
                fprintf(g,'\t\t\t<taxon idref="%s"/>\n', sequences(i).Header);
                fprintf(g,'\t\t\t%s\n', sequences(i).Sequence);
                fprintf(g,'\t\t</sequence>\n');
            end

            
        elseif contains(line, 'insert_dat')
            for i = 1 : length(sequences)
                fprintf(g,'\t\t<taxon id="%s">\n', sequences(i).Header);
                tmp = strsplit(sequences(i).Header, '|');
                tmp2 = strsplit(tmp{3}, '/');
                date = (datenum(tmp{3},'mm/dd/yyyy')-datenum(tmp2{3},'yyyy'))/...
                    (datenum(num2str(str2double(tmp2{3})+1),'yyyy')-datenum(tmp2{3},'yyyy'));
                date = date+str2double(tmp2{3});                
                fprintf(g,'\t\t\t<date value="%.5f" direction="forwards" units="years"/>\n', date);
                fprintf(g,'\t\t\t<attr name="division">\n');
                if contains(sequences(i).Header, 'New_Zealand') || contains(sequences(i).Header, 'Australia')
                    fprintf(g,'\t\t\t\tNewZealand\n');
                else
                    fprintf(g,'\t\t\t\tother\n');
                end
                fprintf(g,'\t\t\t</attr>\n');
                fprintf(g,'\t\t</taxon>\n');
            end
        else
            fprintf(g, line);
        end
    end

end


%% convert xmls to dta
xmls = dir('xmls/h3n2_structured*.xml');
for r = 1 : length(xmls)
    f = fopen(['xmls/' xmls(r).name]);
    id = cell(0,0);
    sequence = cell(0,0);
    state = cell(0,0);
    date = cell(0,0);
    while ~feof(f)
        line = fgets(f);
        if contains(line,'sequence')
            tmp = strsplit(line, '"');
            id{end+1,1} = tmp{4};
            sequence{end+1,1} = strrep(tmp{6}, '/>', '');
        elseif contains(line,'ZIKV|') || contains(line,'ZBR')
            tmp = strsplit(line, '"');
            if contains(line, '>')
                tmp = strtrim(strrep(tmp{end-1}, ',',''));
            else
                tmp = strtrim(strrep(tmp{end}, ',',''));
            end
            tmp = strtrim(strrep(tmp, '">',''));
            tmp = strsplit(tmp, '=');
            if length(tmp)==2
                ind = find(ismember(id,tmp{1}));
                if contains(tmp{2}, '.')
                    date{ind,1} = tmp{2};
                else
                    state{ind,1} = tmp{2};
                end
            end
        end
    end
    fclose('all');
        
end


%% make an unstructured file
xmlfiles = dir('xmls/*.xml');
for i = 1 : length(xmlfiles)
    if contains(xmlfiles(i).name, 'h3n2_structured_0_rep0.xml')
        f = fopen(['xmls/' xmlfiles(i).name]);
        g = fopen(['xmls/' strrep(xmlfiles(i).name, '_structured','_unstructured')], 'w');
        location = strsplit(xmlfiles(i).name,'_');
        while ~feof(f)
            line = fgets(f);
            if contains(line, 'id="migrationIndicator"')
                fprintf(g, strrep(line,'>true<', '>false<'));

            elseif contains(line, 'f="@NeLog.state1"') || contains(line, 'f="@NeLog.state2"') || contains(line,'f="@sigma.val2"') || contains(line, 'LogTransform" f="@migrationConstant.t:EBOV_1')

            elseif contains(line, 'NeScaler2')
                
            elseif contains(line, 'NeScaler2')


            elseif contains(line, 'Prior'' NeLog="@NeLog.state1"')
                while ~contains(line, '</distribution>')
                    line = fgets(f);
                end
            elseif contains(line, 'NeSwapper')
                while ~contains(line, '</operator>')
                    line = fgets(f);
                end
            elseif contains(line, 'Prior''  NeLog="@NeLog.state2"')
                while ~contains(line, '</distribution>')
                    line = fgets(f);
                end

            else
                fprintf(g, line);
            end
        end
        fclose('all');
    end
end

%% make replicates
clear dir
files = dir('xmls/*rep0.xml');

for i = 1:length(files)
    name = ['./xmls/' files(i).name];
    for j = 1 : 2
        f1 = fopen(name);

        g1 = fopen(strrep(name, 'rep0',['rep' num2str(j)]), 'w');

        while ~feof(f1)
            fprintf(g1, fgets(f1));
        end
    end
end
