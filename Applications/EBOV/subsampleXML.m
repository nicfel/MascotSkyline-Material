% removes some percentage of samples randomly
clear

rng(546546543);


system('rm -r xmls');
system('mkdir xmls');


fasta = fastaread('Makona_1610_genomes_2016-06-23.fasta');
use_nr_seqs = 300;

location = cell(0,0);
for i  = 1 : length(fasta)
    tmp = strsplit(fasta(i).Header, '|');
    location{i} = tmp{4};
end


for i = 1 : 2
f = fopen('EBOV.xml');
    if i==1
        rem_seq_indices = sort(randsample(length(fasta),length(fasta)-use_nr_seqs));

        remove_seq = cell(0,0);
        for j = 1 : length(rem_seq_indices)
            remove_seq{end+1,1} = fasta(rem_seq_indices(j)).Header;
        end

        g = fopen('xmls/ebov_proportional_skyline_rep0.xml', 'w');
    else
        sle = randsample(find(ismember(location,'SLE')), use_nr_seqs/3);
        gin = randsample(find(ismember(location,'GIN')), use_nr_seqs/3);
        lbr = randsample(find(ismember(location,'LBR')), use_nr_seqs/3);
        
        
        
        rem_seq_indices = 1:length(fasta);
        rem_seq_indices(ismember(rem_seq_indices,sle)) = [];
        rem_seq_indices(ismember(rem_seq_indices,gin)) = [];
        rem_seq_indices(ismember(rem_seq_indices,lbr)) = [];
        
        remove_seq = cell(0,0);
        for j = 1 : length(rem_seq_indices)
            remove_seq{end+1,1} = fasta(rem_seq_indices(j)).Header;
        end

        g = fopen('xmls/ebov_equal_skyline_rep0.xml', 'w');

    end

    while ~feof(f)
        line = fgets(f);
        if contains(line, 'id="seq_')
            tmp = strsplit(line, '"');
            if sum(ismember(remove_seq,tmp{6}))==0
                fprintf(g, line);
            end
        elseif contains(line, 'dateTrait')
            length(line)
            for i = 1 : length(remove_seq)
                tmp = strsplit(remove_seq{i}, '|');
                line = strrep(line, [remove_seq{i} '=' tmp{end}], '');
            end
            while contains(line, ',,')
                line = strrep(line, ',,',',');
            end
             line = strrep(line, '",','"');

            fprintf(g, line);
        elseif contains(line, 'typeTraitSet')
            length(line)
            for i = 1 : length(remove_seq)
                tmp = strsplit(remove_seq{i}, '|');
                if contains(remove_seq{i}, 'SLE')
                    line = strrep(line, [remove_seq{i} '=SLE'], '');
                elseif contains(remove_seq{i}, 'GIN')
                    line = strrep(line, [remove_seq{i} '=GIN'], '');
                else
                    line = strrep(line, [remove_seq{i} '=LBR'], '');
                end
           end
            while contains(line, ',,')
                line = strrep(line, ',,',',');
            end
             line = strrep(line, '",','"');
            fprintf(g, line);
        else
            fprintf(g, line);
        end
    end
    fclose(g);
end


%% make an unstructured file
for i = 1 :2
    if i==1
        f = fopen('./xmls/ebov_proportional_skyline_rep0.xml');
        g = fopen('./xmls/ebov_proportional_unstructured_rep0.xml', 'w');
    else
        f = fopen('./xmls/ebov_equal_skyline_rep0.xml');
        g = fopen('./xmls/ebov_equal_unstructured_rep0.xml', 'w');
    end
    while ~feof(f)
        line = fgets(f);
        if contains(line, 'SLE GIN LBR') 
            fprintf(g, strrep(line, 'SLE GIN LBR', 'WA dummy'));
        elseif contains(line, '=SLE') || contains(line, '=GIN') || contains(line, '=LBR')
            line = strrep(line, '=SLE','=WA');
            line = strrep(line, '=GIN','=WA');
            line = strrep(line, '=LBR','=WA');
            fprintf(g, line);
        elseif contains(line, '<NeDynamics id="skygrid.state2"')

        elseif contains(line, 'f="@NeLog.state1"') || contains(line, 'f="@NeLog.state2"')

        elseif contains(line, 'RW2') || contains(line, 'RW3')

        elseif contains(line, 'Prior'' NeLog="@NeLog.state1"')
            while ~contains(line, '</distribution>')
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
end

%% make replicates
clear dir
files = dir('xmls/*rep0.xml');

for i = 1:length(files)
    name = ['./xmls/' files(i).name];
    for j =1:2
        f1 = fopen(name);

        g1 = fopen(strrep(name, 'rep0',['rep' num2str(j)]), 'w');

        while ~feof(f1)
            fprintf(g1, fgets(f1));
        end
    end
end
