% removes some percentage of samples randomly
clear

rng(546546543);

fasta = fastaread('Makona_1610_genomes_2016-06-23.fasta');
use_nr_seqs = 350;
rem_seq_indices = sort(randsample(length(fasta),length(fasta)-use_nr_seqs));
remove_seq = cell(0,0);
for i = 1 : length(rem_seq_indices)
    remove_seq{end+1,1} = fasta(rem_seq_indices(i)).Header;
end

system('rm -r xmls');
system('mkdir xmls');

f = fopen('EBOV.xml');
g = fopen('xmls/ebov_skyline_rep0.xml', 'w');

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


%% make an unstructured file
f = fopen('./xmls/ebov_skyline_rep0.xml');
g = fopen('./xmls/ebov_unstructured_rep0.xml', 'w');
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

%% make replicates
for i = 1:2
    f1 = fopen('./xmls/ebov_skyline_rep0.xml');
    f2 = fopen('./xmls/ebov_unstructured_rep0.xml');
    
    g1 = fopen(sprintf('./xmls/ebov_skyline_rep%d.xml', i), 'w');
    g2 = fopen(sprintf('./xmls/ebov_unstructured_rep%d.xml', i), 'w');
    
    while ~feof(f1)
        fprintf(g1, fgets(f1));
    end
    
    while ~feof(f2)
        fprintf(g2, fgets(f2));
    end


end
