
t = readtable("nextstrain_groups_blab_ncov_wa-phylodynamics_metadata.tsv", "FileType","text",'Delimiter', '\t');
% epi_isl = randsample(t.gisaid_epi_isl, 1500);
% f = fopen('gisaid.tsv', 'w');
% 
% for i = 1:length(epi_isl)
%     fprintf(f, '%s\n', epi_isl{i});
% end
% 
% fclose('all')


% t = readtable("nextstrain_groups_blab_ncov_wa-phylodynamics_metadata.tsv", "FileType","text",'Delimiter', '\t');


fasta = fastaread('gisaid_hcov-19_2023_08_16_17.fasta');
ids=cell(0,0);
region=cell(0,0);
for i = 1:length(fasta)
    tmp = strsplit(fasta(i).Header, '|');
    ids{i} = tmp{2};
end

f = fopen('ncov_wa-phylodynamics.json');

vals = zeros(0,0);

while ~feof(f)
    line = fgets(f);
    % if the line contains any of the ids in the cell ids, then...
    if any(contains(line, ids))
        % get the index in ids of the id that matches the line
        name = strsplit(line, '"');
        idx = find(contains(ids, name{4}));
        vals(end+1) = idx;
        disp(length(vals));
        % check the next lines until you find the one that contains the attribute region or county
        c = 0;
        while ~feof(f)
            line = fgets(f);
            if contains(line, '"region"')
                line = fgets(f);
                tmp = strsplit(line, '"');
                region{idx} = tmp{4};
                c=c+1;
            end
            if contains(line, '"county"')
                line = fgets(f);
                tmp = strsplit(line, '"');
                county{idx} = tmp{4};
                c=c+1;                
            end
            if c==2 || contains(line, '"url"')
                break;
            end
        end
    end   
end

%%
nf = fasta;
remove = zeros(0,0);
for i = 1:length(nf)
    tmp = strsplit(nf(i).Header, '|');
    if strcmp(region{i}, 'Washington State')   
        if i>length(county)|| isempty(county{i}) || contains(county{i}, '2')
            remove(end+1) = i;
        else
            nf(i).Header = strcat(nf(i).Header, '|', strrep(county{i}, ' ', ''));
        end
    else    
        if strcmp(region{i}, 'North America')
            nf(i).Header = strcat(nf(i).Header, '|', 'NorthAmerica');
        else
            nf(i).Header = strcat(nf(i).Header, '|', 'restOfWorld');
        end
    end
end

nf(remove) = [];
multialign(nf, 'ScoringMatrix', 'NUC44', 'UseParallel', 'true', 'Verbose', 'true')

%%

delete('dataset.fasta');
% save nf without the remove indices
fastawrite('dataset.fasta', nf)

