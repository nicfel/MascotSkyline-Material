clear
f = fopen('ZIKV.ZiBRA_no3Rio-Pardis_noFlorida-4_Mexico_Chiu.allGenBank.basta.xml');
c = 1;
id = cell(0,0);
while ~feof(f)
    line = fgets(f);
    if contains(line, '<sequence')
        if contains(line, 'Micronesia') ||...
                contains(line, 'Philippines') ||...
                contains(line, 'Thailand') || contains(line, 'Cambodia')
            disp(line)
        else
            tmp = strsplit(line, '"');
            tmp2 = strsplit(tmp{4}, '|');    

            Data(c).Header = [tmp2{1} '|' tmp2{2} '|' tmp2{3}];
            Data(c).Sequence = tmp{6};
            id{c} = tmp{4};
            c = c+1;
        end
    elseif contains(line, id)
        tmp = strsplit(strtrim(line), '=');
        ind = find(ismember(id,tmp{1}));
        Data(ind).Header = [Data(ind).Header '|' strrep(tmp{2}, ',','')];        
    end        
end
fclose(f)
delete('basta_alignment.fasta');
fastawrite('basta_alignment.fasta', Data)