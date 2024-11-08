clear
fasta = fastaread('alignment.fasta');

for i = 1 : length(fasta)
    tmp = strsplit(fasta(i).Header, '|');
    location{i} = tmp{length(tmp)-1};
end

uni_loc = unique(location)';