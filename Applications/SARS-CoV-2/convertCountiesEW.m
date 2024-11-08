rng(6541589)
% read in aligned data
fasta = fastaread('dataset_aligned.fasta');

% vector to keep track of which sequences to keep
keep = false(1, length(fasta));
% loop over all sequences, splitting the header by | and checking if the third group has two -, otherwise discard
for i = 1:length(fasta)
    header = strsplit(fasta(i).Header, '|');
    if length(strsplit(header{3}, '-')) == 3
        disp(header{4})
        % if header{4} says NorthAmerica or restOfWorld, keep it with a probability of 0.3
        if strcmp(header{4}, 'NorthAmerica') || strcmp(header{4}, 'restOfWorld')
            if rand < 0.1
                keep(i) = true;
            else
                keep(i) = false;
            end
        else
            keep(i) = true;
        end    
    end
end

% intiaialize a vector that categorizes the counties in Washington State into western and eastern Washington
west = {'Clallam', 'Clark', 'Cowlitz', 'GraysHarbor', 'Island', 'Jefferson', 'King', 'Kitsap', 'Lewis', 'Mason', 'Pacific', 'Pierce', 'San Juan', 'Skagit', 'Skamania', 'Snohomish', 'Thurston', 'Wahkiakum', 'Whatcom'};
east = {'Adams', 'Asotin', 'Benton', 'Chelan', 'Columbia', 'Douglas', 'Ferry', 'Franklin', 'Garfield', 'Grant', 'Kittitas', 'Klickitat', 'Lincoln', 'Okanogan', 'PendOreille', 'Spokane', 'Stevens', 'WallaWalla', 'Whitman', 'Yakima'};

% copy fasta to a new variable, loop over all the sequences and replace the WA county with the region in the 4th group split by -
fasta_new = fasta(keep);
locs = cell(0,0);
for i = 1:length(fasta_new)
    header = strsplit(fasta_new(i).Header, '|');
    % replace County and WA in header{4} with nothing
    header{4} = strrep(header{4}, 'WA','');
    header{4} = strrep(header{4}, 'County','');
    if ismember(header{4}, west)
        header{4} = 'WesternWashington';
    elseif ismember(header{4}, east)
        header{4} = 'EasternWashington';
    end
    fasta_new(i).Header = strjoin(header, '|');
    locs{end+1} = header{4};
end
disp(unique(locs))
% save fasta_new to new file SARSCoV2_WA.fasta overwriting existing file
delete('SARSCoV2_WA.fasta');
fastawrite('SARSCoV2_WA.fasta', fasta_new);


%% build the beast1 and 2 xmls using mascot_template.xml and dta_template.xml and the sequences in fasta_new
% first build the beast2, mascot xml.
% read in the template
f = fopen('mascot_template.xml');
g = fopen('xmls/SARS2_mascot_rep0.xml', 'w');
% loop over all lines in the template
while ~feof(f)
    line = fgets(f);
    if contains(line, 'insert_sequences')
        % if the line contains insert_sequences, loop over all sequences in fasta_new and write them to the file
        for i = 1:length(fasta_new)
            fprintf(g, ['<sequence id="seq_' fasta_new(i).Header '" spec="Sequence" taxon="' fasta_new(i).Header '" totalcount="4" value="' fasta_new(i).Sequence(250:end-250) '"/>\n']);
        end
    elseif contains(line, 'insert_times')
        tmp = '';
        for i = 1:length(fasta_new)
            time = strsplit(fasta_new(i).Header, '|');
            tmp = [tmp ',' fasta_new(i).Header '=' time{3}];
        end
        fprintf(g, strrep(line, 'insert_times', tmp(2:end)));
    elseif contains(line, 'insert_states')
        tmp = '';
        for i = 1:length(fasta_new)
            state = strsplit(fasta_new(i).Header, '|');
            tmp = [tmp ',' fasta_new(i).Header '=' state{4}];
        end
        fprintf(g, strrep(line, 'insert_states', tmp(2:end)));
    else
        % otherwise just write the line to the file
        fprintf(g, line);
    end
end
fclose(f);
fclose(g);
copyfile('xmls/SARS2_mascot_rep0.xml', 'xmls/SARS2_mascot_rep1.xml');
copyfile('xmls/SARS2_mascot_rep0.xml', 'xmls/SARS2_mascot_rep2.xml');



% Now build the beast1, dta xml, the same way as above
f = fopen('dta_template.xml');
g = fopen('xmls/SARS2_dta_rep0.xml', 'w');
while ~feof(f)
    line = fgets(f);
    if contains(line, 'insert_dat')
        % for each sequence, add sampling time and sampling state as
		% <taxon id="...">
		% 	<date value="..." direction="forwards" units="years"/>
		% 	<attr name="...">
		% 		human
		% 	</attr>
		% </taxon>
        for i = 1:length(fasta_new)
            header = strsplit(fasta_new(i).Header, '|');
            fprintf(g, ['\t\t<taxon id="' fasta_new(i).Header '">\n']);
            % for the date in format yyyy-mm-dd, make it dd/mm/yyyy
            date = strsplit(header{3}, '-');
            fprintf(g, ['\t\t\t<date value="' date{3} '/' date{2} '/' date{1} '" direction="forwards" units="years"/>\n']);
            fprintf(g, ['\t\t\t<attr name="division">\n']);
            fprintf(g, ['\t\t\t\t' header{4} '\n']);
            fprintf(g, ['\t\t\t</attr>\n']);
            fprintf(g, ['\t\t</taxon>\n']);
        end
    elseif contains(line, 'insert_seqs')
        % for each sequennce in the form of 
        %     <sequence>
        %     <taxon idref="..."/>
        %     seq
        % </sequence>
        for i = 1:length(fasta_new)
            fprintf(g, ['<sequence>\n']);
            fprintf(g, ['<taxon idref="' fasta_new(i).Header '"/>\n']);
            fprintf(g, [fasta_new(i).Sequence(250:end-250) '\n']);
            fprintf(g, ['</sequence>\n']);
        end
    elseif contains(line, 'insert_filename')
        fprintf(g, strrep(line, 'insert_filename', 'SARS2_dta_rep0'));
    else
        fprintf(g, line);
    end
end
fclose(f);
fclose(g);
% make 2 copies of the dta xml replacing the filename with rep1 and rep2
for i = 1: 2
    f = fopen('xmls/SARS2_dta_rep0.xml');
    g = fopen(['xmls/SARS2_dta_rep' num2str(i) '.xml'], 'w');    
    while ~feof(f)
        line = fgets(f);
        if contains(line, 'SARS2_dta_rep0')
            fprintf(g, strrep(line, 'SARS2_dta_rep0', ['SARS2_dta_rep' num2str(i)]));
        else
            fprintf(g, line);
        end
    end
end



