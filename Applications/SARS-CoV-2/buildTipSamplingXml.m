rng(6541589)
% read in the fasta file
fasta_new = fastaread('SARSCoV2_WA.fasta');

% loop over 10 replicates
for repeat = 1:10
    % pick 5 members for each group that has in their name
    % EasternWashington, WesternWashington, NorthAmerica or restOfWorld
    % and save their header name
    locations = {'EasternWashington', 'WesternWashington', 'NorthAmerica', 'restOfWorld'};
    tipsToSample=[];
    for i = 1:length(locations)
        idx = find(contains({fasta_new.Header}, locations{i}));
        idx = idx(randperm(length(idx), 5));
        tipsToSample = [tipsToSample {fasta_new(idx).Header}];
    end
    %% build the beast1 and 2 xmls using mascot_template.xml and dta_template.xml and the sequences in fasta_new
    % first build the beast2, mascot xml.
    % read in the template
    f = fopen('mascot_template.xml');
    g = fopen(['xmls/SARS2_mascot_it' num2str(repeat) '_rep0.xml'], 'w');
    % loop over all lines in the template
    while ~feof(f)
        line = fgets(f);
        if contains(line, 'insert_sequences')
            % if the line contains insert_sequences, loop over all sequences in fasta_new and write them to the file
            for i = 1:length(fasta_new)
                fprintf(g, ['<sequence id="seq_' fasta_new(i).Header '" spec="Sequence" taxon="' fasta_new(i).Header '" totalcount="4" value="' fasta_new(i).Sequence(250:end-250) '"/>\n']);
            end
        elseif contains(line, '</state>')
            % add the tipsTOsample as integer parameters with a random state betwee 0 and 3
            for i = 1 :length(tipsToSample)
                fprintf(g, ['\t\t\t<parameter id="' tipsToSample{i} '.sampledState" spec="parameter.IntegerParameter" lower="0" upper="3" name="stateNode">' num2str(randi(4)-1) '</parameter>\n']);
            end
            fprintf(g, line);
        elseif contains(line, '<operator id="MascotWilsonBalding.t:SARSCoV2_WA')
            fprintf(g, line);
            % add the tipsToSample as a list
            for i = 1:length(tipsToSample)
                fprintf(g, ['\t\t<operator id="' tipsToSample{i} 'sampledState.sampling" parameter="@' tipsToSample{i} '.sampledState" spec="beast.base.inference.operator.IntUniformOperator" weight="0.1"/>\n']);
            end
        elseif contains(line, '</mappedMascot>')
            for i = 1:length(tipsToSample)
                fprintf(g, ['\t\t\t\t\t<tipStates idref="' tipsToSample{i} '.sampledState"/>\n']);
            end
            fprintf(g, line);
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
            % replace MASCOT with MascotWithTipSampling, MappedMascot with MappedMascotWithTipSampling and MigrationCountLogger with MigrationCountLoggerWithTipSampling
        elseif contains(line, 'mascot.distribution.Mascot')
            fprintf(g, strrep(line, 'mascot.distribution.Mascot', 'mascot.distribution.MascotWithTipSampling'));
            for i = 1:length(tipsToSample)
                fprintf(g, ['\t\t\t\t\t<tipStates idref="' tipsToSample{i} '.sampledState"/>\n']);
            end
        elseif contains(line, 'mascot.distribution.MappedMascot')
            fprintf(g, strrep(line, 'MappedMascot', 'MappedMascotWithTipSampling'));
        elseif contains(line, 'mascot.logger.MigrationCountLogger')
            fprintf(g, strrep(line, 'MigrationCountLogger', 'MigrationCountLoggerWithTipSampling'));
        elseif contains(line, '<logger id="treelog.t:SARSCoV2_WA"')
            fgets(f);fgets(f);fgets(f);fgets(f);fgets(f);
        elseif contains(line, '<log idref="HyperPrior.hyperNormal-sigma-SkylineNe.EasternWashington.Prior.t:SARSCoV2_WA"/>')
            for i = 1:length(tipsToSample)
                fprintf(g, ['\t\t<log idref="' tipsToSample{i} '.sampledState"/>\n']);
            end
        else
            % otherwise just write the line to the file
            fprintf(g, line);
        end
    end
    fclose(f);
    fclose(g);
    copyfile(['xmls/SARS2_mascot_it' num2str(repeat) '_rep0.xml'], ['xmls/SARS2_mascot_it' num2str(repeat) '_rep1.xml']);
    copyfile(['xmls/SARS2_mascot_it' num2str(repeat) '_rep0.xml'], ['xmls/SARS2_mascot_it' num2str(repeat) '_rep2.xml']);

    % Now build the beast1, dta xml, the same way as above
    f = fopen('dta_template.xml');
    g = fopen(['xmls/SARS2_dta_it' num2str(repeat) '_rep0.xml'], 'w');
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
                if ismember(fasta_new(i).Header,tipsToSample)
                    fprintf(g, ['\t\t\t\t?\n']);
                else
                    fprintf(g, ['\t\t\t\t' header{4} '\n']);
                end
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
            fprintf(g, strrep(line, 'insert_filename', ['SARS2_dta_it' num2str(repeat) '_rep0']));
        else
            fprintf(g, line);
        end
    end
    fclose(f);
    fclose(g);
    % make 2 copies of the dta xml replacing the filename with rep1 and rep2
    for i = 1: 2
        f = fopen(['xmls/SARS2_dta_it' num2str(repeat) '_rep0.xml']);
        g = fopen(['xmls/SARS2_dta_it' num2str(repeat) '_rep' num2str(i) '.xml'], 'w');    
        while ~feof(f)
            line = fgets(f);
            if contains(line, ['SARS2_dta_it' num2str(repeat) '_rep0'])
                fprintf(g, strrep(line, ['SARS2_dta_it' num2str(repeat) '_rep0'], ['SARS2_dta_it' num2str(repeat) '_rep' num2str(i)]));
            else
                fprintf(g, line);
            end
        end
    end
end

