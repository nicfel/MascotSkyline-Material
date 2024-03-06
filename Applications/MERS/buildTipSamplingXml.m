rng(5335249)
% read in the xml to get the sequences
f = fopen('template/MERS_skyline.xml');
ids = cell(0,0);
seq = cell(0,0);
% tipdatesampling all sequences for which the sampling date is not yyyy/mm/ddd
tipdatesampling = cell(0,0);
while ~feof(f)
    line = fgets(f);
    if contains(line, '<sequence')
        tmp = strsplit(line, '"');
        tmp2 = strsplit(tmp{6}, '|');
        tmp3 = strsplit(tmp2{end}, '-');
        if length(tmp3) ~= 3
            tipdatesampling{end+1,1} = tmp{6};
        end            
        ids{end+1,1} = tmp{6};
        seq{end+1, 1} = tmp{10};
    end
end
fclose(f);

% compute the pairwise distances based on the sequences in seq, 
% while ignoring the N an ? characters and only counting the number of
% differences where neither sequence has an N, ?  or -and then divide by 
% the total comparisons to get the fraction of differences
% dist = zeros(length(seq), length(seq));
% for i = 1:length(seq)
%     for j = i+1:length(seq)
%         indices = ~(seq{i} == 'N' | seq{i} == '?' | seq{i} == '-' | seq{j} == 'N' | seq{j} == '?' | seq{j} == '-');
%         dist(i,j) = sum(seq{i}(indices) ~= seq{j}(indices))/sum(indices);
%     end
% end


% make 10 replicates to sample the tip states
for repeat = 1:10
    f = fopen('template/MERS_skyline.xml');

    filename = ['xmls/MERS_mascot_it' num2str(repeat) '_rep0.xml'];

    g = fopen(filename, 'w');

    locations = {'|camel|', '|human|'};
    tipsToSample=[];
    for i = 1:length(locations)
        idx = find(contains(ids, locations{i}));
        idx = idx(randperm(length(idx), 5));
        tipsToSample = [tipsToSample; ids(idx)];
    end


    while ~feof(f)
        line = fgets(f);
        if contains(line, '</state>')
            % add the tipsTOsample as integer parameters with a random state betwee 0 and 3
            for i = 1 :length(tipsToSample)
                fprintf(g, ['\t\t\t<parameter id="' tipsToSample{i} '.sampledState" spec="parameter.IntegerParameter" lower="0" upper="1" name="stateNode">' num2str(randi(2)-1) '</parameter>\n']);
            end
        % replace MASCOT with MascotWithTipSampling, MappedMascot with MappedMascotWithTipSampling and MigrationCountLogger with MigrationCountLoggerWithTipSampling
        elseif contains(line, 'mascot.distribution.Mascot')
            fprintf(g, strrep(line, 'mascot.distribution.Mascot', 'mascot.distribution.MascotWithTipSampling'));
            for i = 1:length(tipsToSample)
                fprintf(g, ['\t\t\t\t\t<tipStates idref="' tipsToSample{i} '.sampledState"/>\n']);
            end
            line = fgets(f);
        elseif contains(line, 'mascot.logger.MigrationCountLogger')
            fprintf(g, strrep(line, 'MigrationCountLogger', 'MigrationCountLoggerWithTipSampling'));
            line = fgets(f);
            fprintf(g, strrep(line, 'MappedMascot', 'MappedMascotWithTipSampling'));
            line = fgets(f);
            for i = 1:length(tipsToSample)
                fprintf(g, ['\t\t\t\t\t<tipStates idref="' tipsToSample{i} '.sampledState"/>\n']);
            end

        elseif contains(line, '<logger id="typedTreelogger.t:MERS_CoV_274"')
            fgets(f);fgets(f);
            line = fgets(f);
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
            for i = 1:length(tipsToSample)
                fprintf(g, ['\t\t<operator id="' tipsToSample{i} 'sampledState.sampling" parameter="@' tipsToSample{i} '.sampledState" spec="beast.base.inference.operator.IntUniformOperator" weight="0.1"/>\n']);
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
            for i = 1:length(tipsToSample)
                fprintf(g, ['\t\t<log idref="' tipsToSample{i} '.sampledState"/>\n']);
            end

            line = fgets(f);

        end
        fprintf(g, line);
    end
    fclose('all');


    % make the DTA xml
    f = fopen('dta_template.xml');
    g = fopen(strrep(filename, 'mascot','dta'), 'w');
    while ~feof(f)
        line = fgets(f);
        if contains(line, 'insert_seqs')
            for i = 1 : length(ids)
                fprintf(g,'\t\t<sequence>\n');
                fprintf(g,'\t\t\t<taxon idref="%s"/>\n', ids{i});
                fprintf(g,'\t\t\t%s\n', seq{i});
                fprintf(g,'\t\t</sequence>\n');
            end

            
        elseif contains(line, 'insert_dat')
            for i = 1 : length(ids)
                tmp = strsplit(ids{i}, '|', 'CollapseDelimiters',false);
                
                tmp2 = strsplit(tmp{4}, '-');
                
                isuncertain=false;
                if length(tmp2)==2
                    isuncertain=true;
                    tmp{4} = [tmp{4} '-15'];                    
                end
                tmp2 = strsplit(tmp{4}, '-');               
                
                
                fprintf(g,'\t\t<taxon id="%s">\n', ids{i});
                if isuncertain
                    fprintf(g,'\t\t\t<date value="%s/%s/%s" direction="forwards" units="years" uncertainty="0.083333333"/>\n', tmp2{3}, tmp2{2}, tmp2{1});
                else
                    fprintf(g,'\t\t\t<date value="%s/%s/%s" direction="forwards" units="years"/>\n', tmp2{3}, tmp2{2}, tmp2{1});
                end
                fprintf(g,'\t\t\t<attr name="division">\n');

                % if id is tip to sample insert ? otherwise insert true state
                if ismember(ids{i}, tipsToSample)
                    fprintf(g,'\t\t\t\t?\n');
                else
                    fprintf(g,'\t\t\t\t%s\n',tmp{3});
                end
                fprintf(g,'\t\t\t</attr>\n');
                fprintf(g,'\t\t</taxon>\n');
            end
        
        elseif contains(line, 'insert_filename')
            tmp = strrep(filename, 'xmls/', '');
            fname = strrep(tmp, 'mascot','dta');
            fprintf(g, strrep(line, 'insert_filename', strrep(fname, '.xml','')));

        % block that adds tip sampling
        elseif contains(line, '</treeModel>')
            % for each tip, add a block of the form         
            % <leafHeight taxon="DQ164186_Cb_2002">
            % <parameter id="age(DQ164186_Cb_2002)"/>
            % </leafHeight>
    
            for i = 1: length(tipdatesampling)
                fprintf(g, ['\t\t\t<leafHeight taxon="' tipdatesampling{i} '">\n']);
                fprintf(g, ['\t\t\t\t<parameter id="age(' tipdatesampling{i} ')"/>\n']);
                fprintf(g, ['\t\t\t</leafHeight>\n']);
            end
            fprintf(g, line);
        elseif contains(line, '</operators>')
            % add operators of the form 		
            % <uniformOperator weight="0.1">
			% <parameter idref="age(DQ164186_Cb_2002)"/>
            % </uniformOperator>

            for i = 1 : length(tipdatesampling)
                fprintf(g, ['\t\t<uniformOperator weight="0.1">\n']);
                fprintf(g, ['\t\t\t<parameter idref="age(' tipdatesampling{i} ')"/>\n']);
                fprintf(g, ['\t\t</uniformOperator>\n']);
            end
            fprintf(g, line);
        elseif contains(line, '<rateStatistic idref="geo.meanRate"/>')
            % add logger of the form 			<parameter idref="age(DQ164186_Cb_2002)"/>
              fprintf(g, line);
            for i = 1 : length(tipdatesampling)
                fprintf(g, ['\t\t<parameter idref="age(' tipdatesampling{i} ')"/>\n']);
            end
            
        
        else
            fprintf(g, line);
        end
    end

end


