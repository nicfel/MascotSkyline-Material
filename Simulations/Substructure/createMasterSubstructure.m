clear
% set the random number generator
rng(1);

% open the template
f = fopen('Rates_master.xml','r');   

% cell to save the lines
temp_lines = cell(0,0);  

% prevents accidental running
system('rm -r master');
system('mkdir master');


% while there are still lines
while ~feof(f) 
    % read line by line
    temp_lines{end+1,1} = fgets(f);   
end

% close the template xml    
fclose(f); 
states = 50;
interval_length = 1;
time_intervals = 0:interval_length:5;
S = 1;




g = fopen('rates.txt','w');
fprintf(g, 'run');

for a = 1 : states
    fprintf(g, '\tR0.state%d', a-1);
end

for a = 1 : states
    for b = 1 : states
        if a~=b
            fprintf(g, '\tf_migrationRatesSkyline.state%d_to_state%d', a-1, b-1);
        end
    end
end
fprintf(g, '\tpopSize0');
fprintf(g, '\tpopSize1');
fprintf(g, '\tsamplingProportion0');
fprintf(g, '\tsamplingProportion1');

fprintf(g, '\n');


% compute birth values for R0's and k
unin_rate = 52;

scalers = zeros(0,0);
while S <= 200
    filename = sprintf('Substructure_S%d_master',S);        
    start = normrnd(0,1,1,states);       
    fprintf(g, '%d', S);

    R0 = ones(states, 1)*10;
    R0(1) = 1.5;
    R0(2:20) = 5;

    m_tot_scaler = 0.05;

    migration_rates = zeros(states,states);
    migration_rates(1,:) = exprnd(m_tot_scaler, states,1);

    for a = 1 : states
        fprintf(g, '\t%f', R0(a));
    end
    
    
    for a = 1 : states
        for b = 1 : states
            if a~=b
                fprintf(g, '\t%f', migration_rates(a,b));
            end
        end
    end

    fname = sprintf('./master/%s.xml',filename);   % set the file name
    p = fopen(fname,'w');
    
    
    popSize = ones(1,states)*5;
    popSize(2:20) = 20;
    popSize(1) = 10000;
    sampling_proportion = ones(1, states)*min(S, 100)/100;
    correction = min(abs(S-200),100)/100;
    sampling_proportion(2:20) = 0.25*correction;
    sampling_proportion(1) = 0.025*correction;

    for a = 1 : states
        fprintf(g, '\t%f', popSize(a));
    end

    for a = 1 : states
        fprintf(g, '\t%f', sampling_proportion(a));
    end
    fprintf(g, '\n');
    
    fprintf('%d %s\n', S, sprintf('%d ',popSize));
    skip_lines = 0;
    offset=rand(1,states);

    for l = 1 : length(temp_lines)
        if skip_lines>0
            skip_lines = skip_lines-1;
        elseif ~isempty(strfind(temp_lines{l},'insert_infection'))
            for a = 1 : states
                if a==1
                    rate_string = sprintf('%f', (sin(0*2*pi+offset(a))+4)/5*R0(a)*unin_rate/popSize(a));
                    for time=0.1:0.1:10
                        rate_string = [rate_string, sprintf(',%f:%f', (sin(time*2*pi+offset(a))+4)/5*R0(a)*unin_rate/(popSize(a)), time)];
                    end    
                else
                    rate_string = sprintf('%f', (sin(0*2*pi+offset(a))+1)/2*R0(a)*unin_rate/popSize(a));
                    % rate_string = sprintf('%f', R0(a)*unin_rate/popSize(a));
                    for time=0.1:0.1:10
                        rate_string = [rate_string, sprintf(',%f:%f', (sin(time*2*pi+offset(a))+1)/2*R0(a)*unin_rate/(popSize(a)), time)];
                        % rate_string = [rate_string, sprintf(',%f:%f', R0(a)*unin_rate/(popSize(a)), time)];
                    end    
                end
                
                fprintf(p, '\t\t\t<reaction spec=''Reaction'' rate="%s">\n',rate_string);
                fprintf(p, '\t\t\t\tI[%d]:1 +%dS[%d]:1 -> %dI[%d]:1\n',a-1, 1,a-1,2,a-1);
                fprintf(p, '\t\t\t</reaction>\n');
            end
        elseif ~isempty(strfind(temp_lines{l},'insert_recovery'))
            
            for i = 1 : states
                fprintf(p, '\t\t\t\t<reaction spec=''Reaction'' rate="%f">\n',unin_rate-unin_rate*sampling_proportion(i));
                fprintf(p, '\t\t\t\t\tI[%d]:1 -> R[%d]\n',i-1,i-1);
                fprintf(p, '\t\t\t\t</reaction>\n');
            end

        elseif ~isempty(strfind(temp_lines{l},'insert_sampling'))            
            for i = 1 : states                
                fprintf(p, '\t\t\t\t<reaction spec=''Reaction'' rate="%f">\n',unin_rate*sampling_proportion(i));
                fprintf(p, '\t\t\t\t\tI[%d]:1 -> sample[%d]\n',i-1,i-1);
                fprintf(p, '\t\t\t\t</reaction>\n');
            end
        elseif ~isempty(strfind(temp_lines{l},'insert_waning'))
            delay = [0,0,0,0,0];
            for i = 1 : states
                fprintf(p, '\t\t\t\t<reaction spec=''Reaction'' rate="%f">\n',0*unin_rate);
                fprintf(p, '\t\t\t\t\tR[%d] -> S[%d]\n',i-1,i-1);
                fprintf(p, '\t\t\t\t</reaction>\n');
            end

        elseif ~isempty(strfind(temp_lines{l},'insert_start'))
            fprintf(p,'%s',strrep(temp_lines{l},'insert_start', num2str(0)));
        elseif  ~isempty(strfind(temp_lines{l},'insert_migration'))
            counter = 1;
            for a = 1 : states
                for b = 1 : states
                    if a ~= b
                        rate_string = sprintf('%f', (sin(0*2*pi+offset(b))+1)/2*migration_rates(a,b)/popSize(b));
                        for time=0.1:0.1:10
                            rate_string = [rate_string, sprintf(',%f:%f', (sin(time*2*pi+offset(b))+1)/2*migration_rates(a,b)/popSize(b), time)];
                        end    

                        fprintf(p, '\t\t\t\t<reaction spec=''Reaction'' rate="%s">\n',rate_string);
                        fprintf(p, '\t\t\t\t\tI[%d]:1 +S[%d]:1 -> I[%d]:1 + I[%d]:1\n',a-1,b-1,a-1,b-1);
                        fprintf(p, '\t\t\t\t</reaction>\n');
                        counter = counter + 1;
                    end
                end
            end
        elseif contains(temp_lines{l}, 'insert_population_size')
            for a = 1 : states  
                fprintf(p, '\t\t\t<populationSize spec=''PopulationSize'' size=''%d''>\n', popSize(a));
                fprintf(p, '\t\t\t\t<population spec=''Population'' type=''@S'' location="%d"/>\n', a-1);
                fprintf(p, '\t\t\t</populationSize>\n');
                % fprintf(p, '\t\t\t<populationSize spec=''PopulationSize'' size=''%d''>\n', floor(popSize(a)/5));
                % fprintf(p, '\t\t\t\t<population spec=''Population'' type=''@I'' location="%d"/>\n', a-1);
                % fprintf(p, '\t\t\t</populationSize>\n');
                % fprintf(p, '\t\t\t<populationSize spec=''PopulationSize'' size=''%d''>\n', popSize(a));
                % fprintf(p, '\t\t\t\t<population spec=''Population'' type=''@R'' location="%d"/>\n', a-1);
                % fprintf(p, '\t\t\t</populationSize>\n');
            end
        elseif contains(temp_lines{l}, '<inheritancePostProcessor spec=''LineageSampler''>')
            skip_lines = 7;
        elseif contains(temp_lines{l}, 'insert_samples')
            fprintf(p,'%s',strrep(temp_lines{l},'insert_samples', num2str(0)));
        elseif ~isempty(strfind(temp_lines{l},'insert_dimension'))
            fprintf(p,'%s',strrep(temp_lines{l},'insert_dimension',num2str(states)));
        elseif contains(temp_lines{l}, 'insert_filename.tree') && S>800
            fprintf(p,'%s',strrep(temp_lines{l},'insert_filename',[filename '.original']));
        elseif ~isempty(strfind(temp_lines{l},'insert_filename'))
            fprintf(p,'%s',strrep(temp_lines{l},'insert_filename',filename));
        else
            fprintf(p,'%s',temp_lines{l});  % print line unchanged
        end
    end
    fclose(p); %close file again
    
    cd('master');
    [~,~] = system(sprintf('/Applications/BEAST\\ 2.7.6/bin/beast %s', sprintf('Substructure_S%d_master.xml',S)));
    % system(sprintf('/Applications/BEAST\\ 2.7.4/bin/beast %s', sprintf('SIR_S%d_master.xml',S)))
    cd('..');
    
    a=dir(sprintf('master/Substructure_S%d_master.tree',S));
    if a.bytes > 10000 && a.bytes < 7000000
        % Open the tree file for reading and writing
        tree_filename = sprintf('master/Substructure_S%d_master.tree',S);
        tree_file = fopen(tree_filename, 'r');        
        % Read the entire tree file into a cell array, line by line
        tree_lines = {};
        while ~feof(tree_file)
            tree_lines{end+1} = fgets(tree_file);
        end
        fclose(tree_file);
        
        % Open the tree file for writing (this will overwrite the original)
        tree_file = fopen(tree_filename, 'w');
        
        % Loop through each line and replace location values
        for l = 1:length(tree_lines)
            line = tree_lines{l};
            
            % Replace location="%d" for d between 2 and 30 with location="2"
            for d = 1:19
                line = strrep(line, sprintf('location="%d"', d), 'location="0"');
            end
            
            % Replace location="%d" for d between 31 and 40 with location="3"
            for d = 20:states-1
                line = strrep(line, sprintf('location="%d"', d), 'location="1"');
            end
            
            % Write the modified line back to the file
            fprintf(tree_file, '%s', line);
        end
        fclose(tree_file);

        tree_file = fopen(strrep(tree_filename,'.tree','.background.tree'), 'w');
        % Loop through each line and replace location values
        for l = 1:length(tree_lines)
            line = tree_lines{l};
            
            % Replace location="%d" for d between 2 and 30 with location="2"
            for d = 1:19
                line = strrep(line, sprintf('location="%d"', d), 'location="1"');
            end
            
            % Replace location="%d" for d between 31 and 40 with location="3"
            for d = 20:states-1
                line = strrep(line, sprintf('location="%d"', d), 'location="2"');
            end
            
            % Write the modified line back to the file
            fprintf(tree_file, '%s', line);
        end
        fclose(tree_file);


        S = S + 1;
    end
    

end
fclose('all')

