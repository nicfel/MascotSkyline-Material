clear
% set the random number generator
rng(1);
khjljjkllkj

% open the template
f = fopen('Rates_master.xml','r');   

% cell to save the lines
temp_lines = cell(0,0);  

% prevents accidental running
dasjkdlasaslkd
system('rm -r master');
system('mkdir master');


% while there are still lines
while ~feof(f) 
    % read line by line
    temp_lines{end+1,1} = fgets(f);   
end

% close the template xml    
fclose(f); 

states = 2;
interval_length = 1;
time_intervals = 0:interval_length:5;

S = 1;

migrates = [repmat(5,200,1);repmat(25,200,1);repmat(5,100,1);repmat(25,100,1);
    repmat(5,100,1);repmat(25,100,1);repmat(5,100,1);repmat(25,100,1)];
samples = [repmat([repmat(250,100,1);repmat(500,100,1)],2,1);...
                repmat(250,200,1);repmat(50,200,1);repmat(4000,200,1)];

r0_random = [repmat(false,400,1);repmat(true,200,1);repmat(false,400,1)];
recovery_rate = [repmat(0,600,1);repmat(1,200,1);repmat(0,200,1)];



g = fopen('rates.txt','w');
fprintf(g, 'run');
fprintf(g, '\tavg_migration');
fprintf(g, '\tsamples');
fprintf(g, '\trandom_r0');

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
while S <= length(r0_random)
    filename = sprintf('SIR_S%d_master',S);
        
    start = normrnd(0,1,1,states);
    disp(start)
        
    fprintf(g, '%d', S);
    fprintf(g, '\t%f', migrates(S));
    fprintf(g, '\t%f', samples(S));
    fprintf(g, '\t%d', r0_random(S));


    % samples R0
    if r0_random(S)
        R0 = lognrnd(0.2,1,100000,1);
    else
        R0 = ones(states, 1)*1.5;
    end
    m_tot_scaler = migrates(S);

    migration_rates = exprnd(m_tot_scaler, states,states);
    
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
    
    
    popSize = floor(rand(states,1)*99500+ 500);
    if recovery_rate(S)==1
        sampling_proportion = repmat(betarnd(1,150, 1,1), 1, states);
    else
        sampling_proportion = ones(1, states);
    end

    for a = 1 : states
        fprintf(g, '\t%f', popSize(a));
    end

    for a = 1 : states
        fprintf(g, '\t%f', sampling_proportion(a));
    end
    fprintf(g, '\n');


    
    fprintf('%d %s\n', S, sprintf('%d ',popSize));
    skip_lines = 0;
    for l = 1 : length(temp_lines)
        if skip_lines>0
            skip_lines = skip_lines-1;
        elseif ~isempty(strfind(temp_lines{l},'insert_infection'))
            for a = 1 : states
                offset=rand;
                rate_string = sprintf('%f', (sin(0*2*pi+offset)+4)/5*R0(a)*unin_rate/popSize(a));
                % rate_string = sprintf('%f', R0(a)*unin_rate/popSize(a));
                for time=0.1:0.1:10
                    rate_string = [rate_string, sprintf(',%f:%f', (sin(time*2*pi+offset)+4)/5*R0(a)*unin_rate/(popSize(a)), time)];
                    % rate_string = [rate_string, sprintf(',%f:%f', R0(a)*unin_rate/(popSize(a)), time)];
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
            fprintf(p,'%s',strrep(temp_lines{l},'insert_start', num2str(randi(states)-1)));
        elseif  ~isempty(strfind(temp_lines{l},'insert_migration'))
            counter = 1;
            for a = 1 : states
                for b = 1 : states
                    if a ~= b
                        fprintf(p, '\t\t\t\t<reaction spec=''Reaction'' rate="%f">\n',migration_rates(a,b));
                        fprintf(p, '\t\t\t\t\tI[%d]:1 -> I[%d]:1\n',a-1,b-1);
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
        elseif contains(temp_lines{l}, '<inheritancePostProcessor spec=''LineageSampler''>') && recovery_rate(S)==1
            skip_lines = 7;
        elseif contains(temp_lines{l}, 'insert_samples')
            fprintf(p,'%s',strrep(temp_lines{l},'insert_samples', num2str(samples(S))));
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
    [~,~] = system(sprintf('/Applications/BEAST\\ 2.7.4/bin/beast %s', sprintf('SIR_S%d_master.xml',S)));
    % system(sprintf('/Applications/BEAST\\ 2.7.4/bin/beast %s', sprintf('SIR_S%d_master.xml',S)))
    cd('..');
    
    if S<801
        a=dir(sprintf('master/SIR_S%d_master.tree',S));
    else
        a=dir(sprintf('master/SIR_S%d_master.original.tree',S));
    end
    
    if a.bytes>100  && a.bytes < 7000000
        S = S + 1;
    end

end
fclose('all')

