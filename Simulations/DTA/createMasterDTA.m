%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% creates master files from the template with the rates being linear
% combinations of covariates plus and error term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

% set the random number generator
rng(1);

% open the template
f = fopen('Rates_master.xml','r');   

% cell to save the lines
temp_lines = cell(0,0);  

system('rm -r master');
system('mkdir master');


% while there are still lines
while ~feof(f) 
    % read line by line
    temp_lines{end+1,1} = fgets(f);   
end

% close the template xml    
fclose(f); 

states = 3;
interval_length = 1;
time_intervals = 0:interval_length:5;

nr_reps = 100;
S = 1;

g = fopen('rates.txt','w');
fprintf(g, 'run');
for a = 1 : states
    for b = 1 : states
        if a~=b
            fprintf(g, '\tf_migrationRatesSkyline.state%d_to_state%d', a-1, b-1);
        end
    end
end
fprintf(g, '\n');


% compute birth values for R0's and k

unin_rate = 25;
x = 1;%[1:100];



scalers = zeros(0,0);
while S <= nr_reps
    filename = sprintf('DTA_S%d_master',S);
        
    start = normrnd(0,1,1,states);
    disp(start)
        
    fprintf(g, '%d', S);
    
    R0 = lognrnd(1,0.2,states,1);

        
    % sample the Ne and the migration scalers
    if S>50
        m_tot_scaler = 5;
    else
        m_tot_scaler = 1;
    end
    Ne_tot_scaler = 1;
    migration_rates = exprnd(m_tot_scaler, states,states);
%     migration_rates = migration_rates - diag(diag(migration_rates));
% 
%     while mean(migration_rates)<2.5
%         migration_rates = exprnd(m_tot_scaler, states,states);
%         migration_rates = migration_rates - diag(diag(migration_rates));
%     end
    
    for a = 1 : states
        for b = 1 : states
            if a~=b
                fprintf(g, '\t%f', migration_rates(a,b));
            end
        end
    end
    fprintf(g, '\n');

%     
%     migration = ones(states,states,1);    
%     migration = migration * m_tot_scaler;
       
    % check if the maximal rates are below a 5
%     subplot(nr_reps/5,5,S)
%     for a = 1 : states
%         plot(time_intervals,log(Ne(a,1:end))); hold on
%     end

%     plot(time_intervals, zeros(size(time_intervals)))
%     legend('1','2')


    fname = sprintf('./master/%s.xml',filename);   % set the file name
    p = fopen(fname,'w');

    Ne_nr = 1;
    m_nr = 1;
    
    if S<11
        popSize = 2000 * ones(states, 1);
        sampling = unin_rate*0.05*exprnd(1,states,1);
    end
    
    fprintf('%d %s\n', S, sprintf('%d ',popSize));

    for l = 1 : length(temp_lines)
        if ~isempty(strfind(temp_lines{l},'insert_infection'))
            for a = 1 : states
                offset = rand*2*pi;
                period = rand*0.5;
                
                for i = 1 : length(x)
                    rate_string = sprintf('%f', R0(a)*unin_rate/popSize(a));

                    for time=0.1:0.1:10
                        rate_string = [rate_string, sprintf(',%f:%f', R0(a)*unin_rate/popSize(a), time)];
                    end
                    
                    
                    fprintf(p, '\t\t\t<reaction spec=''Reaction'' rate="%s">\n',rate_string);
                    fprintf(p, '\t\t\t\tI[%d]:1 +%dS[%d]:1 -> %dI[%d]:1\n',a-1, x(i),a-1,x(i)+1,a-1);
                    fprintf(p, '\t\t\t</reaction>\n');
%                     fprintf(g, '\t\t</reactionGroup>\n');
                end
            end
        elseif ~isempty(strfind(temp_lines{l},'insert_recovery'))
            
            for i = 1 : states
                fprintf(p, '\t\t\t\t<reaction spec=''Reaction'' rate="%f">\n',unin_rate-sampling(i));
                fprintf(p, '\t\t\t\t\tI[%d]:1 -> R[%d]\n',i-1,i-1);
                fprintf(p, '\t\t\t\t</reaction>\n');
            end

        elseif ~isempty(strfind(temp_lines{l},'insert_sampling'))            
            for i = 1 : states
                offset = rand*2*pi;
                period = rand;

                
%                 rate_string = sprintf('%f', (sin(offset)+2)/3*0.1*unin_rate);
% 
%                 for time=0.1:0.1:10
%                     rate_string = [rate_string, sprintf(',%f:%f', (sin(time/period*pi+offset)+2)/3*0.01*unin_rate, time)];
%                 end

                if S<11
                    rate_string = num2str(0.05*unin_rate*10000/popSize(i));
                else
                    rate_string = num2str(0.5*unin_rate*10000/popSize(i));
                end
                rate_string = num2str(sampling(i));
                
                fprintf(p, '\t\t\t\t<reaction spec=''Reaction'' rate="%s">\n',rate_string);
                fprintf(p, '\t\t\t\t\tI[%d]:1 -> sample[%d]\n',i-1,i-1);
                fprintf(p, '\t\t\t\t</reaction>\n');
            end
        elseif ~isempty(strfind(temp_lines{l},'insert_waning'))
            delay = [0,0,0,0,0];
            for i = 1 : states
                fprintf(p, '\t\t\t\t<reaction spec=''Reaction'' rate="%f">\n',0);
                fprintf(p, '\t\t\t\t\tR[%d] -> S[%d]\n',i-1,i-1);
                fprintf(p, '\t\t\t\t</reaction>\n');
            end


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
            end
        elseif ~isempty(strfind(temp_lines{l},'insert_dimension'))
            fprintf(p,'%s',strrep(temp_lines{l},'insert_dimension',num2str(states)));
        elseif ~isempty(strfind(temp_lines{l},'insert_filename'))
            fprintf(p,'%s',strrep(temp_lines{l},'insert_filename',filename));
        else
            fprintf(p,'%s',temp_lines{l});  % print line unchanged
        end
    end
    fclose(p); %close file again
    
    cd('master');
    [~,~] = system(sprintf('/Applications/BEAST\\ 2.6.4/bin/beast %s', sprintf('SIR_S%d_master.xml',S)));
    cd('..');
    
    a=dir(sprintf('master/SIR_S%d_master.tree',S));
    
    if a.bytes>100  
        S = S + 1;
    end

end
fclose('all')

%%
% fclose('all');
% cd('master');
% files = dir('*.xml');
% for i = 1 : length(files)
%     system(sprintf('/Applications/BEAST\\ 2.6.4/bin/beast %s', files(i).name))
% end
% cd('..');

