clear

% set the random number generator
rng(10);

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
interval_length = 0.01;
time_intervals = 0:interval_length:10;

S = 1;

g = fopen('rates.txt','w');
fprintf(g, 'run');
for i = 1 : states
    fprintf(g, '\tCapacity.state%d', i-1);
end
for i = 1 : states
    fprintf(g, '\tCarryingProportion.state%d', i-1);
end
for i = 1 : states
    fprintf(g, '\tGrowthRate.state%d', i-1);
end
for a = 1 : states
    for b = 1 : states
        if a~=b
            fprintf(g, '\tf_migrationRatesSkyline.state%d_to_state%d', a-1, b-1);
        end
    end
end
fprintf(g, '\n');


scalers = zeros(0,0);
while S <= 100
    filename = sprintf('Logistic_S%d_master',S);
        
    Ne = zeros(states,length(time_intervals));
    growth = normrnd(-0.5,0.25,1,states);    
    carrying = normrnd(0,1,1,states);
    percentage_carry = rand(1,states)*0.75+0.25;
   
%    capacity.getValue()/(1 + (1-cP.getValue())/cP.getValue() * Math.exp(-t*growthRate.getArrayValue()))
    for a = 1 : states
        for i = 1 : length(time_intervals)
            Ne(a,i) = exp(carrying(a))/(1+ (1-percentage_carry(a))/ percentage_carry(a)*exp(-growth(a)*time_intervals(i)));
        end
    end
        
    % sample the Ne and the migration scalers
    m_tot_scaler = 0.5;
    
    migration_rates = exprnd(0.5, states,states);
    

    
    migration = ones(states,states,length(time_intervals));
    
    counter = 1;
    for a = 1 : states
        for b = 1 : states
            if a ~= b
                m_in = 'abcdef';
                for j = 1 : length(time_intervals)
                    migration(b,a,j) =...
                        Ne(a,j)*exp(migration_rates(a,b))/Ne(b,j);
                    
                end
                counter = counter + 1;
           end
        end
    end
    
    
    % check if the maximal rates are below a 5
    if max(migration(:)) < 100000 && max(1./(Ne(:))) < 100000
        
        fprintf(g, '%d', S);    
        for i = 1 : states
            fprintf(g, '\t%f', carrying(i));
        end
        for i = 1 : states
            fprintf(g, '\t%f', percentage_carry(i));
        end

        for i = 1 : states
            fprintf(g, '\t%f', growth(i));
        end

        
        for a = 1 : states
            for b = 1 : states
                if a~=b
                    fprintf(g, '\t%f', migration_rates(a,b));
                end
            end
        end
        fprintf(g, '\n');

            
        for a = 1 : states
            plot(time_intervals,Ne(a,1:end)); hold on
        end
        legend('1','2','3')

        S = S + 1;
        
        fname = sprintf('./master/%s.xml',filename);   % set the file name
        p = fopen(fname,'w');



        % print all the glm data to file

        ri = randi(states,40,1);
        for i = 1 : states
            sample_nr(i) =  sum(ri==i);
        end

        while min(sample_nr)==0
            ri = randi(states,50,1);
            for i = 1 : states
                sample_nr(i) =  sum(ri==i);
            end
        end


        sample_nr = sample_nr*20;

        Ne_nr = 1;
        m_nr = 1;

        for l = 1 : length(temp_lines)
            if ~isempty(strfind(temp_lines{l},'insert_coalescent'))
                for i = 1 : states
                    ne_in = 'abcdef';
                    for j = 1 : length(time_intervals)
                        ne_in = [ne_in ',' num2str(1/(2*Ne(i,j))) ':' num2str(time_intervals(j))];
                    end
                    ne_in = strrep(ne_in, 'abcdef,','');
                    fprintf(p, '\t\t\t\t<reaction spec=''Reaction'' rate="%s">\n',ne_in);
                    fprintf(p, '\t\t\t\t\t2L[%d]:1 -> L[%d]:1\n',i-1,i-1);
                    fprintf(p, '\t\t\t\t</reaction>\n');
                end
            elseif  ~isempty(strfind(temp_lines{l},'insert_migration'))
                counter = 1;
                for a = 1 : states
                    for b = 1 : states
                        if a ~= b
                            m_in = 'abcdef';
                            for j = 1 : length(time_intervals)
                                m_in = [m_in ',' num2str(...
                                    migration(a,b,j)...
                                    ) ':' num2str(time_intervals(j))];
                            end
                            m_in = strrep(m_in, 'abcdef,','');
                            fprintf(p, '\t\t\t\t<reaction spec=''Reaction'' rate="%s">\n',m_in);
                            fprintf(p, '\t\t\t\t\tL[%d]:1 -> L[%d]:1\n',a-1,b-1);
                            fprintf(p, '\t\t\t\t</reaction>\n');
                            counter = counter + 1;
                        end
                    end
                end
            elseif  ~isempty(strfind(temp_lines{l},'insert_samples'))
                fprintf(p,'\t\t\t<lineageSeedMultiple spec="MultipleIndividuals" copies="1" time="0">\n');
                fprintf(p,'\t\t\t\t<population spec="Population" type="@L" location="0"/>\n');
                fprintf(p,'\t\t\t</lineageSeedMultiple>\n');

                for i = 1 : states
                    rest_samples = sample_nr(i);
                    next_samples = 1;
                    while rest_samples > 0
                        time = 5*rand;
                        fprintf(p,'\t\t\t<lineageSeedMultiple spec="MultipleIndividuals" copies="%d" time="%.4f">\n',next_samples, time);
                        fprintf(p,'\t\t\t\t<population spec="Population" type="@L" location="%d"/>\n',i-1);
                        fprintf(p,'\t\t\t</lineageSeedMultiple>\n');
                        rest_samples = rest_samples - next_samples;
                        next_samples = 1;
                    end
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
    end
    
end
fclose('all')

% fclose('all');
% cd('master');
% files = dir('*.xml');
% for i = 1 : length(files)
%     system(sprintf('/Applications/BEAST\\ 2.6.4/bin/beast %s', files(i).name))
% end
% cd('..');

