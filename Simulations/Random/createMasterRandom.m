%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% creates master files from the template with the rates being linear
% combinations of covariates plus and error term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
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
time_intervals = 0:1:50;

% set the random number generator
rng(1)

S = 1;


scalers = zeros(0,0);
while S <= 1
    filename = sprintf('Rates_S%d_master',S);
        
    Ne_log = zeros(states*length(time_intervals),1);
    
    val = [-4,-2,3]
    
    Ne_log(1) = 3;
    Ne_log(2) = 2;
    Ne_log(3) = -3;
    
    for a = 1 : states
        for i = 2 : length(time_intervals)
            Ne_log(states*(i-1) + a) = Ne_log(states*(i-2) + a)+val(a)*0.05;
        end
%         Ne_log(a:states:end) = Ne_log(a:states:end)-mean(Ne_log(a:states:end));
    end
    
    
    % sample the Ne and the migration scalers
    m_tot_scaler = 0.5;
    Ne_tot_scaler = 1;
    
    migration = ones(states*(states-1)*length(time_intervals),1);
    
    migration = migration * m_tot_scaler;
    Ne = exp(Ne_log) * Ne_tot_scaler;
    for a = 1 : states
        plot(Ne(a:states:end)); hold on
    end
    legend('1','2','3')
    
    counter = 1;
    for a = 1 : states
        for b = 1 : states
            if a ~= b
                m_in = 'abcdef';
                for j = 1 : length(time_intervals)
                    migration((states*(states-1))*(j-1) + counter) =...
                        Ne(states*(j-1) + b)*migration((states*(states-1))*(j-1) + counter)/Ne(states*(j-1) + a);
                    
                end
                counter = counter + 1;
           end
        end
    end
    
    
    % check if the maximal rates are below a 5
    if max(migration) < 100000 && max(1./(Ne)) < 100000
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
                        ne_in = [ne_in ',' num2str(1/(2*Ne(states*(j-1) + i))) ':' num2str(time_intervals(j))];
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
                                    migration((states*(states-1))*(j-1) + counter)...
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
                        time = 25*rand;
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
