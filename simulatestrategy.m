function [temp_fit, successful,keep]=simulatestrategy(replicates,...
    new_strategy,mort,cancer_danger,targetsize,withextmort,...
    extmortthreshold,celldeath,fitness_thr1,fitness_thr2,nof_onco_steps)
    mat_times=zeros(1,replicates); % when each stategy matures
    death_times=zeros(1,replicates); % when each strategy dies
    death_cause=zeros(1,replicates); % why each strategy dies
    fitness_vals_temp=zeros(replicates,length(mort)-1); % all fitness values
    cell_dist=zeros(1,replicates*3); % how many of which type of cells there were at the end
    killed_dist=zeros(1,replicates*3); % why cells were killed throughout the simulation
    successful=true;keep=false;
    for i=1:replicates
        % simulate the life of the organism and calculate fitness at the
        % end
        [M_res,LS_res,DC_res,fitness,cell_types,killed_cells]=...
            onelife(new_strategy,mort,cancer_danger,targetsize,...
            withextmort,extmortthreshold,celldeath,nof_onco_steps);
        mat_times(i)=M_res;death_times(i)=LS_res;death_cause(i)=DC_res;
        ind=i*3-2;
        cell_dist(1,ind:ind+2)=cell_types;
        killed_dist(1,ind:ind+2)=killed_cells;
        fitness_vals_temp(i,:)=fitness;
        current_mean=sum(fitness_vals_temp)/i;
        % if the current mean is below thr1 after three trials, stop
        % simulating, but if it is above thr2 at that point, keep that
        % strategy in the population (explained further in ga.m)
        if current_mean <= fitness_thr1 && i > 3
            if current_mean >= fitness_thr2
                keep=true;
            end
            successful=false;mean_fitness=current_mean;
            temp_fit=horzcat(new_strategy,mean_fitness,fitness_vals_temp',...
                mat_times,death_times,death_cause,cell_dist,killed_dist);
            break;
        end
    end
    if successful ~= false 
        keep=true;
        mean_fitness=calculatemeanfitness(fitness_vals_temp,replicates);
        temp_fit=horzcat(new_strategy,mean_fitness,fitness_vals_temp',...
            mat_times,death_times,death_cause,cell_dist,killed_dist);
    end        
end