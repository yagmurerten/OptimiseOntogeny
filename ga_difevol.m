% this is the main genetic algorithm function 

function ga_difevol(extmort,targetsize,cancer_danger,cellmort,replicate_ind,c_thr1,c_thr2)
    
    celldeath=1; % whether we want cells to die or not    
    nof_onco_steps=4; % this is =4 in here because that is the dimension of the matrix
    % but this means K=3, i.e. here number of onco. steps are from 1 to 4
    % in the text they are from 0 to 3
    
    % the maximum apoptosis threshold: this goes from 0 to 3, like K in the
    % text
    alt=nof_onco_steps-1;
    
    % bounds for the number of differentiation steps
    min_dif_layer=3;max_dif_layer=50; 
    
    % we have this many traits evolving
    length_strategy=7;
    
    % traits can evolve within these bounds
    strategy_bounds=[0.0, 0.0, 1.0, 1.0;... % prob asym
                  0.0, 0.0, 1.0, 1.0;... % prob dif
                  5, 5, 500, 500;... % tel length
                  0,0,alt,alt;... % apop_thr
                  0.0, 0.0, 1.0, 1.0;... % apop_percent
                  min_dif_layer,min_dif_layer,max_dif_layer,max_dif_layer;... % nof layer
                  0.01,0.01,200,200]; % div_prop

    if ~celldeath           
        cellmort=0.; % without cell death
        cancer_danger=0;
    end
           
    withextmort=0;extmortthreshold=0.001; 
    % switching stochastic extrinsic mortality on/off
    % if off (=0) the organism lives until it reaches to t at which 
    % the probability of it being alive is extmortthreshold
    % at that point we stop simulating

    replicates=10; % how many times we try a strategy
    
	reshape(strategy_bounds,length_strategy,4);dlmwrite('bounds.txt',strategy_bounds);
      
    % maximum # rounds at (stop ga if it did not converge until
    % this)
    rounds=5000;

    % # of stategies in the population
    popsize=100;
    
    num_candidate_strategies=10;
    candidate_indices=round(linspace(1,popsize,num_candidate_strategies));
    % pick these individuals from the population as potential gametes   

    filename1=strcat('allval_',num2str(extmort),'_',num2str(targetsize),...
            '_',num2str(c_thr2),'_',num2str(replicate_ind),'.txt');
        % the new 'best' strategies are added here
    filename2=strcat('bestval_',num2str(extmort),'_',num2str(targetsize),...
            '_',num2str(c_thr2),'_',num2str(replicate_ind),'.txt');
    convergencethreshold1=c_thr2;stable1=false;   
    
    % new ga from scratch
    if c_thr1==0       
        index=0;
        % this has the first set of strategies that are tried when forming the
        % initial population via the function found in shotgun.m
        filename5=strcat('shotgun_',num2str(extmort),'_',num2str(targetsize),'_',num2str(replicate_ind),'.txt');
        lastbest=0; 
        % shotgun is a function that randomly samples popsize potential
        % strategies and simulates them to initiate a population for the
        population=shotgun(extmort,cellmort,cancer_danger,targetsize,popsize,withextmort,replicates,...
            extmortthreshold,0,'',celldeath,strategy_bounds,nof_onco_steps);
        tmp=sortrows(population,-(length_strategy+1));
        dlmwrite(filename5, population);    
        total_fit=sum(population(:,length_strategy+1));
        beststrategy=tmp(1,:);
        dlmwrite(filename2, horzcat(0,beststrategy),'-append');
    else
        % read population
        % this is to continue from an already ran genetic algorithm
        filetoread=strcat('population_',num2str(extmort),'_',num2str(targetsize),...
            '_',num2str(c_thr1),'_',num2str(replicate_ind),'.txt');
        population = importdata(filetoread);
        tmp=sortrows(population,-(length_strategy+1));
        total_fit=sum(population(:,length_strategy+1));
        beststrategy=tmp(1,:);
        % read the indices of the population
        indexfile=strcat('indices',num2str(c_thr1),'_',num2str(replicate_ind),'.txt');
        A=importdata(indexfile);
        % find the index of the best
        lastbest=A(1,1);index=A(1,2);
        dlmwrite(filename2, horzcat(lastbest,beststrategy),'-append');
        % read all the values and get the index
    end
    
    % do ga if actually there are fit strategies
    % because for some parameter combinations, there might not be any
    % strategy that works, and then there is no point trying to optimise it
    if total_fit > 0
        while index<rounds && stable1==false
            index=index+1;
            % sort the population of strategies according to their mean fitness 
            temp_str=sortrows(population,-(length_strategy+1));

            % pick the strategies that have zero fitness and shuffle them around so
            % that we don't always pick the same zero-fitness strategy as a parent
            t_ind=length(temp_str)-length(temp_str(temp_str(:,length_strategy+1)<=0))+1;
            temp_t=temp_str(t_ind:end,:);
            temp_t =  temp_t(randperm(end),:);
            temp_str(t_ind:end,:)=temp_t;

            % these are the potential parent strategies 
            gametes=temp_str(candidate_indices,1:length_strategy);  

            % we calculate the relative fitness of these strategies
            gametes_f=temp_str(candidate_indices,1:length_strategy+1);
            total_fit=sum(gametes_f(:,length_strategy+1));
            fitnesses=gametes_f(:,length_strategy+1)/total_fit;

            % the ones with higher relative fitness are more likely to be picked as
            % the first 'parent strategy' (or here res_strategy)
            ind=ddists(fitnesses,1);
            res_strategy=gametes(ind,:);

            % create a mutant of the first parent
            new_strategy=mutate_strategy(res_strategy,strategy_bounds); 

            % recombine the new strategy mutated strategy with a random potential
            % parent
            new_strategy=recombine_strategy(new_strategy,gametes,ind);

            % each strategy is checked for replicate many times and its mean
            % fitness is compared to the 'best' and the 'worst' of the population
            % if it is better than the 'best' after 10 trials, it becomes the new
            % 'best' strategy
            % if it is not better than the 'best', but still has a higher mean
            % fitness than the 'worst', it is added to the population of potential
            % parents
            fitness_thr1=temp_str(1,length_strategy+1);
            fitness_thr2=temp_str(popsize,length_strategy+1);

            [temp_fit, successful, keep]=simulatestrategy(replicates,...
                new_strategy,[cellmort,extmort],cancer_danger,targetsize,...
                withextmort,extmortthreshold,celldeath,fitness_thr1,fitness_thr2,nof_onco_steps);
            dlmwrite(filename1, horzcat(index,temp_fit),'-append');

            if keep
                population=sortrows(population,-(length_strategy+1));
                population(popsize,:)=temp_fit;
                if successful
            % is our strategy better than our old 'best' strategies?
                    newfitness=temp_fit(length_strategy+1);
                    oldfitness=population(:,length_strategy+1);
                    isBetter=size(oldfitness(oldfitness<newfitness),1);
                    if isBetter==popsize-1
                        lastbest=index;
                        dlmwrite(filename2, horzcat(index,temp_fit),'-append');
                    end
                end
            end

            % if we haven't found a new 'best' for this many rounds (but found at
            % least one since the genetic algorithm started) we stop the genetic
            % algorithm
            if index-lastbest>=convergencethreshold1 
                stable1=true;
            end
        end
    end
    % the population at the end of n rounds
    filename6=strcat('population_',num2str(extmort),'_',num2str(targetsize),...
            '_',num2str(c_thr2),'_',num2str(replicate_ind),'.txt');
    dlmwrite(filename6,population);
    indices=horzcat(lastbest,index);
    dlmwrite(strcat('indices',num2str(c_thr2),'_',num2str(replicate_ind),'.txt'),indices);
end

