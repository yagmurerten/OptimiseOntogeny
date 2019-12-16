% function for random sampling of the parameter space to start the genetic 
% algoritm or to look for general trends within the parameter space we sample

function [fitness_vals]=shotgun(extmort,cellmort,cancer_danger,targetsize,...
   rounds,withextmort,replicates,extmortthreshold, sg_only, filename,celldeath,strategy_bounds,nof_onco_steps)
    length_strategy=size(strategy_bounds,1);
    writeouts=10; % dimension of the matrix HARDCODED
    % meaning that if we decide to output anything else from the simulation
    % this NEEDS TO BE CHANGED
    
    % if sg_only, it is without the subsequent genetic algorithm and we
    % just want to write out fitness values, instead of keeping them in the
    % memory (so that we can actually simulate a very large population)
    % if we want to create a population for the genetic algortihm, then we
    % need to keep them in the memory 
    if ~sg_only
        fitness_vals=zeros(rounds,length_strategy+length(extmort)+replicates*writeouts);
    end
    % creating the set of strategies to be tested, randomly within the bounds
    % we do this at the start to prevent calling random numbers all the time
    probs=strategy_bounds(:,2)+rand(length_strategy,rounds).*(strategy_bounds(:,3)-strategy_bounds(:,2));
    probs(3,:)=round(probs(3,:)); % telomere length
    probs(6,:)=round(probs(6,:)); % differentiation steps
    % these two *always* need to be rounded, as they determine the
    % dimensions of a matrix
    
    % these are to prevent simulating strategies if they fall below a
    % certain fitness threshold, but we don't want to use this option if we
    % are doing shotgun, relevant for ga function and explained further
    % there
    fitness_thr1=0;fitness_thr2=0;
    for index=1:rounds   
        strategy=probs(:,index)';
        [temp_fit, ~, ~]=simulatestrategy(replicates,...
            strategy,[cellmort,extmort],cancer_danger,targetsize,...
            withextmort,extmortthreshold,celldeath,fitness_thr1,fitness_thr2,nof_onco_steps);
        if sg_only
            dlmwrite(filename,temp_fit,'-append')
        else
            fitness_vals(index,1:end)=temp_fit;
        end
    end
       
end