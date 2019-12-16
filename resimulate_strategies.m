% the function to resimulate strategies evolved for a given body size in
% different body sizes

function resimulate_strategies(index)
    % note that this requires an input files derived from the raw data
    % generated in baseline simulations, that contains the optimised
    % strategies
    t1 = datetime('now','Format','dd-MMM-yyyy HH:mm:ss.SSS');
    s=second(t1,'secondofday');
    seed=s*index; 
    rng(seed);
    seedfilename=strcat('seed_',num2str(index),'.txt');
    dlmwrite(seedfilename,seed);

    A = importdata('strategies_0.01_1-4_300.txt');
    extmort= A(index,4);cellmort=A(index,2);cancer_danger=A(index,3);
    bodysizes=unique(A(1:end,5)); % bodysizes from the simulation
    
    % HARDCODED, need change these in case it is necessary
    celldeath=1;withextmort=0;extmortthreshold=0.001;
    fitness_thr1=-1;fitness_thr2=0;replicate=A(index,1);replicates=10;
    nof_onco_steps=4;
    original_bs = A(index,5);
    strategy=A(index,6:end);

    filename=strcat('resimulated_strategies',num2str(index),'.txt');

    for i=1:length(bodysizes)
       targetsize=bodysizes(i);
       [temp_fit, ~, ~]=simulatestrategy(replicates,...
                strategy,[cellmort,extmort],cancer_danger,targetsize,...
                withextmort,extmortthreshold,celldeath,fitness_thr1,fitness_thr2,nof_onco_steps);
        temp_fit2=horzcat(replicate,cellmort,cancer_danger,extmort,original_bs,targetsize,temp_fit);
        dlmwrite(filename,temp_fit2,'-append'); 
    end
   
end

