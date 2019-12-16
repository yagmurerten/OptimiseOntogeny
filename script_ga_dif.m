% this is the function ran via the python script in the cluster
% function ga.m can otherwise be ran individually to check out different
% things

% I added the seed_index in the code for the random number generator,
% because for some reason Matlab does weird parallelisation within itself
% and if things are ran at the same time with rng('shuffle'), one might end
% up with the same series of random numbers

function script_ga_dif(extmort,targetsize,cancer_danger,cellmort,seed_index,num_reps_b,num_reps_e, c_thr1,c_thr2)
    t1 = datetime('now','Format','dd-MMM-yyyy HH:mm:ss.SSS');
    s=second(t1,'secondofday');
    seed=num_reps_b*s*seed_index; 
    seedfilename=strcat('seed_',num2str(c_thr2),'.txt');
    dlmwrite(seedfilename,seed);
    % number of parallel optimisation rounds
    num_reps=num_reps_e-num_reps_b+1;
    s_index=1; % this changes for every replicate ran
    % for more robust random number parallelisation
    for j=num_reps_b:num_reps_e
        % we are creating this, so that in case there are 2+ replicates,
        % the random numbers are independent from each other
        stream=RandStream.create('mrg32k3a','NumStreams',num_reps,'Seed',seed,'StreamIndices',s_index);   
        RandStream.setGlobalStream(stream);
        ga_difevol(extmort,targetsize,cancer_danger,cellmort,j,c_thr1,c_thr2)
        s_index=s_index+1;
    end
end