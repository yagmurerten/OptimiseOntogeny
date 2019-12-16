% the function to swap optimised strategies trait by trait, in order to see
% if the traits are coadapted

function disrupt_strategy(index,bs_index)
    % note that this requires an input files derived from the raw data
    % generated in baseline simulations, that contains the optimised
    % strategies
    t1 = datetime('now','Format','dd-MMM-yyyy HH:mm:ss.SSS');
    s=second(t1,'secondofday');
    seed=s*index*bs_index; 
    rng(seed);

    A = importdata('strategies_0.01_1-4_300.txt');
    % the strategies optimised 300 rounds at extmort 0.01 for 4 replicates
    % previously

    extmort= A(1,4);cellmort=A(1,2);cancer_danger=A(1,3);
    bodysizes=unique(A(1:end,5)); % bodysizes from the simulation
    
    % HARDCODED, need change these in case it is necessary
    celldeath=1;withextmort=0;extmortthreshold=0.001;
    fitness_thr1=-1;fitness_thr2=0;replicates=10;
    nof_onco_steps=4;

    % the seven traits that will be swapped
    traits=[string('prob_a'),string('prob_d'),string('tel'),...
        string('onco_thr'),string('onco_prob'),string('dif_l'),string('div_prop')];
    trait=traits(index);
    subset_bs=bodysizes(bs_index);
    filename=strcat('data_disrupt/mixed_',trait,'_',num2str(subset_bs),'.txt');

    % get the four strategies optimised for that body size
    subset_A=A(A(:,5)==subset_bs,:);

    for j=1:size(subset_A)
        focal_rep_index=j;
        focal_strategy=subset_A(focal_rep_index,:);
        other_strategies=subset_A(subset_A(:,1)~=focal_rep_index,:);
        f_str=focal_strategy(6:12);
        % get the mean fitness with the optimised version as a control
        [temp_fit, ~, ~]=simulatestrategy(replicates,...
                        f_str,[cellmort,extmort],cancer_danger,subset_bs,...
                        withextmort,extmortthreshold,celldeath,fitness_thr1,fitness_thr2,nof_onco_steps);
        temp_fit2=horzcat(focal_rep_index,focal_rep_index,cellmort,cancer_danger,extmort,subset_bs,temp_fit);
        dlmwrite(filename,temp_fit2,'-append');

        for i=1:size(other_strategies)
               % get the 'disrupted' versions by mixing in the 'trait'
               % we are looking at
               dist_str=other_strategies(i,6:12);
               dist_str_ind=other_strategies(i,1);
               mixed_str=get_trait(f_str,dist_str,trait);
               [temp_fit, ~, ~]=simulatestrategy(replicates,...
                        mixed_str,[cellmort,extmort],cancer_danger,subset_bs,...
                        withextmort,extmortthreshold,celldeath,fitness_thr1,fitness_thr2,nof_onco_steps);
                temp_fit2=horzcat(focal_rep_index,dist_str_ind,cellmort,cancer_danger,extmort,subset_bs,temp_fit);
                dlmwrite(filename,temp_fit2,'-append');
        end
    end
end

