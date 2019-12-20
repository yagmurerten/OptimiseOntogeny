function [mat_time,lifespan,deathcause,fitness,cell_types,killed_cells]=onelife(strategy,mort,...
    cancerdanger,maturitycriterion,withextmort,extmortthreshold,celldeath,nofonco)

% strategy: 
% 1 - probability of asymmetric division; P

% 2 - conditional probability of differentiation GIVEN there was a
%   symmetric division (= given the option of asymmetric division was 
%   not used); Q

% 3 - telomere length: how many telomere layers does the organism
%   use, dim(1) of the N matrix; Hayflick limits, H

% 4: apoptosis threshold; A

% 5: apoptosis percent; S

% 6: number of differentiation levels; T

% 7: division propensities: relative use of different 
%   cells for division, compared with the previous level; X

%%%%%%%%%

% mort(1): per-cell death rate (v)
% mort(2) whole-organism death rate, both expressed per time step (mu)

%%%%%%%%%

% maturitycriterion is the target number of fully differentiated cells

%%%%%%%%%

% N will contain:
% rows = H-telomeres (1 = longest, strategy(3) = so short you won't divide 
%   any more)
% columns = K-damage level (1 = no damage, nof_oncosteps = so damaged you're
%   recorded dead because of cancer [assuming 1 such cell is considered 
%   bad enough])
% layers = T-differentiation level (1 = stem cell, T = fully differentiated 
%   cell)

H=strategy(3); % how many telomere-related layers there will be, 
               % i.e. how many times a cell can divide (Hayflick limit)             
K=nofonco; % how many steps until cancer
apoptosis_thr=strategy(4); % A
apoptosis_percent=strategy(5); % S

% this implements noise around damage detection
noisesd=0.5;

% to count how many cells were killed due to each cause
killedcells_turn=0;killedcells_apop=0;killedcells_H=0;

k=0.05; % parameter for the growth trajectory

minpercenttissue=0.8; 
% the organism needs to have at least 80% of its target terminally differentiated cells
% after maturation - related to death cause #5
% intituitvely, we would not hit to this limit as before this, the organism
% would die due to lack of stem cells
% unless the apoptosis rate or baseline cell death rate is perhaps really high

T=strategy(6); 
% differentiation layers; how many steps there are from the stem cell to tissue, 
% assuming the first layer consists of the most stem cells and the last
% layer consists of full differentiated cells (with no division
% possibility)

divisionpropensities=strategy(7).^linspace(1,T-1,T-1);
% what is the division propensity of each differentiation layer compared to
% the most stem cells (which has the propensity 1). We assume that this is
% 0 for the last layer
% i.e. this is 1/straregy(7) for all the layer, compared to the next
% one, except for the terminally differentiated layer

N=zeros([H K T]); 

N(1,1,1)=1; % the zygote has been created! 

bodysize=N(1,1,1); % number of total cells is the 'body size'
tissue=0; % number of cells in totally differentiated state (tissues)
dividibles=[N(1,1,1) zeros(1,T-2)]; % these cells can divide
% note - code requires at least 1 intermediate step between stem and tissue
alive=1;mature=0;t=0;mat_time=0;

extmort=mort(2);
endtimemort=ceil(-log(extmortthreshold)/min(extmort));
% probability that the organism lives until this is = extmortthreshold (0.001
% is the default) given its extrinsic mortality

pertfitness=0;fitness_vals=zeros(1,endtimemort);

while alive
    t=t+1;
    % we calculate fitness at the beginning of each round
    fitness_vals(1,t)=pertfitness;
    if withextmort % if we want extrinsic mortality to be applied stochastically througout sim. time
        if rand(1)<extmort
            % whole organism dies of extrinsic causes
            alive=0; deathcause=1; 
            break;
        end
    end
    if ~mature && alive % not mature yet, let's grow
         % this kind of leads to vanB. growth        
        tissue=sum(sum(N(:,:,end)));
        % how many cells are needed
        total_divisions=ceil(k*(maturitycriterion-tissue));
        
        % if the stem cells are depleted before maturation, it is
        % dead --> it won't mature anyway
        if sum(dividibles)==0 
            alive=0; deathcause=2;
            break;
        end
        if alive && total_divisions > 0 
            N=divisionsandmutations(N,strategy,total_divisions,...
                dividibles,divisionpropensities,cancerdanger,H,K);
        end
        % now check if maturity has been reached
        tissue=sum(sum(N(:,:,end)));
        % maturitycriterion=target number of differentiated cells that
        % means the organism has matured
        if tissue>=maturitycriterion
            % now it is mature!
            mature=1;mat_time=t;
            % update your info on total body size
            bodysize=sum(sum(sum(N,3)));bodysize_at_maturity=bodysize;
        end
    else % what to do when mature
        bodysize=sum(sum(sum(N,3)));tissue=sum(sum(N(:,:,end)));
        pertfitness=tissue/bodysize;        
        % counteract the loss of cells that may have happened
        total_divisions=ceil(k*(maturitycriterion-tissue));
        
        % this is also a little bit harsh -> you die immediately as soon as
        % you don't have stem cells
        % but might also be reasonable considering the fact that some
        % tissues (e.g. blood)should be renewed every day and rely on the
        % stem cells being present
        % also: you will keep losing your tissue, then it becomes a matter
        % of time when you die, i.e. you cannot live a long time without
        % any stem cell anyway, but of course another idea is that you can
        % divide your stem cells to replenish        
        if sum(dividibles)==0 
            alive=0; deathcause=8;
            break;
        elseif total_divisions > 0 % do this part only if you need divisions
            N=divisionsandmutations(N,strategy,total_divisions,...
                dividibles,divisionpropensities,cancerdanger,H,K);
        end
    end
    
    % if you managed to survive the first two death causes:

    % now some cells also die randomly (baseline cell death, sort of like 
    % cell extrinsic mortality)
    if celldeath
        deadcells=binornd(N,mort(1)); N=N-deadcells;
        killedcells=sum(sum(sum(deadcells)));
        killedcells_turn=killedcells_turn+killedcells;
    end
    % additionally, max-telomere cells always die (note: this means they 
    % won't contribute to body size which is computed later, see below)
    killedcells_H=killedcells_H+sum(sum(N(H,:,:)));
    N(H,:,:)=0;

    % apoptosis clears out some cells perceived as above the damage
    % threshold, but there is some noise (implemented above)
    if celldeath
        killed=apoptosis(N,noisesd,apoptosis_thr,apoptosis_percent,H,K,T);
        N=N-killed;
        killedcells_apop=killedcells_apop+sum(sum(sum(killed)));
    end
    N(N<0)=0;
    
    % now check if you've got cancer, even 1 cell in
    % the O (maximum onco-state) category kills you!
    cancercells=sum(sum(N(:,end,:)));
    if cancercells>0
        alive=0; deathcause=3;
        break;
    end

    % update your info on dividibles
    for i=1:length(dividibles)
        dividibles(i)=sum(sum(N(1:H-1,1:K,i)));
    end

    % update your info on total body size and tissue size
    bodysize=sum(sum(sum(N,3)));tissue=sum(sum(N(:,:,end)));

    % death causes related to cell numbers = losing all cells, 
    % imbalance between differentiated cells and stem cells

    if bodysize==0 % the organism actually vanished (this cause should only be applied if not mature yet)
        alive=0; deathcause=4; 
        break;
    end
  
    % Then once you mature and reach your 'final' body size, 
    % you cannot go below a certain percent of that; i.e. you need a
    % certain number of tissue cells to maintain your functions
    
    if mature && tissue < minpercenttissue*maturitycriterion
        alive=0; deathcause=5;
        break;
    end  

    if t>=ceil(endtimemort)
        alive=0; 
        if ~mature
            deathcause=6;
            break;
        else
            deathcause=7; 
            % this is when the organism does not die but 
            % the probability it survived beyond this point is 0.001
            break;
        end 
    end
end
lifespan=t;
cell_types=horzcat(dividibles(1), sum(dividibles(2:end)), tissue);
killed_cells=horzcat(killedcells_turn,killedcells_H,killedcells_apop);
fitness=getfitnessval(fitness_vals,extmort);