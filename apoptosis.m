function outflow=apoptosis(N,noisesd,threshold,apercent,H,K,T)

deltaN_outflow=0*N; 
% to go oncogenic layer by layer to calculate apoptosis per layer
% take one layer, get detection
% calculate whether these cells that are detected are killed
% remove from N

 for i=1:K % we go through each ongogenic layer
    mutsteps=i-1; % this many mutational steps happened already
    % i.e. if the cell is at oncogenic layer dimension 1, it had zero
    % mutations
    layertocheck=N(:,i,:);
    sum_cells=sum(sum(layertocheck)); 
    % we calculate the number of cells that have this much damage
    % then draw # of cells many epsilon (noise) values from a normal
    % distribution, or as implemented here, calculate how the perceived
    % damage is distributed around the mean of actual damage with standard
    % deviation noisesd
    % 1) how many cells are above the apoptosis threshold?
    % 2) how many of them we will kill? (apercent)
    killed=round(sum(normrnd(mutsteps,noisesd,[1,sum_cells])>=threshold)*apercent);
    % if there are more cells in that layer, we randomly pick which
    % specific ones are killed
    if killed > 0
        killed_cells=reshape(ddistss(layertocheck(:)',killed),H,1,T);
        deltaN_outflow(:,i,:)=killed_cells;
    end
 end
outflow=deltaN_outflow;
end