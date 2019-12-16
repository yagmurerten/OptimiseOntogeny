% get the mean fitness of a strategy over n replicates for different
% extrinsic mortality values
function fitness=calculatemeanfitness(fitnessmatrix,replicates)
    scaled=(1/replicates).*fitnessmatrix;
    fitness=sum(scaled,1);
end