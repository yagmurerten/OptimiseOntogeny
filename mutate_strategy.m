function new_strategy=mutate_strategy(unmutatedstrategy,strategy_bounds)

length_strategy=length(unmutatedstrategy);
% new random strategy
probs=strategy_bounds(:,2)+rand(length_strategy,1).*(strategy_bounds(:,3)-strategy_bounds(:,2));
% how far the resident strategy is from the new random strategy
difference=probs'-unmutatedstrategy;
% percent we move towards the new random strategy
% calculated with a gamma distribution, truncated to be between 0 and 1 as
% it is a percent change
pd=makedist('gamma','a',0.25,'b',0.1);pd_t=truncate(pd,0,1);
percentchange=random(pd_t,1,length_strategy);
new_strategy=unmutatedstrategy+difference.*percentchange;
% things that need to be integers are rounded
new_strategy(3)=round(new_strategy(3));
new_strategy(6)=round(new_strategy(6));
end