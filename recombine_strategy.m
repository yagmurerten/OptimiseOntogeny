% the new recombination function

function next_strategy=recombine_strategy(new_strategy,gametes,ind)
    length_strategy=length(new_strategy);
    midpoint=randi([1,length_strategy-1],1); % random mid-point for crossover
    lenpop=size(gametes);
    temp_parent_index=ind;
    while temp_parent_index==ind
        temp_parent_index=randi([1 lenpop(1)],1,1);
    end % we want to pick a different parent compared to the first parent
    temp_parent=gametes(temp_parent_index,:);
    % do the cross-over
    
    % is this random pick actually necessary, seeing that the mid-point is already
    % random?
    if rand(1,1) > 0.5
        parent1=temp_parent;parent2=new_strategy;
    else
        parent1=new_strategy;parent2=temp_parent;
    end
    next_strategy=horzcat(parent1(1:midpoint),parent2(midpoint+1:end));
end