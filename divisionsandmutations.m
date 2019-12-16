function N=divisionsandmutations(N,strategy,target,dividibles,divisionpropensities,cancerdanger,T,O)

deltaN_outflow=0*N; deltaN_inflow=0*N;
diffprob=strategy(2);
asymprob=strategy(1);
% now divide the divisions into differentiation categories as dictated by the numbers in strategy(5) onwards
% the following will contain total divisions that differ in differentiatedness
% I divide the target into two, because each division will result in two
% cells
nof_divisions=ddistss(dividibles.*divisionpropensities,target);
if sum(nof_divisions~=0)
    for i=1:length(nof_divisions) % i refers to the layer at hand; we already know the number of divisions to be performed by this "layer"
        potentialcells=N(1:end-1,:,i);
        divisions=min(potentialcells,reshape(ddistss(potentialcells(:)',nof_divisions(i)),T-1,O));
        asymmetricdivisions=binornd(divisions,asymprob);
        symmetricdivisions=divisions-asymmetricdivisions;
        symmdiffdivisions=binornd(symmetricdivisions,diffprob);
        symmnodiffdivisions=symmetricdivisions-symmdiffdivisions;
        deltaN_outflow(1:T-1,1:O,i)=-divisions; % these ones disappear (and appear, doubled, in the inflow)
        % they potentially move up one layer (differentiate) or stay
        % put; but regardless, they move down along the T axis. Note
        % that movements along the O axis come a bit later
        deltaN_inflow(2:T,1:O,i+1)=deltaN_inflow(2:T,1:O,i+1)+2*symmdiffdivisions+asymmetricdivisions;
        deltaN_inflow(2:T,1:O,i)=deltaN_inflow(2:T,1:O,i)+2*symmnodiffdivisions+asymmetricdivisions;
        alreadyDivided=sum(sum(sum(deltaN_inflow)));
        enoughCells=alreadyDivided >= sum(nof_divisions);
        if enoughCells
            break;
        end
    end
    % then the mutations; any cell throughout the 'inflow' can be mutated. Note that this is done once 
    % the whole inflow otherwise is known.
    % Note also that restricting this to deltaN means that you need to divide in order to mutate. 
    % One could use the same structure on N itself, if one wanted to relax this assumption.
    % Note also that I'm assuming that if a rightmost (already cancerous) cell mutates, it still is cancerous...
    mutations=binornd(deltaN_inflow,cancerdanger);
    sum_mut=sum(sum(sum(mutations)));
    deltaN_inflow(:,1:O-1,:)=deltaN_inflow(:,1:O-1,:)-mutations(:,1:O-1,:);
    deltaN_inflow(:,2:O,:)=deltaN_inflow(:,2:O,:)+mutations(:,1:O-1,:);
end
% now ready to update N
N=N+deltaN_outflow+deltaN_inflow; % note the '+' in outflow; we've kept it negative above, so shouldn't use it negatively here