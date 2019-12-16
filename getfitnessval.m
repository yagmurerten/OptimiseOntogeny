% function to calculate fitness taking 
% the extrinsic mortality and per-time-step fitness into account
function expected_fit_onelife=getfitnessval(fitnessvals,extmort)
    t=0:length(fitnessvals)-1; % starting from age 0 makes step 1 individuals always be alive (no time to have died yet)
    expected_fit_onelife=sum(fitnessvals.*exp(-extmort.*t));
end