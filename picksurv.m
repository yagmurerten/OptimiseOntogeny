function ind=picksurv(S,n)
% ind=picksurv(S,n)
% Picks n survived individuals if their survivals relate like S
% Picks all individuals if n>=length(S) - even if some S are 0
if length(S)<n 
	ind=1:length(S); 
else
	p=rand(size(S));
	indsuccess=S./p;
	[~,ind]=sort(-indsuccess);
	ind=sort(ind(1:n));
end

