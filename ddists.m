function ind=ddists(p,n)
% sample=ddists(p,n): discrete distribution sample
% n individuals sampled from a discrete distribution given in p
% 'sample' gives the individual numbers, in random order
% E.g. p = [0.2 0.6 0.2], n=4 => ddists(p,n) = (e.g.) [2 1 2 3]
% (Note: p need not sum up to 1)

ind=ones([1 n]);
pp=cumsum(p/sum(p));
rnd=rand([1 n]); 
for i=1:length(pp)
	ind=ind+(rnd>pp(i));
end