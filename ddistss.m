function count=ddistss(p,n)
% sample=ddistss(p,n): discrete distribution sample sums
% n individuals sampled from a discrete distribution given in p
% 'sample' is a vector that gives the number of different types obtained
% e.g.: sample=ddistss([0.3 0.5 0.2],10) may give [3 5 2] or perhaps [2 6 2] or...
count=ones([1 length(p)]);
ind=ddists(p,n);
for i=1:length(p)
	count(i)=sum(ind==i);
end
