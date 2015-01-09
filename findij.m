% find the i and j of a certain site
% i and j are parameters n=2^(i-1)(2j+1)
% this function is for HN3_WL_ToyModel.m
function [i,j] = findij(n)
i=1;
while mod(n,2)==0
     n=n/2;
     i=i+1;
end
j=(n-1)/2;
end
