function Permute = gen_randInitialPoint(randNum)
%GEN_RANDINITIALPOINT Summary of this function goes here
%   Detailed explanation goes here

load initialPoints
n=size(NodeMMrate005,1); 
rng('shuffle');
out1 = randperm(n);
Permute=out1(1:randNum); 
end
