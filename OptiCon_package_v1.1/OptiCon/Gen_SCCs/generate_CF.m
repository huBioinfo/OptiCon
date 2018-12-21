function cf = generate_CF(tempCR)
%GEN_CF 
%   Detailed explanation goes here

n=size(tempCR,1);
cf=zeros(n,n); 
for i=1:n
    ind=find(tempCR(i,:));
    len=length(ind);
    cf(i,1:len)=ind;
end
    
end

