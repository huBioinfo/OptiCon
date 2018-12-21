function [RNexp_ph,count] = gen_RNexp(GeneID,dm_ph)
%GEN_RNEXP 
%   Detailed explanation goes here
nrow=size(GeneID,1);  
ncol=size(dm_ph,2);  
RNexp_ph=zeros(nrow,ncol);
count=0;    

for i=1:nrow
    vector=strcmp(GeneID(i),dm_ph.RowNames);
    index=find(vector);
    if ~isempty(index)  
        RNexp_ph(i,:)=dm_ph(index,:);  
        count=count+1;
    end
end

end

