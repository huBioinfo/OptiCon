function [RNsamSco,count] = RNsamScore(GeneID,microsamID,microsamSco)
%RNSAMSCORE 
%   Detailed explanation goes here
n=size(GeneID,1);
RNsamSco=[];
count=0; 
for i=1:n
    %------------progress bar--------------%
    %fprintf('gene %d.\n',i);
    %--------------------------------------%
    vector=strcmp(GeneID(i),microsamID);
    index=find(vector);
    if isempty(index)  
        RNsamSco=[RNsamSco;0];
    elseif microsamSco(index)>0.05  %-cutoff for differential expression
        RNsamSco=[RNsamSco;0];
    else
        RNsamSco=[RNsamSco;-log10(microsamSco(index))]; 
        count=count+1;
    end
end
        

end

