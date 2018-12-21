function adjMatrix = gen_adjM(GeneID,kegg_RN)
%GEN_ADJM 
%   Detailed explanation goes here
nvet=size(GeneID,1);
adjMatrix=zeros(nvet,nvet);
nedge=size(kegg_RN,1);
for i=1:nedge
    vector1=strcmp(kegg_RN(i,1),GeneID);
    source=find(vector1);
    vector2=strcmp(kegg_RN(i,2),GeneID);
    target=find(vector2);
    adjMatrix(source,target)=1;
end


end

