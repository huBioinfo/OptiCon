function [down_relation,down_weight,down_distance] = downstr_weight(adjMatrix,dys_g,RNexp_ph)
%DOWNSTR_WEIGHT
%   dys_g:dysregulated gene.


down_relation=[];
down_weight=[];
k=dys_g;
k1=dys_g;
while 1
    %------------progress bar--------------%
    %fprintf('k_size %d.\n',length(k));
    %--------------------------------------%
    if isempty(k)
        break;
    else
        v=k(1); 
        [relation,weight,k,k1] = gen_downrelation(adjMatrix,v,k,k1,dys_g,RNexp_ph);
        down_relation=[down_relation;relation]; 
        down_weight=[down_weight;weight];
        continue;
    end  
end
down_distance=-1*log2(down_weight)+1;
end

