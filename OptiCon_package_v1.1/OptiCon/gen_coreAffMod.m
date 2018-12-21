function [Imp,SP] = gen_coreAffMod(adjMatrix,dys_g,RNexp_ph)
%GEN_COREAFFEMOD 
%   Detailed explanation goes here
[down_relation,down_weight,down_distance] = downstr_weight(adjMatrix,dys_g,RNexp_ph); 
if ~isempty(down_relation) 
    [DownNodes,weigDownNet] = gen_weigDownNet(down_relation,down_distance);
    weigDownNetSpar=sparse(weigDownNet); 
    [Imp,SP] = com_imp(weigDownNetSpar,DownNodes,dys_g);
else
    Imp=inf; 
    SP={dys_g};
end

end

