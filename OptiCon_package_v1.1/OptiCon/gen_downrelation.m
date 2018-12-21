function [relation,weight,k,k1] = gen_downrelation(G,v,k,k1,dys_g,RNexp_ph)
%GEN_DOWNRELATION

relation=[];
weight=[];
a=find(G(v,:)==1);
if ~isempty(a)
    for i=1:length(a)
        relation=[relation;[v,a(i)]];
        w=comp_w(v,a(i),dys_g,RNexp_ph);
        weight=[weight;w];     
        
        if ~ismember(a(i),k1)
            k=[k,a(i)]; 
            k1=[k1,a(i)];
        end
    end
end
k(1)=[]; 
end

