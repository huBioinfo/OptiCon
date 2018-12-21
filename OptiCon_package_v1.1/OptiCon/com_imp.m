function [Imp,SP] = com_imp(weigDownNetSpar,DownNodes,dys_g)
%COM_IMP
%   Detailed explanation goes here
u=find(DownNodes==dys_g); 
[d,pred]=dijkstra_sp(weigDownNetSpar,u);

Imp=1./d;

n=size(Imp,1);
SP=cell(n,1);
for i=1:n
    SP{i}=gen_shortp(i,DownNodes,pred);
end

end

