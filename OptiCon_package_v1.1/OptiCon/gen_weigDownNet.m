function [DownNodes,weigDownNet] = gen_weigDownNet(down_relation,down_distance)
%GEN_WEIGHTDOWNNET
%   Detailed explanation goes here
DownNodes=[down_relation(:,1);down_relation(:,2)];
DownNodes=unique(DownNodes);

n=size(DownNodes,1);
weigDownNet=zeros(n,n);
nedge=size(down_relation,1);
for i=1:nedge
    index_s=find(DownNodes==down_relation(i,1)); 
    index_t=find(DownNodes==down_relation(i,2));
    weigDownNet(index_s,index_t)=down_distance(i);
end

end

