function shortp = gen_shortp(i,DownNodes,pred)
%GEN_SHORTP
%  
shortp=i;
while pred(i)~=0
    shortp=[pred(i),shortp]; 
    i=pred(i);
end

np=length(shortp);  
for m=1:np     
    shortp(m)=DownNodes(shortp(m));
end
    

end

