function cfinf = generate_cfinf(cf,RNsamSco_ph,RNcis)
%GENERATE_DFINF
%   Detailed explanation goes here
n=size(cf,1);
cfinf=zeros(n,n);  
for i=1:n
    %------------progress bar--------------%
    %fprintf('Gene: %d.\n',i);
    %--------------------------------------%
    index=[];
    j=1; 
    while cf(i,j)~=0
        index=[index,cf(i,j)];  
        
        if RNsamSco_ph(cf(i,j))~=0   
            index=[index,RNcis{cf(i,j)}]; 
        end
        
        j=j+1;  
    end
    index=unique(index);  
    len=length(index);
    cfinf(i,1:len)=index;
end
    
end

