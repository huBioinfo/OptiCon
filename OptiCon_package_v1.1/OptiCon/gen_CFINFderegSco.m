function cfinf_deregSco = gen_CFINFderegSco(cfinf,RNsamSco_ph,RN_deregSum_hepa)
%GEN_cfinfDEREGSCORE 
%   Detailed explanation goes here
n=size(cfinf,1);
cfinf_deregSco=zeros(n,1);  
for i=1:n
    j=1;
    while cfinf(i,j)~=0  
        cfinf_deregSco(i)=cfinf_deregSco(i)+RNsamSco_ph(cfinf(i,j));   
        j=j+1;
    end
    cfinf_deregSco(i)=cfinf_deregSco(i)/RN_deregSum_hepa;  %convert to fraction.
end

      
end

