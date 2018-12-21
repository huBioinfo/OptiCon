function [ICV_valueTable,ICV_SPtable] = gen_cis(RNsamSco_ph,adjMatrix,RNexp_ph)
%GEN_CIS Summary of this function goes here
%   Detailed explanation goes here
n=size(RNsamSco_ph,1);  
%RNcis=cell(n,1); 
ICV_valueTable=cell(n,1);
ICV_SPtable=cell(n,1);
for i=1:n
    %------------progress bar--------------%
    fprintf('Gene: %d.\n',i);
    %--------------------------------------%
  
    if RNsamSco_ph(i)~=0   
            dys_g=i;
            [Imp,SP] = gen_coreAffMod(adjMatrix,dys_g,RNexp_ph);

            ICV_valueTable{i}=Imp;
            ICV_SPtable{i}=SP;
    end
   
end


end

