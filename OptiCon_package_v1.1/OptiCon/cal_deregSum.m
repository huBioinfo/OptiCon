function [effect_deregSum,uni_effectInd] = cal_deregSum(temp,cfinf,RNsamSco_ph,final_tarCFINF)
%CAL_DEREGSUM 


index=find(cfinf(temp,:));
tempCFINF=cfinf(temp,index); 
effectInd=[final_tarCFINF;tempCFINF'];
uni_effectInd=unique(effectInd);

neffect=size(uni_effectInd,1);  
effect_deregSum=0;
for q=1:neffect
    effect_deregSum=effect_deregSum+RNsamSco_ph(uni_effectInd(q));
end

end

