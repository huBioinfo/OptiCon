function [best_score,best_tarAdd,best_tarAddCFINF,flag] = get_CF(cfinf,cfinf_deregSco,best_score,best_tarAdd,best_tarAddCFINF,final_tar,final_tarCFINF,RNsamSco_ph,RN_deregSum,RN_WeakInd,flag)
%GET_CF Summary of this function goes here
%   Detailed explanation goes here

controlDereg_Ind=find(cfinf_deregSco);
%search_Ind=controlDereg_Ind;

left=setdiff(controlDereg_Ind,final_tar);  

%nlef=size(left,1);
for p=1:size(left,1)  %-greedy search from ~5959 nodes.
    
      temp=left(p);  %--try to add "temp", which is a GeneID
      %---------desired effect--------------%
%       [effect_deregSum,uni_effectInd] = cal_deregSum(temp,cfinf,RNsamSco_ph,final_tarCFINF);
%       
%       index=find(cfinf(temp,:));
%       tempCFINF=cfinf(temp,index);
      
      tempCFINF_all=cfinf(temp,:);
      tempCFINF=tempCFINF_all(tempCFINF_all~=0); %-create new variable might improve performance.
     
      effectInd=[final_tarCFINF;tempCFINF'];  %-merge OCRs.
      uni_effectInd=unique(effectInd);

%       neffect=size(uni_effectInd,1);  
%       effect_deregSum=0;
%       for q=1:neffect
%           effect_deregSum=effect_deregSum+RNsamSco_ph(uni_effectInd(q));
%       end
      effect_deregSum=sum(RNsamSco_ph(uni_effectInd));  %-vectorized computation.
      DerRat=effect_deregSum/RN_deregSum;  %between 0--1
      
      %-----------undesired effect-------------%
      eff_WeakInd=intersect(uni_effectInd,RN_WeakInd);  
      WeakRat=length(eff_WeakInd)/length(RN_WeakInd);     %between 0--1
      
      %--------objective function--------------------------%
      DWscore=DerRat-WeakRat;
      
      %-----greedy search-------%
      if DWscore>best_score  
              best_score=DWscore;
              best_tarAdd=temp;
              %ind=find(cfinf(best_tarAdd,:));
              %best_tarAddCFINF=cfinf(best_tarAdd,ind);
              
              best_tarAddCFINF_all=cfinf(best_tarAdd,:);
              best_tarAddCFINF=best_tarAddCFINF_all(best_tarAddCFINF_all~=0);  %-create new variable might improve performance.
              
              flag=1;
      end
      
end  %-end of for.

end
