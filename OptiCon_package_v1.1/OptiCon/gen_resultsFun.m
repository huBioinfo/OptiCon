function [] = gen_resultsFun(initial)
%GEN_RESULTSFUN Summary of this function goes here
%   Detailed explanation goes here

best_score=0;
best_tarAdd=[];
best_tarAddCFINF=[];
best_m=[];
oldBest_score=0;
final_tar=[];
final_tarCFINF=[];
final_m=[];
num=0;
final_score=[];
final_CFINFtable=[];
incre_rate=[];

%----Set initial search point-----------%
load InitialPoints
num=1;
final_tar=[final_tar;NodeMMrate(initial,1)];  
final_m=[final_m;NodeMMrate(initial,2)];  
best_score=NodeMMrate(initial,3);

str3=['load CFINF_hepa_' num2str(final_m) ];
eval(str3);
str4=['cfinf=cfinf_hepa_' num2str(final_m) ';'];
eval(str4);
index=find(cfinf(final_tar,:));
final_tarCFINF=cfinf(final_tar,index);
final_tarCFINF=final_tarCFINF';

oldBest_score=best_score;
final_score=[final_score;oldBest_score]; 
final_CFINFtable=[final_CFINFtable;{final_tarCFINF}];
incre_rate=[incre_rate;(best_score-0)/0];

load sam_hepa

while 1

for j=0:999
%     %------------progress bar--------------%
%     fprintf('CF %d.\n',j); %-this cost time
%     %--------------------------------------%

%     str3=['load CFINF_hepa_' num2str(j) ];
%     eval(str3);
    CFINF_file=strcat('CFINF_hepa_',num2str(j));
    load(CFINF_file); %-checked. better than eval.
    
    str4=strcat('cfinf=cfinf_hepa_',num2str(j),';');
    eval(str4);  %-seems no alternative way.
    
    str5=strcat('cfinf_deregSco=cfinf_deregSco_hepa_',num2str(j),';'); %determine search space
    eval(str5);  %-seems no alternative way.
    
    flag=0;
    [best_score,best_tarAdd,best_tarAddCFINF,flag]=get_CF(cfinf,cfinf_deregSco,best_score,best_tarAdd,best_tarAddCFINF,final_tar,final_tarCFINF,RNsamSco_hepa,RN_deregSum_hepa,RN_WeakInd_hepa,flag);
    if flag==1 
        best_m=j;
    end
    clearvars -except j best_score best_tarAdd best_tarAddCFINF best_m oldBest_score final_tar final_tarCFINF final_m num final_score final_CFINFtable incre_rate RNsamSco_hepa RN_deregSum_hepa RN_WeakInd_hepa initial
end

if (best_score-oldBest_score)/oldBest_score>=0.05 
    incre_rate=[incre_rate;(best_score-oldBest_score)/oldBest_score];
    oldBest_score=best_score;
    final_tar=[final_tar;best_tarAdd]; 
    num=num+1;
    %------------progress bar--------------%
    fprintf('targetNo %d.\n',num);
	fix(clock)
    %--------------------------------------%
    final_tarCFINF=[final_tarCFINF;best_tarAddCFINF'];
    final_tarCFINF=unique(final_tarCFINF);
    
    final_m=[final_m;best_m]; 
    final_score=[final_score;oldBest_score];
    final_CFINFtable=[final_CFINFtable;{final_tarCFINF}]; 
else
    break;
end
end  %end of while.

str9=['save finTherapTar_relax' num2str(initial) ' oldBest_score final_tar final_tarCFINF final_m num final_score final_CFINFtable incre_rate'];
eval(str9);

end

