clear
load SynergyScore
synScoRand_summary=zeros(size(s1Rand_summary,1),1);
for i=1:size(s1Rand_summary,1)
    synScoRand_summary(i)=((s1Rand_summary(i)-min_s1)/(max_s1-min_s1))*((s2Rand_summary(i)-min_s2)/(max_s2-min_s2)); %Normalized before multiplying.
end
save SynergyScore synScoRand_summary -append

clear
load SynergyScore
nCombTar=size(synergySco_sorted,1);
synergySco_sortedP=ones(nCombTar,1);
nRand=size(synScoRand_summary,1);
for i=1:nCombTar
    %------------progress bar--------------%
    fprintf('combTar %d.\n',i);
    %--------------------------------------%
    sco_largerCount=0;
    for j=1:nRand
        if synScoRand_summary(j)>=synergySco_sorted(i)
            sco_largerCount=sco_largerCount+1;
        end
    end
    synergySco_sortedP(i)=sco_largerCount/nRand;
end
save SynergyScore synergySco_sortedP -append

clear
load SynergyScore
[h, crit_p, adj_ci_cvrg, synergySco_sortedP_adj]=fdr_bh(synergySco_sortedP,0.05,'pdep','yes');
%-'pdep' is BH-adjusted method.
save SynergyScore synergySco_sortedP_adj -append

clear
load SynergyScore
OptCombTar_sorted_Sym=[];
OptSynergySco_sorted=[];
OptSynergySco_sortedP=[];
for i=1:size(synergySco_sortedP_adj,1)
    if synergySco_sortedP_adj(i)<=0.05  %-pValue cutoff
        OptCombTar_sorted_Sym=[OptCombTar_sorted_Sym;combTar_sorted_Sym(i,:)];
        OptSynergySco_sorted=[OptSynergySco_sorted;synergySco_sorted(i)];
        OptSynergySco_sortedP=[OptSynergySco_sortedP;synergySco_sortedP_adj(i)];
    end
end
save SynergyScore OptCombTar_sorted_Sym OptSynergySco_sorted OptSynergySco_sortedP -append


