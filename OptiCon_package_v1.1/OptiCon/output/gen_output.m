clear
identifyMultiBest_considerFinal

clear
%------------progress bar--------------%
fprintf('Computing mutation scores...\n');
%--------------------------------------%
cosmic_validation

clear
%------------progress bar--------------%
fprintf('Computing crosstalk scores...\n');
%--------------------------------------%
gen_interactDensity

%------------progress bar--------------%
fprintf('Computing synergy scores...\n');
%--------------------------------------%
cmp_SynScore

%------------progress bar--------------%
fprintf('Computing p-values of synergy scores...\n');
%--------------------------------------%
%-Normalization using Min-Max method.
clear
D=dir('randSynScore*.mat');
s1Rand_summary=[];
s2Rand_summary=[]; %-same row order with s1Rand_summary.
for i=1:size(D,1)
    fileName=D(i).name;
    load(fileName);
    s1Rand_summary=[s1Rand_summary;s1Rand];
    s2Rand_summary=[s2Rand_summary;s2Rand];
    clearvars -except D i s1Rand_summary s2Rand_summary
end
max_s1=max(s1Rand_summary);
min_s1=min(s1Rand_summary);
max_s2=max(s2Rand_summary);
min_s2=min(s2Rand_summary);
save SynergyScore s1Rand_summary s2Rand_summary max_s1 min_s1 max_s2 min_s2 -append

clear
load SynergyScore
synergySco_orig=cell(size(s1,1),1);
synergySco_origBest=zeros(size(s1,1),1);
origOCRbestIndex=zeros(size(s1,1),1); %--origOCRbestIndex stores the index of OCR combinations for each OCN pair, same outer order with combTar_OCR.
for i=1:size(s1,1)
    %------------progress bar--------------%
    fprintf('combTar_SynScoNorm %d.\n',i);
    %--------------------------------------%
    synergySco_orig{i}=zeros(size(s1{i},1),1);  %-initiate as "0".
    for j=1:size(s1{i},1)
       synergySco_orig{i}(j)=((s1{i}(j)-min_s1)/(max_s1-min_s1))*((s2{i}(j)-min_s2)/(max_s2-min_s2)); %Normalized before mutiplying.
    end
    [synergySco_origBest(i),origOCRbestIndex(i)]=max(synergySco_orig{i}); %--synergySco_origBest has same order with origOCRindex.
end
save SynergyScore synergySco_orig synergySco_origBest origOCRbestIndex -append

%--Below is the sorted version of synergySco_origBest.
clear
load SynergyScore
load MultiBest_considerFinal   %--synergySco_origBest also has same order with combTar.
[synergySco_sorted,combTar_index]=sort(synergySco_origBest,'descend');
combTar_sorted=combTar(combTar_index,:); %--combTar_sorted has same order with synergySco.
combTar_sorted_multiOCR=combTar_OCR(combTar_index,:);
sorted_OCRbestIndex=origOCRbestIndex(combTar_index,:);
combTar_sorted_OCR=cell(size(combTar_sorted_multiOCR,1),1);
for i=1:size(combTar_sorted_multiOCR,1)
    combTar_sorted_OCR{i}=combTar_sorted_multiOCR{i}(sorted_OCRbestIndex(i),:);
end
save SynergyScore synergySco_sorted combTar_sorted combTar_sorted_OCR combTar_index -append

clear
load SynergyScore
load ../RNgeneID 
combTar_sorted_Sym=[];
for i=1:size(combTar_sorted,1)
    En1=GeneID(combTar_sorted(i,1));
    vector1=strcmp(En1,symbol2entrez_Integ(:,2));
    index1=find(vector1);
    Sym1=symbol2entrez_Integ(index1,1);

    En2=GeneID(combTar_sorted(i,2));
    vector2=strcmp(En2,symbol2entrez_Integ(:,2));
    index2=find(vector2);
    Sym2=symbol2entrez_Integ(index2,1);

    combTar_sorted_Sym=[combTar_sorted_Sym;[Sym1,Sym2]];
end
save SynergyScore combTar_sorted_Sym -append

cmp_SynScore_EmpP

%----------------Write output to files------------------------%
clear
load SynergyScore
fid=fopen('OCN_pairs.txt','w');
fprintf(fid, '%s\t%s\t%s\t%s\t%s\n', 'OCN_1', 'OCN_2', 'SynergyScore', 'pValue', 'BH-adjusted_pValue');
for i=1:size(combTar_sorted_Sym,1)
   fprintf(fid, '%s\t%s\t%f\t%f\t%f\n', combTar_sorted_Sym{i,1},combTar_sorted_Sym{i,2},synergySco_sorted(i),synergySco_sortedP(i),synergySco_sortedP_adj(i));
end
fclose(fid);

%------------progress bar--------------%
fprintf('\n\n\nAll analysis is done! Please refer to "OCN_pairs.txt" for identified synergistic key regulators.\n');
%--------------------------------------%


