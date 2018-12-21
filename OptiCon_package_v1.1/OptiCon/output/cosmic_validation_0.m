%---This is to compare OCR with recurrent mutant "cancer" genes in a specific cancer type of interest------%

%--import cancer gene census-----V79--568 genes---%
fid1=fopen('CancerCensus_En.txt');
Census_En=textscan(fid1,'%s');
Census_En=Census_En{1};
fclose(fid1);  
save cosmicValid Census_En

%---import recurrently mutated genes in a specific cancer type------%
clear
fid2=fopen('RecurMutant_entrez.txt'); 
RecurMutant_En=textscan(fid2,'%s');
RecurMutant_En=RecurMutant_En{1};
fclose(fid2);  
save cosmicValid RecurMutant_En -append

%---identify cancer genes recurrently mutated genes in a specific cancer type----%
clear
load cosmicValid
RecurCensus_En=intersect(Census_En,RecurMutant_En);
save cosmicValid RecurCensus_En -append
clear
load ../RNgeneID
load cosmicValid
knownMu_inNet=intersect(RecurCensus_En,GeneID);
knownMu_inNet_Sym=[];
n=size(knownMu_inNet,1);
for i=1:n  
     vector=strcmp(knownMu_inNet(i),symbol2entrez_Integ(:,2));
     index=find(vector);
     Sym=symbol2entrez_Integ(index,1);
     knownMu_inNet_Sym=[knownMu_inNet_Sym;Sym];
end
save cosmicValid knownMu_inNet_Sym knownMu_inNet -append

exit
