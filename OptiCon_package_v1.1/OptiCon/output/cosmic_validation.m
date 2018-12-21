%---This is to compare OCR with recurrent mutant "cancer" genes------%

%--import cancer gene census-----V79--568 genes---%
fid1=fopen('CancerCensus_En.txt');
Census_En=textscan(fid1,'%s');
Census_En=Census_En{1};
fclose(fid1);  
save cosmicValid Census_En

%---import recurrently mutated genes in a specific cancer------%
clear
fid2=fopen('RecurMutant_entrez.txt'); 
RecurMutant_En=textscan(fid2,'%s');
RecurMutant_En=RecurMutant_En{1};
fclose(fid2);  
save cosmicValid RecurMutant_En -append

%---identify cancer genes recurrently mutated genes in a specific cancer----%
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

clear
load MultiBest_considerFinal
load ../RNgeneID
load cosmicValid
n=size(combTar_OCR,1);
combTar_OCRmu3=cell(n,1);
for i=1:n
    %------------progress bar--------------%
    fprintf('combTar %d.\n',i);
    %--------------------------------------%
    ocr_comb=combTar_OCR{i};
    ocr_muOverlap=zeros(size(ocr_comb,1),2);  %-initiate as "0".
    for j=1:size(ocr_comb,1)
        use_CFINFEn_1=GeneID(ocr_comb{j,1});
        MU_overlap_1=intersect(use_CFINFEn_1,knownMu_inNet);
        ocr_muOverlap(j,1)=length(MU_overlap_1);
        
        use_CFINFEn_2=GeneID(ocr_comb{j,2});
        MU_overlap_2=intersect(use_CFINFEn_2,knownMu_inNet);
        ocr_muOverlap(j,2)=length(MU_overlap_2);
    end
    combTar_OCRmu3{i}=ocr_muOverlap; 
end
save cosmicValid combTar_OCRmu3 -append

