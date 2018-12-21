function [] = gen_randSynScoreFun(jobNum)
%GEN_RANDSYNSCOREFUN Summary of this function goes here
%   Detailed explanation goes here

load ../RNgeneID
load cosmicValid
nComb=2; %-may change. How many nodes we randomly selected.
MMnum=1000;  %-may change. Here use 1000 SCCs.
randTimes=50000;  %-may change. 50,000 for each job. There're 200 jobs.
densityRand=ones(randTimes,1);
overlapMuRand=ones(randTimes,1);
n=length(GeneID);
correct=0;
correct2=0;
for kkk=1:randTimes
    randOCR=cell(nComb,1);  %-initiate
    %------------progress bar--------------%
    fprintf('kkk:%d.\n',kkk);
    %--------------------------------------%
    rng('shuffle');
    out1 = randperm(n);
    randNode=out1(1:nComb);
    randNode=randNode';

    rng('shuffle');
    a=rand(nComb,1); 
    b=MMnum.*a; 
    randMM=floor(b);  
    %-above will generate same Ind.It's correct for SCCs.
    
    %--extract OCR of randNode-------%
    for i=1:nComb
        str3=['load ../CFINF_hepa_' num2str(randMM(i)) ];
        eval(str3);
        str4=['cfinf=cfinf_hepa_' num2str(randMM(i)) ';'];
        eval(str4);
        index=find(cfinf(randNode(i),:));
        randOCR{i}=cfinf(randNode(i),index)';
        clearvars -except nComb n randMM randNode randOCR i kkk MMnum densityRand adjMatrix correct randTimes GeneID knownMu_inNet overlapMuRand jobNum correct2
    end
    
    %--compute cosmic overlap enrichment p-value------------------------------%
    ocr_1=randOCR{1};
    ocr_2=randOCR{2};
    ocr_1_En=GeneID(ocr_1);
    ocr1_MU_overlap=length(intersect(ocr_1_En,knownMu_inNet));
    ocr1_MU_p=1-hygecdf(ocr1_MU_overlap-1,length(GeneID),length(knownMu_inNet),length(ocr_1_En));
        
    ocr_2_En=GeneID(ocr_2);
    ocr2_MU_overlap=length(intersect(ocr_2_En,knownMu_inNet));
    ocr2_MU_p=1-hygecdf(ocr2_MU_overlap-1,length(GeneID),length(knownMu_inNet),length(ocr_2_En));
    
    if ocr1_MU_p<1.0E-20  % very few can hold.
            ocr1_MU_p=1.0E-20;
    end
    if ocr2_MU_p<1.0E-20  % very few can hold.
            ocr2_MU_p=1.0E-20;
    end
    
    overlapMuRand(kkk)=((-log10(ocr1_MU_p))*(-log10(ocr2_MU_p)))^(1/2);
    
    
    %-----compute interaction density between OCRs---------%
    ocr_1=randOCR{1};
    ocr_2=randOCR{2};
       
    ocr_overlap=intersect(ocr_1,ocr_2); %-checked.
    if ~isempty(ocr_overlap)
        overlapSubnet=comp_OCRoverlapEdge(ocr_overlap,adjMatrix);
        overlapEdgeNum=size(overlapSubnet,1); %-edges for density--will be computed twice.
    else
        overlapEdgeNum=0;
    end
        
    ocr_interact=comp_OCRinteract(ocr_1,ocr_2,adjMatrix);  %-ocr_1 and ocr_2 must not be empty.
    ocr_interactNum=size(ocr_interact,1);
    ocr_interactReal=unique(ocr_interact,'rows');
    ocr_interactRealNum=size(ocr_interactReal,1);
        
    if ocr_interactRealNum==(ocr_interactNum-overlapEdgeNum)  %if true, it means "ocr_interactReal" and Num are correct!
        correct=correct+1;
    end
    
    %-from ocr1-->ocr2---%
    part1=(length(ocr_1)-length(ocr_overlap))*length(ocr_2);
    part2=length(ocr_overlap)*(length(ocr_2)-length(ocr_overlap));
    
    %-from ocr2-->ocr1---%
    part3=(length(ocr_2)-length(ocr_overlap))*length(ocr_1);
    part4=length(ocr_overlap)*(length(ocr_1)-length(ocr_overlap));
    
    %-overlap carefully--%
    part5=length(ocr_overlap)*(length(ocr_overlap)-1); %-do not consider self-loop of overlapped genes, since there's no selfloop in large network.
    
    allPossibleNum=part1+part2+part3+part4+part5;
        
    if allPossibleNum==(length(ocr_1)*length(ocr_2)*2-(length(ocr_overlap)*length(ocr_overlap))-length(ocr_overlap)) %-the last item is #of selfloop.
        correct2=correct2+1;
    end
    
    densityRand(kkk)=ocr_interactRealNum/allPossibleNum;
    
end
%--change the name for consistency.
s1Rand=overlapMuRand;
s2Rand=densityRand;

str9=['save randSynScore' num2str(jobNum) ' s1Rand s2Rand correct correct2'];
eval(str9);

end

