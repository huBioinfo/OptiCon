clear
fix(clock)
load MultiBest_considerFinal
load ../RNgeneID
n=size(combTar_OCR,1);
%combTar_OCRoverlap=cell(n,1); 
combTar_OCRoverlapSize=cell(n,1);
%combTar_OCRinteract=cell(n,1);
combTar_OCRinteractDensity=cell(n,1);
correctFlag=0;
correctFlag2=0;
for i=1:n
    %------------progress bar--------------%
    fprintf('combTar %d.\n',i);
    %--------------------------------------%
    ocr_comb=combTar_OCR{i};
    ocr_overlap=cell(size(ocr_comb,1),1);
    ocr_overlapSize=zeros(size(ocr_comb,1),1);
    %ocr_interactReal=cell(size(ocr_comb,1),1);
    ocr_interactDensity=zeros(size(ocr_comb,1),1);
    for j=1:size(ocr_comb,1)
        %------------progress bar--------------%
        fprintf('ocr_comb %d.\n',j);
        %--------------------------------------%
        ocr_1=ocr_comb{j,1};
        ocr_2=ocr_comb{j,2};
        
        ocr_overlap{j}=intersect(ocr_1,ocr_2); 
        ocr_overlapSize(j)=length(ocr_overlap{j});
        if ocr_overlapSize(j)~=0
            overlapSubnet=comp_OCRoverlapEdge(ocr_overlap{j},adjMatrix);
            overlapEdgeNum=size(overlapSubnet,1);
        else
            overlapEdgeNum=0;
        end
        
        ocr_interact=comp_OCRinteract(ocr_1,ocr_2,adjMatrix);
        ocr_interactNum=size(ocr_interact,1);
        ocr_interactReal=unique(ocr_interact,'rows');
        ocr_interactRealNum=size(ocr_interactReal,1);
        
        if ocr_interactRealNum~=(ocr_interactNum-overlapEdgeNum) 
            correctFlag=1;  %if never goes to 1, then correct!
        end
        
        %-from ocr1-->ocr2---%
        part1=(length(ocr_1)-ocr_overlapSize(j))*length(ocr_2);
        part2=ocr_overlapSize(j)*(length(ocr_2)-ocr_overlapSize(j));
 
        %-from ocr2-->ocr1---%
        part3=(length(ocr_2)-ocr_overlapSize(j))*length(ocr_1);
        part4=ocr_overlapSize(j)*(length(ocr_1)-ocr_overlapSize(j));
        
        %-overlap carefully--%
        part5=ocr_overlapSize(j)*(ocr_overlapSize(j)-1);
       
        allPossibleNum=part1+part2+part3+part4+part5;
        
        if allPossibleNum~=(length(ocr_1)*length(ocr_2)*2-(ocr_overlapSize(j)*ocr_overlapSize(j))-ocr_overlapSize(j))
            correctFlag2=1;  %if never goes to 1, then correct!
        end
        
        ocr_interactDensity(j)=ocr_interactRealNum/allPossibleNum;
    end
    %combTar_OCRoverlap{i}=ocr_overlap;
    combTar_OCRoverlapSize{i}=ocr_overlapSize;
    %combTar_OCRinteract{i}=ocr_interactReal;
    combTar_OCRinteractDensity{i}=ocr_interactDensity;
end
save interactDensity combTar_OCRoverlapSize combTar_OCRinteractDensity correctFlag correctFlag2
fix(clock)
