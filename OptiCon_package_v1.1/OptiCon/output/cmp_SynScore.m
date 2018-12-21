clear
load cosmicValid
load ../RNgeneID
load MultiBest_considerFinal
s1=cell(size(combTar_OCRmu3,1),1); %--combTar_OCRmu3 and combTar_OCR have same order.
for i=1:size(combTar_OCRmu3,1)
    %------------progress bar--------------%
    fprintf('combTar_OCRmu3 %d.\n',i);
    %--------------------------------------%
    ocr_muOverlap=combTar_OCRmu3{i};
    ocr_comb=combTar_OCR{i};
    s1{i}=zeros(size(ocr_muOverlap,1),1);  %-initiate as "0".
    for j=1:size(ocr_muOverlap,1)
        p1=1-hygecdf(ocr_muOverlap(j,1)-1,length(GeneID),length(knownMu_inNet),length(ocr_comb{j,1}));
        p2=1-hygecdf(ocr_muOverlap(j,2)-1,length(GeneID),length(knownMu_inNet),length(ocr_comb{j,2}));
        %-The max of p1 or p2 is 1 which the overlap is zero. The min is set to
        %1.0E-20
        if p1<1.0E-20 
            p1=1.0E-20;
        end
        if p2<1.0E-20 
            p2=1.0E-20;
        end
        
        s1{i}(j)=((-log10(p1))*(-log10(p2)))^(1/2);
    end
end
save SynergyScore s1

clear
load interactDensity
s2=combTar_OCRinteractDensity; 
save SynergyScore s2 -append
