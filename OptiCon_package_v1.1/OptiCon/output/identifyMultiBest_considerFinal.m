clear
D=dir('../finTherapTar_*.mat'); 
nSolution=size(D,1);
optTar=[];
optTar_MM=[];
optTar_freq=[];
for i=1:nSolution
    %------------progress bar--------------%
    fprintf('Relax %d.\n',i);
    %--------------------------------------%
    fileName=D(i).name;
    load(strcat('../',fileName)); 
    optTar=[optTar;final_tar];
    optTar_MM=[optTar_MM;final_m]; %--optTar has same order with optTar_MM.
    clearvars -except D i nSolution optTar optTar_freq optTar_MM
end
optTar_uniq=unique(optTar); 
for i=1:size(optTar_uniq,1)
    count=0;
    for j=1:size(optTar,1)
        if isequal(optTar_uniq(i),optTar(j)) 
            count=count+1;
        end
    end
    freq=count/nSolution;
    optTar_freq=[optTar_freq;freq]; %-optTar_freq has same order with optTar_uniq.
end
save MultiBest_considerFinal optTar_uniq optTar_freq optTar optTar_MM nSolution

clear
FreqCutoff=0.06;  %-This is a cutoff with False Discovery rate (FDR) of 0.05 using a null distribution.
save MultiBest_considerFinal FreqCutoff -append

%---identify recurrent OCNs.--sorted version of optTar_uniq-----%
clear
load MultiBest_considerFinal
[optTar_freq_sorted,ix]=sort(optTar_freq,'descend');
optTar_uniq_sorted=optTar_uniq(ix); %-Do so because optTar_freq has same order with optTar_uniq.
n=size(optTar_uniq_sorted,1);
recurTar_freq=[];
recurTar=[];
for i=1:n
    if optTar_freq_sorted(i) > FreqCutoff 
        recurTar_freq=[recurTar_freq;optTar_freq_sorted(i)];
        recurTar=[recurTar;optTar_uniq_sorted(i)];
    end
end
save MultiBest_considerFinal recurTar_freq recurTar -append

%--extract symbol of recurTar---%
clear
load MultiBest_considerFinal
load ../RNgeneID
recurTar_Sym=[];
n=size(recurTar,1);
for i=1:n  
     En=GeneID(recurTar(i)); 
     vector=strcmp(En,symbol2entrez_Integ(:,2));
     index=find(vector);
     Sym=symbol2entrez_Integ(index,1);
     recurTar_Sym=[recurTar_Sym;Sym];  
end
save MultiBest_considerFinal recurTar_Sym -append

%----identify corresonding MM of recurrent OCNs--------------------%
clear
load MultiBest_considerFinal
n=size(recurTar,1);
recurTar_MM=cell(n,1);
for i=1:n
    for j=1:size(optTar,1)
        if isequal(recurTar(i),optTar(j))
                recurTar_MM{i}=[recurTar_MM{i};optTar_MM(j)];
        end  %-recurTar_MM has same order with recurTar.
    end
end 
save MultiBest_considerFinal recurTar_MM -append

%----extract corresponding OCR (OCRs) of recurrent OCNs------%
clear
load MultiBest_considerFinal
n=size(recurTar,1);
recurTar_OCR=cell(n,1);
for i=1:n
    %------------progress bar--------------%
    fprintf('recurTar %d.\n',i);
    %--------------------------------------%
    ocr_mm=recurTar_MM{i};
    ocr_mm=unique(ocr_mm); 
    ocr=cell(size(ocr_mm,1),1);
    for j=1:size(ocr_mm,1) 
        str3=['load ../CFINF_hepa_' num2str(ocr_mm(j)) ];
        eval(str3);
        str4=['cfinf=cfinf_hepa_' num2str(ocr_mm(j)) ';'];
        eval(str4);
        index=find(cfinf(recurTar(i),:));
        ocr{j}=cfinf(recurTar(i),index)';
        clearvars -except ocr ocr_mm i j n recurTar_OCR recurTar_MM recurTar
    end
    recurTar_OCR{i}=ocr; %--recurTar_OCR, recurTar_MM and recurTar have same order.
end
save MultiBest_considerFinal recurTar_OCR -append


%----generate all possible combination of OCNs and corresponding OCRs--------%
clear
nk=2; 
load MultiBest_considerFinal
combTar=combnk(recurTar,nk);
nComb=size(combTar,1);
combTar_OCR=cell(nComb,1);
for i=1:nComb 
   ind1=find(recurTar==combTar(i,1)); 
   ind2=find(recurTar==combTar(i,2));
   ocr=[];
   for pp=1:size(recurTar_OCR{ind1},1)
       for qq=1:size(recurTar_OCR{ind2},1)
           ocr=[ocr;[recurTar_OCR{ind1}(pp),recurTar_OCR{ind2}(qq)]];
       end
   end
   combTar_OCR{i}=ocr; %--combTar_OCR has same order with combTar.
end
save MultiBest_considerFinal combTar combTar_OCR -append 
