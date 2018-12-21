clear
fix(clock)
load RNgeneID
load sam_hepa
ngenes=size(GeneID,1);
initialOptiRat=zeros(ngenes,1000);
initialDerRat=zeros(ngenes,1000);
initialWeakRat=zeros(ngenes,1000);

topFracInd=ceil(ngenes*1000*0.0001);  %top (p<0.01%) as initial points.

for j=0:999
    %------------progress bar--------------%
    %fprintf('CF %d.\n',j);
    %--------------------------------------%
    str3=['load CFINF_hepa_' num2str(j) ];
    eval(str3);
    str4=['cfinf=cfinf_hepa_' num2str(j) ';'];
    eval(str4);
    str5=['cfinf_deregSco=cfinf_deregSco_hepa_' num2str(j) ';'];
    eval(str5);
    
    initialDerRat(:,(j+1))=cfinf_deregSco;
    
    for i=1:ngenes
        index=find(cfinf(i,:));  % index must have at least one value,i.e. itself.
        tempCFINF=cfinf(i,index);
        tempCFINF=tempCFINF';
        
        tempWeakInd=intersect(tempCFINF,RN_WeakInd_hepa);
        initialWeakRat(i,(j+1))=length(tempWeakInd)/length(RN_WeakInd_hepa);
        
    end
    clearvars -except j i initialDerRat initialWeakRat topFracInd ngenes RN_WeakInd_hepa initialOptiRat
end

for p=1:ngenes
    for q=1:1000
        initialOptiRat(p,q)=initialDerRat(p,q)-initialWeakRat(p,q);
    end
end

[bb,ix]=sort(initialOptiRat(:),'descend');
top=size(find(bb>bb(topFracInd)),1);
NodeMMrate=zeros(top,3);
[rowInd,colInd]=ind2sub(size(initialOptiRat),ix(1:top));

for kkk=1:top
    NodeMMrate(kkk,:)=[rowInd(kkk),(colInd(kkk)-1),initialOptiRat(rowInd(kkk),colInd(kkk))];
end
save InitialPoints initialOptiRat NodeMMrate
fix(clock)



