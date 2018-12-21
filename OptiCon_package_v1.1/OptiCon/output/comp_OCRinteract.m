function interact = comp_OCRinteract(ocr_1,ocr_2,adjMatrix)
%comp_OCRinteract
%   Detailed explanation goes here
interact=[];
n_1=length(ocr_1);
n_2=length(ocr_2);
for i=1:n_1 
    
    for j=1:n_2 
        if adjMatrix(ocr_1(i),ocr_2(j))==1 
            source=ocr_1(i);
          
            terminal=ocr_2(j);
            
            interact=[interact;[source,terminal]];
        end
    end
end

i=0;
j=0;

for i=1:n_2 
    
    for j=1:n_1
        if adjMatrix(ocr_2(i),ocr_1(j))==1 
            source=ocr_2(i);
         
            terminal=ocr_1(j);
       
            interact=[interact;[source,terminal]];
        end
    end
end

end

