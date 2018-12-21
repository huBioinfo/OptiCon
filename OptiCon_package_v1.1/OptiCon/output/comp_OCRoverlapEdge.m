function subnet = comp_OCRoverlapEdge(coninfed,adjMatrix)
%TARGETS_CIRNET
%   Detailed explanation goes here
subnet=[];
n=size(coninfed,1);
for i=1:n 
    for j=1:n
        if adjMatrix(coninfed(i),coninfed(j))==1
            source=coninfed(i);
      
            terminal=coninfed(j);
   
            subnet=[subnet;[source,terminal]];
        end
    end
end

end

