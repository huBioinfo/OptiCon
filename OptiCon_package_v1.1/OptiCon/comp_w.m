function w = comp_w(v,aa,dys_g,RNexp_ph)
%COMP_W


e=0.18;  %-This is the Epsilon for lung cancer data. This should be changed based on the number of samples.

gc=RNexp_ph(dys_g,:); 
ga=RNexp_ph(v,:);
gb=RNexp_ph(aa,:);

if isempty(find(gc))
    Pear1=0;
    Pear2=0;
else  
 
    if isempty(find(ga))
        Pear1=0;
    else
        R1=corrcoef(gc',ga');
        Pear1=abs(R1(1,2));
    end
    if isempty(find(gb))
        Pear2=0;
    else
        R2=corrcoef(gc',gb');
        Pear2=abs(R2(1,2));
    end
end

w=(Pear1+Pear2)/2;
if w<e
    w=e;
end

end

