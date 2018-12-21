%-----------Import SCCs to Matlab variables---------------%
clear
fix(clock)
nrun=1000;
for i=0:(nrun-1)
    %------------progress bar--------------%
    fprintf('Importing SCC %d.\n',i);
    %--------------------------------------%
    tempCRfile=sprintf('tempCRMatrix_%d.txt',i);
    tempCR=load(tempCRfile);
    str1=['cf_' num2str(i) '=generate_CF(tempCR);'];
    eval(str1);
    str2=['save CF_' num2str(i) ' cf_' num2str(i) ];
    eval(str2);
    clearvars -except nrun i
end
fix(clock)

%---Import GeneID in the network--%  
clear
fid1=fopen('GeneIDindex.txt');
GeneID=textscan(fid1,'%s');
GeneID=GeneID{1};
fclose(fid1);
save RNgeneID GeneID

%---Generate adjacent matrix of the network-----%
clear
load RNgeneID
fid2=fopen('MyGeneNetwork.txt');
integ=textscan(fid2,'%s %s','delimiter','\t');
fclose(fid2);
a=integ{1};
b=integ{2};
integ_RN=[a,b];
adjMatrix = gen_adjM(GeneID,integ_RN);
save RNgeneID integ_RN adjMatrix -append

clear
fid3=fopen('symbol2entrez_Integ.txt');
symbolEntrez=textscan(fid3,'%s %s','delimiter','\t');
fclose(fid3);
c=symbolEntrez{1};
d=symbolEntrez{2};
symbol2entrez_Integ=[c,d];
save RNgeneID symbol2entrez_Integ -append

