%-----------Import gene expression data-----------------%
clear
%------------progress bar--------------%
fprintf('Importing gene expression data.\n');
%--------------------------------------%
load RNgeneID
FpkmCorrect_dm=bioma.data.DataMatrix('File', 'GeneExpression.txt');
[RNexp_hepa,count] = gen_RNexp(GeneID,FpkmCorrect_dm);
save RNexp RNexp_hepa

%----Import deregulation score (DScore)-----------------%
clear
%------------progress bar--------------%
fprintf('Importing deregulation score (DScore).\n');
%--------------------------------------%
fid4=fopen('DScore.txt');
microsam_hepa=textscan(fid4,'%s %f','delimiter','\t');  
fclose(fid4);
microsamID_hepa=microsam_hepa{1};
microsamSco_hepa=microsam_hepa{2};
save micro_hepa microsamID_hepa microsamSco_hepa

clear
load micro_hepa
load RNgeneID
[RNsamSco_hepa,count] = RNsamScore(GeneID,microsamID_hepa,microsamSco_hepa);
save sam_hepa GeneID RNsamSco_hepa

clear
load sam_hepa
RN_deregSum_hepa=0;  
n=size(RNsamSco_hepa,1);
for i=1:n
        RN_deregSum_hepa=RN_deregSum_hepa+RNsamSco_hepa(i);
end

RN_WeakInd_hepa=[];  
for j=1:n
    if RNsamSco_hepa(j)==0  
        RN_WeakInd_hepa=[RN_WeakInd_hepa;j];
    end
end
save sam_hepa RN_deregSum_hepa RN_WeakInd_hepa -append 

%--------------Compute indirect control value (ICV)--------------%
clear
%------------progress bar--------------%
fprintf('Computing indirect control value (ICV).\n');
%--------------------------------------%
fix(clock)
load sam_hepa
load RNgeneID
load RNexp

[ICV_valueTable,ICV_SPtable] = gen_cis(RNsamSco_hepa,adjMatrix,RNexp_hepa);
save cis_hepa ICV_valueTable  
save ciscis ICV_SPtable -v7.3 
fix(clock)

%-----------------Generate indirect control region----------------------------%
clear
load cis_hepa ICV_valueTable
load ciscis 
n=size(ICV_valueTable,1); 
RNcis03=cell(n,1);  
for i=1:n
    %------------progress bar--------------%
    %fprintf('Gene: %d.\n',i);
    %--------------------------------------%
    index3=[];
    if ~isempty(ICV_valueTable{i}) 
        Imp=ICV_valueTable{i};
        SP=ICV_SPtable{i}; 
        
        nimp=length(Imp); 
        for k=1:nimp
            if Imp(k)>=0.3 
               index3=[index3,SP{k}]; 
            end
        end

        index3=unique(index3);
        RNcis03{i}=index3;
    end
    
end
save cis_hepa RNcis03 -append

%---------------------Generate control region (CR) of each gene---------------%
clear
%------------progress bar--------------%
fprintf('Generating control region.\n');
%--------------------------------------%
fix(clock)
nrun=1000;
for i=0:(nrun-1)
    %------------progress bar--------------%
    fprintf('SCC %d.\n',i);
    %--------------------------------------%
    load sam_hepa
    load cis_hepa RNcis03
    str6=['load CF_' num2str(i) ];
    eval(str6);
 
    str7=['cfinf_hepa_' num2str(i) '=generate_cfinf(cf_' num2str(i) ',RNsamSco_hepa,RNcis03);'];
    eval(str7);
 
    str8=['cfinf_deregSco_hepa_' num2str(i) '=gen_CFINFderegSco(cfinf_hepa_' num2str(i) ',RNsamSco_hepa,RN_deregSum_hepa);'];
    eval(str8);
    %--save
    str9=['save CFINF_hepa_' num2str(i) ' cfinf_hepa_' num2str(i) ' cfinf_deregSco_hepa_' num2str(i) ];
    eval(str9);

    clearvars -except nrun i
end
fix(clock)
%---------Compute initial points for greedy search-------------%
clear
%------------progress bar--------------%
fprintf('Computing initial points for greedy search.\n');
%--------------------------------------%
gen_initialPoints

%------------progress bar--------------%
fprintf('\n\nStep 1 finished. You can start Step 2.\n');
%--------------------------------------%



