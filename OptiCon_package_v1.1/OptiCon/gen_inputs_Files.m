clear
load InitialPoints
NumberOfRun=size(NodeMMrate,1);
fid=fopen('NumberOfRun.txt','w');
fprintf(fid,'%d',NumberOfRun);
fclose(fid);

