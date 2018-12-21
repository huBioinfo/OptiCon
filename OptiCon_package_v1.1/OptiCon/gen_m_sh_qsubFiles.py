import os, sys

fid = open('./NumberOfRun.txt', 'r')
NumberOfRun=fid.readlines()

outputDir = "./mFiles/"

for i in range(int(NumberOfRun[0])):
        if True:
                f = file(outputDir + "gen_relaxResults_" + str(i+1) + ".m", "w")
                f.write("clear" + "\n")
                f.write("fix(clock)" + "\n")
                f.write("initial=" + str(i+1) + ";" + "\n")
                f.write("gen_resultsFun(initial);" + "\n")
                f.write("fix(clock)" + "\n")
                f.close

bash1 = "#!/bin/bash"
bash2 = "#$ -o "
bash3 = "#$ -j y"
bash5 = "#$ -l h_vmem=16G"    #configure the memory setting
bash6 = "#$ -l mem_free=16G"  #configure the memory setting
bash4 = "#$ -cwd"
command = "/sw/matlab/el6/current/bin/matlab -nodisplay -nodesktop -nosplash -r "      #one should specify the absolute path of the executable Matlab program in your HPC.

source_forShFiles = "./mFiles/"
output_forShFiles = "./shFiles/"

DirList_forShFiles = os.listdir(source_forShFiles)
for i in range(0, len(DirList_forShFiles)):
        name_forShFiles = DirList_forShFiles[i].split(".")[0]
        if True:
                f = file(output_forShFiles + name_forShFiles + "_job.sh", "w")
                f.write(bash1 + "\n")
                f.write("\n")
                f.write(bash2 + name_forShFiles + "_log.log" + "\n")
                f.write("\n")
                f.write(bash3 + "\n")
                f.write("\n")
                f.write(bash5 + "\n")
                f.write("\n")
                f.write(bash6 + "\n")
                f.write("\n")
                f.write(bash4 + "\n")
                f.write("\n")
                f.write(command + name_forShFiles)
                f.close


source_forQSUB = "./shFiles/"

f = file("QSUB.sh", "w")
DirList_forQSUB = os.listdir(source_forQSUB)
for i in range(0, len(DirList_forQSUB)):

        name_forQSUB = DirList_forQSUB[i].split(".")[0]
        if True:
                f.write("qsub " + name_forQSUB + ".sh" + "\n")
                f.close


