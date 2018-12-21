import os, sys

outputDir = "./mFiles_randNull/"

for i in range(1, 201):
        if True:
                f = file(outputDir + "gen_randSynScore_" + str(i) + ".m", "w")
                f.write("clear" + "\n")
                f.write("fix(clock)" + "\n")
                f.write("jobNum=" + str(i) + ";" + "\n")
                f.write("gen_randSynScoreFun(jobNum);" + "\n")
                f.write("fix(clock)" + "\n")
                f.close


bash1 = "#!/bin/bash"
bash2 = "#$ -o "
bash3 = "#$ -j y"
bash5 = "#$ -l h_vmem=16G"    #configure the memory setting
bash6 = "#$ -l mem_free=16G"  #configure the memory setting
bash4 = "#$ -cwd"
command = "/sw/matlab/el6/current/bin/matlab -nodisplay -nodesktop -nosplash -r "      #one should specify the absolute path of the exeutable Matlab program in your HPC.

source_forShFiles = "./mFiles_randNull/"
output_forShFiles = "./shFiles_randNull/"

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


source_forQSUB = "./shFiles_randNull/"

f = file("QSUB_randNull.sh", "w")
DirList_forQSUB = os.listdir(source_forQSUB)
for i in range(0, len(DirList_forQSUB)):

        name_forQSUB = DirList_forQSUB[i].split(".")[0]
        if True:
                f.write("qsub " + name_forQSUB + ".sh" + "\n")
                f.close


