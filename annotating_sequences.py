# Importations
#from subprocess import call
import os

def prokka (listing, input_path, output_path, fichier) :
    for elements in listing :
        os.system(" prokka --force " + output_path + "/" + elements + "/Spades/contigs.fasta "
                  + " --prefix " + elements
                  + " --outdir " + output_path + "/" + elements + "/Prokka")
        fichier.write("\n       " + elements + " has been annotated")

def prokka_unicycler (listing, input_path, output_path, fichier) :
    for elements in listing :
        os.system(" prokka --force " + output_path + "/" + elements + "/Unicycler/assembly.fasta "
                  + " --prefix " + elements
                  + " --outdir " + output_path + "/" + elements + "/Prokka")
        fichier.write("\n       " + elements + " has been annotated")

def prokka_reference (listing, input_path, output_path, fichier) :
    for element in listing[0] :
        os.system(" prokka --force " + input_path + "/" + element + ".fasta "
                  + " --prefix " + element
                  + " --outdir " + output_path + "/" + element + "/Prokka")
        fichier.write("\n       " + element + " has been annotated")
    for element in listing[1] :
        if len(listing[1]) != 0 :
            os.system(" prokka --force " + input_path + "/" + element + ".fas "
                      + " --prefix " + element
                      + " --outdir " + output_path + "/" + element + "/Prokka")
            fichier.write("\n       " + element + " has been annotated")
    for element in listing[2] :
        if len(listing[1]) != 0:
            os.system(" prokka --force " + input_path + "/" + element + ".fa "
                      + " --prefix " + element
                      + " --outdir " + output_path + "/" + element + "/Prokka")
            fichier.write("\n       " + element + " has been annotated")
    return



# prokka : calls prokka
# --force : rewrite on a existing file
# --prefix : name of the output files
# --kingdom : transparent
# --outdir : output directory

