import os
import shutil
import tempfile
from zipfile import ZipFile


def preprocess(file_dir: str) -> str:
    '''Removes all irrelevant information from the given genome file'''
    
    with open(file_dir, 'r') as file:
        with open(file_dir + ".copy", 'w') as file1:
            for line in file.readlines():
                if not line.startswith('>'): #all irrelevant lines start with >
                    file1.write(line)
    return file_dir + ".copy"


def main():
    '''Entry point for program'''
    
    basedir = os.getcwd() + "/"
    genomes_dir = basedir + "genomes"
    genomes_raw_dir = basedir + "genomes_raw"
   
    # notify user if no genomes_raw folder
    if not os.path.isdir(genomes_raw_dir):
        print("ERROR: The genomes_raw folder does not exist!")
        return
    
    # reset the genomes folder
    if os.path.isdir(genomes_dir):
        shutil.rmtree(genomes_dir)
     
    os.makedirs(basedir + "genomes") #create the genomes folder
    zipdir = tempfile.mkdtemp() #create a temporary folder to unzip our genomes
    
    # iterate through all files in the folder
    for file in os.listdir(genomes_raw_dir):
        file_clipped = file[:-4]
        
        # unzip the folder
        with ZipFile(genomes_raw_dir + "/" + file, 'r') as zipfile:
            zipfile.extractall(zipdir + "/" + file) #ensure no two unzipped folders overlap
        genome_location = zipdir + f"/{file}/ncbi_dataset/data/{file_clipped}/"
        
        # retrieve the organism name
        organism_name = "default"
        with open(zipdir + f"/{file}/ncbi_dataset/data/assembly_data_report.jsonl") as file: #open the file containing the organism name
            for line in file.readlines():
                if "organismName" in line:
                    substr = line[line.find("organismName") + 15:]
                    organism_name = substr[:substr.find('"')].replace(" ", "_").replace("/", "_") #extract and perform input validation
        
        
        # iterate through all files in the extracted location (the file ends with .fna but the name could be anything)
        folder_loc = ""
        if not os.path.isdir(genomes_dir + "/" + organism_name):
            os.makedirs(genomes_dir + "/" + organism_name)
            folder_loc = genomes_dir + "/" + organism_name
        else:
            i = 1
            while os.path.isdir(genomes_dir + "/" + organism_name + "_" + str(i)):
                i += 1
            os.makedirs(genomes_dir + "/" + organism_name + "_" + str(i))
            folder_loc = genomes_dir + "/" + organism_name + "_" + str(i)
        for tfile in os.listdir(genome_location):
            
            # if this isn't the file we're looking for, skip
            if tfile[tfile.rfind('.'):] != ".fna":
                continue
            
            genome_file = preprocess(genome_location + tfile) #preprocess the file before copying
            shutil.copyfile(genome_file, folder_loc + "/" + tfile) #copy the file into the output
    

if __name__ == "__main__":
    main()
