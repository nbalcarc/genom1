import os
import shutil
import tempfile
from zipfile import ZipFile
import numpy as np


def preprocess(file_dir: str) -> str:
    '''Removes all irrelevant information from the given genome file'''
    
    with open(file_dir, 'r') as file:
        with open(file_dir + ".copy", 'w') as file1:
            for line in file.readlines():
                if not line.startswith('>'): #all irrelevant lines start with >
                    file1.write(line)
    return file_dir + ".copy"

def preprocess1(file_dir: str) -> str:
    '''Removes all irrelevant information from the given genome file and compresses the data'''
    
    accumulator: str = ""
    all_lines: str = ""
    
    # retrieve the file data
    with open(file_dir, 'r') as file:
        all_lines = file.readlines()
    
    # start preprocessing before compression
    for line in all_lines:
        if not line.startswith('>'): #all irrelevant lines start with >
            accumulator += line.strip() #push the line to the accumulator
        else:
            accumulator += ">" #compress the irrelevant line into a unique character
            accumulator += "_" * (0 if len(accumulator) % 3 == 0 else 3 - len(accumulator) % 3) #pad until divisible by 3
            
    #accumulator = accumulator[1:].replace("\n", "") #remove all newlines and the first >
    accumulator = accumulator[3:]
    accumulator += "_" * (0 if len(accumulator) % 3 == 0 else 3 - len(accumulator) % 3) #have the file be a length divisible by 3, pad with _
    
    # start compression
    #compressed = np.empty(int(len(accumulator) / 3), np.uint8)
    compressed = np.zeros(int(len(accumulator) / 3), np.uint8)
    
    for i in range(int(len(accumulator) / 3)):
        substr = accumulator[i*3:i*3+3]
        #print(len(substr))
        print(i)
        match substr[0]: #don't need to cover > case because it's all 0's
            case "_": #only padding
                compressed[i] += 1
                continue
            case "A":
                compressed[i] += 64 #increment 64 for at least one char
            case "C":
                compressed[i] += 64
                compressed[i] += 16 #this is a C
            case "G":
                compressed[i] += 64
                compressed[i] += 32 #this is a G
            case "T":
                compressed[i] += 64 
                compressed[i] += 48 #this is a T/U
            case "U":
                compressed[i] += 64 
                compressed[i] += 48 #this is a T/U
        
        match substr[1]: #don't need to cover < case
            case "_":
                compressed[i] += 1
                continue
            case "A":
                compressed[i] += 64
            case "C":
                compressed[i] += 64
                compressed[i] += 4
            case "G":
                compressed[i] += 64
                compressed[i] += 8
            case "T":
                compressed[i] += 64
                compressed[i] += 12
            case "U":
                compressed[i] += 64
                compressed[i] += 12
                
        match substr[2]:
            case "_":
                compressed[i] += 1
                continue
            case "A":
                compressed[i] += 64
            case "C":
                compressed[i] += 64
                compressed[i] += 1
            case "G":
                compressed[i] += 64
                compressed[i] += 2
            case "T":
                compressed[i] += 64
                compressed[i] += 3
            case "U":
                compressed[i] += 64
                compressed[i] += 3
                
    #compressed = compressed.tobytes() #convert from uint8's to bytes
    compressed.tofile(file_dir + ".copy")       
    
    # write to the output file
    #with open(file_dir + ".copy", "w") as file:
    #    file.write(accumulator)
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
        if not os.path.isdir(genomes_dir + "/" + organism_name): #case 1: no genome with the given name exists yet
            os.makedirs(genomes_dir + "/" + organism_name)
            folder_loc = genomes_dir + "/" + organism_name
        else:
            i = 1
            while os.path.isdir(genomes_dir + "/" + organism_name + "_" + str(i)): #case 2: genome already exists, come up with a unique name
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
