import os
from IonPickerEngine import deisotope

fileslist = []

for file in os.listdir("."):
    if file.endswith(".csv") and "output" not in file and 'ions' not in file:
        fileslist.append(file)
        print "processing", file

for file in fileslist:
    output_fname = 'output' + file
    print "output file = ",output_fname

    new_compound_list = deisotope(file, output_file = output_fname,max_error=10, check_gaps=True)
    
