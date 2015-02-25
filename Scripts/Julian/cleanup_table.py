import sys, re

def clean_string(string):
    # remove unnecessary leading/trailing characters
    strips = [",",";",'"'," "]
    while string[0] in strips or string[-1] in strips:
        for element in strips:
            string = string.strip(element)
    # remove multiple whitespaces
    string = re.sub(' +',' ', string)
    return string

def main(filename):
    print_file = []
    with open(filename,"r") as open_file:
        for line in open_file:
            line = line.strip("\n")
            cleaned_line = [clean_string(x) for x in line.split("\t")]
            print_file.append(cleaned_line)
    dot_pos = filename.rfind(".")
    cleaned_filename = filename[:dot_pos]+"_cleaned."+filename[dot_pos+1:]
    with open(cleaned_filename,"w") as write_file:
        for line in print_file:
            write_file.write("\t".join(line)+"\n")
    
if __name__ == "__main__":
    main(sys.argv[1])