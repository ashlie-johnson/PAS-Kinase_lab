"""
This program takes as input 2 files which contain sequences for 2 different proteins of the same species (for example,
yeast UGP1 and yeast PAS-kinase).

To choose the input files to be read:
    Input files must be FAASTA-aligned .txt files (I used this website to align sequences: https://www.ebi.ac.uk/Tools/msa/muscle/
    Save the .txt files to be read inside the same folder as this python project
    Select "Edit Configurations"
    Inside the "Parameters" box, enter the exact titles of the files, including .txt, separated by a space.
    Run the program

The user is given 2 options at the start of the program:
1. check for missing species--this will create an output file titled "missing_species.txt", which is a list of every
species name that is listed in one file but not the other.
2. create an output file with pairs of all possible combinations of isoforms--this will create
an output file titled "paired_sequences.txt", which contains all possible pairings of sequences for direct coupling
analysis. (for example, some_species isoform 1 from file 1 will be paired with some_species isoforms 2, 3, and 4 from
file 2, etc)

To select an option, type either "1" or "2" in the terminal, then press enter. Rerun the program to select a different
option.

*****Please note that if either of these options is selected more than once, the file names remain the same and the
previous data will be overwritten with the new data.
(A potential future modification to the code would be to let the user create a new file name every time they run the
program)

Created by Ashlie Johnson, ashlie.johnson07@gmail.com
"""


import sys

file_name_1 = str(sys.argv[1])
file_name_2 = str(sys.argv[2])

matches = []

def main():
    with open(file_name_1) as file1:
        with open(file_name_2) as file2:

            lines1 = file1.readlines()
            lines2 = file2.readlines()

            f1sequenceDict = fillDict(lines1)
            f2sequenceDict = fillDict(lines2)

    choice = input(
        "Would you like to 1. check for missing species in files or 2. create an output file with pairs of all possible combinations of isoforms?")

    if choice == "1":
        findMissing(True, lines1, lines2, f1sequenceDict)

    elif choice == "2":
        findMissing(False, lines1, lines2, f1sequenceDict)
        pairSequences(f1sequenceDict, f2sequenceDict)

    else:
        print("Invalid entry")




def findMissing(showMissing, lines1, lines2, f1sequenceDict):   #creates new output file with missing species names
    with open("missing_species.txt", "w") as out_file:
            for line in lines2:
                if ">" in line:
                    species = removeFormatting(line)
                    f2newSpecies = removeIsoform(species)
                    f2newSpecies.sort()

                    for file1species in f1sequenceDict.keys():
                        f1wordList = file1species.split(" ")

                        if f1wordList == f2newSpecies:
                            matches.append(f1wordList)


            if showMissing:
                out_file.write("\n")
                out_file.write("                    File 1 species missing from file 2:")
                out_file.write("\n")

            for line in lines1:
                if ">" in line:
                    species = removeFormatting(line)
                    f1species = removeIsoform(species)
                    f1species.sort()

                    if f1species not in matches:
                        if showMissing:
                            displayedSpecies= line.replace(">", "")
                            out_file.write(displayedSpecies)

            if showMissing:
                out_file.write("\n\n\n")
                out_file.write("                    File 2 species missing from file 1:")
                out_file.write("\n")


            for line in lines2:
                if ">" in line:
                    species = removeFormatting(line)
                    f2species = removeIsoform(species)
                    f2species.sort()

                    if f2species not in matches:
                        if showMissing:
                            displayedSpecies= line.replace(">", "")

                            out_file.write(displayedSpecies)



def pairSequences(f1sequenceDict, f2sequenceDict):      #creates new output file with sequences paired
    with open("paired_sequences.txt", "w") as out_file:
        allSpecies = f1sequenceDict.keys()

        for species1 in allSpecies:
            for species2 in allSpecies:
                if species1.split(" ") in matches and species2.split(" ") in matches:
                    newSpecies1 = removeNumbers(species1)
                    newSpecies2 = removeNumbers(species2)
                    if newSpecies1 == newSpecies2:
                        out_file.write(">")
                        out_file.write(species1)
                        out_file.write("\n")
                        out_file.write(f1sequenceDict[species1])
                        out_file.write("\n")

                        out_file.write(">")
                        out_file.write(species2)
                        out_file.write("\n")
                        out_file.write(f2sequenceDict[species2])
                        out_file.write("\n")


def removeFormatting(line):
    line = line.lower()
    line = line.replace("isoform", "")
    line = line.replace(">", "")
    line = line.replace("\n", "")
    line = line.replace("(", "")
    line = line.replace(")", "")
    line = line.replace("_", " ")
    line = line.replace("  ", " ")
    line = line.replace(".", "")
    line = line.replace('\u2028', "")
    words = line.split(" ")

    return words


def removeIsoforms(words):
    for word in words:
        if word.isnumeric:
            words.remove(word)
    return words


def removeIsoform(words):
    newSpecies = []
    for word in words:
        if any(str.isdigit(c) for c in word):
            if "x" in word or "X" in word:
                word = word.replace("x", "")
                word = word.replace("X", "")
                newSpecies.append(word)
        elif word != "isomer":
            newSpecies.append(word)
    return newSpecies


def removeNumbers(species):
    newSpecies = []
    wordList = species.split(" ")
    for word in wordList:
        if word.isalpha():
            newSpecies.append(word)
    return newSpecies


def fillDict(lines):
    sequenceDict = {}
    speciesString = ""
    sequence = ""

    for line in lines:
        if ">" in line:
            if sequence != "":
                sequenceDict[speciesString] = sequence
                sequence = ""
            species = removeFormatting(line)
            newSpecies = removeIsoform(species)
            newSpecies.sort()
            speciesString = ' '.join(newSpecies)
            sequenceDict[speciesString] = ""
        else:
            sequence += line

    return sequenceDict


main()
