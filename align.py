import sys

def readInput(file1):
    f = open(file1, "r")
    name = f.readline()[1:-1]
    seq = ""
    for line in f:
        if (not "[" in line) and (not "]" in line):
            seq = seq + line[:-1]
        else:
            name = name + line[:-1]
    return seq, name

# read scoring matrix of either protein or dna
def MatrixRead(file):
    f = open(file, "r")
    f.readline()
    f.readline()
    matrix = []
    initials = f.readline().split()
    for line in f:
        row = line.split()
        if len(row) == 1:
            gapPenalty = row[0]
        else:
            row.pop(0)
            matrix.append(row)
    return initials, matrix, gapPenalty

# create table
def tableConstruction(length1, length2):
    table = []
    for i in range(length1):
        col = []
        for j in range(length2):
            col.append([0,""])
        table.append(col)

    # add origin arrows for first row and col
    for i in range(1,len(seq1)+1):
        table[i][0][1] = "T"
    for j in range(1,len(seq2)+1):
        table[0][j][1] = "L"

    return table

# global sequence alignment function
def globalSeqAlignment(Matrix, seq1, seq2):
    initials = Matrix[0]
    matrix = Matrix[1]
    gapPen = int(Matrix[2])

    table = tableConstruction(len(seq1)+1, len(seq2)+1)

    # add gap penalty in blobal sequence alignment
    for i in range(1,len(seq1)+1):
        table[i][0][0] = table[i-1][0][0] + gapPen
    for j in range(1,len(seq2)+1):
        table[0][j][0] = table[0][j-1][0] + gapPen
    
    #Filling out table, seq1 is verticle, seq2 is horizontal
    for i in range(1,len(seq1)+1):
        for j in range(1,len(seq2)+1):
            # add match/mismatch score
            matching = table[i-1][j-1][0] + int(matrix[initials.index(seq1[i-1])][initials.index(seq2[j-1])])
            # add gap penalty to top
            topGap = table[i-1][j][0] + gapPen
            # add gap penalty to bot
            botGap = table[i][j-1][0] + gapPen

            # add the max score and origin arrow to current cell
            value = max(matching, topGap, botGap)
            origin = ""
            if value == matching:
                origin = "TL"
            elif value == topGap:
                origin = "T"
            elif value == botGap:
                origin = "L"
            table[i][j] = [value, origin]
    # for row in table:
    #     print(row)

    seq1Char, seq2Char = Backtracking(table, seq1, seq2)
    return seq1Char, seq2Char, 1, table[len(seq1)][len(seq2)][0]

# function to backtrack based on the table and 2 sequences
def Backtracking(table, seq1, seq2):
    aligning = True
    i = len(seq1)
    j = len(seq2)
    seq1Char = ""
    seq2Char = ""

    # while cell is no the first one on the table, trace back through i and/or j
    while aligning:
        curr = table[i][j]
        if curr[1] == "":
            aligning = False
        elif curr[1] == "TL":
            seq1Char = seq1[i-1] + seq1Char
            seq2Char = seq2[j-1] + seq2Char
            i = i-1
            j = j-1
        elif curr[1] == "T":
            seq1Char = seq1[i-1] + seq1Char
            seq2Char = "_" + seq2Char
            i = i-1
        elif curr[1] == "L":
            seq2Char = seq2[j-1] + seq2Char
            seq1Char = "_" + seq1Char
            j = j-1
    return seq1Char, seq2Char

# function to get the sequence without gaps from alignment, used to find the small sequences in the original sequences
def getSeqNoGap(alignedSeq1, alignedSeq2):
    seq1NoGap = ""
    seq2NoGap = ""
    for i in range(len(alignedSeq1)):
        if alignedSeq1[i] != "_":
            seq1NoGap = seq1NoGap + alignedSeq1[i]
        if alignedSeq2[i] != "_":
            seq2NoGap = seq2NoGap + alignedSeq2[i]
    return seq1NoGap, seq2NoGap

# function to print general output, used for global and semi-global sequence alignment
def printOut(alignedSeq1, alignedSeq2, seq1, seq2, score, outFile, name1, name2):
    f = open(outFile, "w")
    identities = 0

    for i in range(len(alignedSeq1)):
        if alignedSeq1[i] == alignedSeq2[i]:
            identities += 1
    seq1NoGap, seq2NoGap = getSeqNoGap(alignedSeq1, alignedSeq2)
    start1 = seq1.index(seq1NoGap) + 1
    start2 = seq2.index(seq2NoGap) + 1

    partialSeq1 = alignedSeq1
    partialSeq2 = alignedSeq2
    while len(partialSeq1) > 60:
        partialSeq1NoGap, partialSeq2NoGap = getSeqNoGap(partialSeq1[0:60], partialSeq2[0:60])
        end1 = len(partialSeq1NoGap) + start1 - 1
        end2 = len(partialSeq2NoGap) + start2 - 1
        out1 = name1 + ":  " + str(start1) + " " + partialSeq1[0:60] + " " + str(end1)
        out2 = name2 + ":  " + str(start2) + " "  + partialSeq2[0:60] + " " + str(end2)
        
        f.write(out1 + "\n")
        f.write(out2 + "\n")
        f.write("\n")
        
        start1 = len(partialSeq1NoGap) + start1
        start2 = len(partialSeq2NoGap) + start2
        partialSeq1 = partialSeq1[60:]
        partialSeq2 = partialSeq2[60:]

    out1 = name1 + ":  " + str(start1) + " " + partialSeq1 + " " + str(len(partialSeq1) + start1 - 1)
    out2 = name2 + ":  " + str(start2) + " "  + partialSeq2 + " " + str(len(partialSeq2) + start2 - 1)
    f.write(out1 + "\n")
    f.write(out2 + "\n")
    f.write("\n")

    f.write("Score: " + str(score) + "\n")
    f.write("Identities: " + str(identities) + "/" + str(len(alignedSeq1)) + " (" + str(round((identities/len(alignedSeq1))*100)) + "%)" + "\n")
    f.close()

# function to print only local alignment output, return one instance of alignment file output as string
def printOutLocal(alignedSeq1, alignedSeq2, seq1, seq2, score, name1, name2):
    f = ""
    identities = 0

    for i in range(len(alignedSeq1)):
        if alignedSeq1[i] == alignedSeq2[i]:
            identities += 1
    seq1NoGap, seq2NoGap = getSeqNoGap(alignedSeq1, alignedSeq2)
    start1 = seq1.index(seq1NoGap) + 1
    start2 = seq2.index(seq2NoGap) + 1

    partialSeq1 = alignedSeq1
    partialSeq2 = alignedSeq2
    while len(partialSeq1) > 60:
        partialSeq1NoGap, partialSeq2NoGap = getSeqNoGap(partialSeq1[0:60], partialSeq2[0:60])
        end1 = len(partialSeq1NoGap) + start1 - 1
        end2 = len(partialSeq2NoGap) + start2 - 1
        out1 = name1 + ":  " + str(start1) + " " + partialSeq1[0:60] + " " + str(end1)
        out2 = name2 + ":  " + str(start2) + " "  + partialSeq2[0:60] + " " + str(end2)
        
        f = f + out1 + "\n" + out2 + "\n" + "\n"
        
        start1 = len(partialSeq1NoGap) + start1
        start2 = len(partialSeq2NoGap) + start2
        partialSeq1 = partialSeq1[60:]
        partialSeq2 = partialSeq2[60:]

    f = f + name1 + ":  " + str(start1) + " " + partialSeq1 + " " + str(len(partialSeq1) + start1 - 1) + "\n" + name2 + ":  " + str(start2) + " "  + partialSeq2 + " " + str(len(partialSeq2) + start2 - 1) + "\n" + "\n" + "Score: " + str(score) + "\n" + "Identities: " + str(identities) + "/" + str(len(alignedSeq1)) + " (" + str(round((identities/len(alignedSeq1))*100)) + "%)" + "\n"
    return f

# semi-global sequence alignment function
def semiglobalSeqAlignment(Matrix, seq1, seq2):
    initials = Matrix[0]
    matrix = Matrix[1]
    gapPen = int(Matrix[2])

    table = tableConstruction(len(seq1)+1, len(seq2)+1)

    #Filling out table, seq1 is verticle, seq2 is horizontal
    for i in range(1,len(seq1)+1):
        for j in range(1,len(seq2)+1):
            matching = table[i-1][j-1][0] + int(matrix[initials.index(seq1[i-1])][initials.index(seq2[j-1])])
            if j == len(seq2):
                topGap = table[i-1][j][0]
            else:
                topGap = table[i-1][j][0] + gapPen
            if i == len(seq1):
                botGap = table[i][j-1][0]
            else:
                botGap = table[i][j-1][0] + gapPen
            value = max(matching, topGap, botGap)
            origin = ""
            if value == matching:
                origin = "TL"
            elif value == topGap:
                origin = "T"
            elif value == botGap:
                origin = "L"
            table[i][j] = [value, origin]
    # for row in table:
    #     print(row)
    seq1Char, seq2Char = Backtracking(table, seq1, seq2)

    return seq1Char, seq2Char, 1, table[len(seq1)][len(seq2)][0]

# local sequence alignment function
def localSeqAlignment(Matrix, seq1, seq2):
    initials = Matrix[0]
    matrix = Matrix[1]
    gapPen = int(Matrix[2])

    table = tableConstruction(len(seq1)+1, len(seq2)+1)

    maxVal = 0
    maxValIndexes = []
    #Filling out table, seq1 is verticle, seq2 is horizontal
    for i in range(1,len(seq1)+1):
        for j in range(1,len(seq2)+1):
            matching = table[i-1][j-1][0] + int(matrix[initials.index(seq1[i-1])][initials.index(seq2[j-1])])
            if matching < 0:
                matching = 0
            topGap = table[i-1][j][0] + gapPen
            if topGap < 0:
                topGap = 0
            botGap = table[i][j-1][0] + gapPen
            if botGap < 0:
                botGap = 0
            value = max(matching, topGap, botGap)
            origin = ""
            if value == matching:
                origin = "TL"
            elif value == topGap:
                origin = "T"
            elif value == botGap:
                origin = "L"
            table[i][j] = [value, origin]
            if value > maxVal:
                maxVal = value
                maxValIndexes.clear()
                maxValIndexes.append([i,j])
            elif value == maxVal:
                maxValIndexes.append([i,j])
    # for row in table:
        # print(row)
    
    # list of alignments with same highest score
    list = []
    for i, j in maxValIndexes:
        sequence1, sequence2, score = BacktrackingLocal(i,j, table)
    list.append([sequence1, sequence2,score])
    
    return list

# Separate backtracking for local seq alignment
def BacktrackingLocal(row, col, table):
    aligning = True
    i = row
    j = col
    seq1Char = ""
    seq2Char = ""

    # while cell is no the first one on the table, trace back through i and/or j
    while aligning:
        curr = table[i][j]
        if curr[1] == "" or curr[0] == 0:
            aligning = False
            break
        elif curr[1] == "TL":
            seq1Char = seq1[i-1] + seq1Char
            seq2Char = seq2[j-1] + seq2Char
            i = i-1
            j = j-1
        elif curr[1] == "T":
            seq1Char = seq1[i-1] + seq1Char
            seq2Char = "_" + seq2Char
            i = i-1
        elif curr[1] == "L":
            seq2Char = seq2[j-1] + seq2Char
            seq1Char = "_" + seq1Char
            j = j-1
    for i in range(len(seq1Char)):
        if seq1Char[i] == "_" or seq2Char[i] == "_":
            sequence1 = seq1Char[i+1:]
            sequence2 = seq2Char[i+1:]
        else:
            sequence1 = seq1Char
            sequence2 = seq2Char
            break
        #is start needed anywhere
    return sequence1, sequence2, table[row][col][0]

if __name__ == "__main__":
    protein = True
    alType = "g"

    # Read args
    # i = sequence 1
    # j = sequence 2
    # p = protein (T/F)
    # atype = alignment type [G(lobal)/S(emi-local)/L(ocal)]
    # o = output file name
    if "-i" not in sys.argv or "-j" not in sys.argv or "-p" not in sys.argv or "-atype" not in sys.argv or "-o" not in sys.argv:
        print("missing arguments")
    for i in range(1,len(sys.argv),2):
        if "-" in sys.argv[i+1]:
            print("invalid input")
            break
        if sys.argv[i] == "-i":
            f1 = sys.argv[i+1]
            seq1, name1 = readInput(f1)
            name1 = name1.split()[0]
        elif sys.argv[i] == "-j":
            f2 = sys.argv[i+1]
            seq2, name2 = readInput(f2)
            name2 = name2.split()[0]
        elif sys.argv[i] == "-o":
            outFile = sys.argv[i+1]
            if ".txt" not in sys.argv[i+1]:
                print("invalid output file name")
        elif sys.argv[i] == "-p":
            if sys.argv[i+1] == "F":
                protein = False
            elif sys.argv[i+1] != "T":
                print("invalid -p")
        elif sys.argv[i] == "-atype":
            if sys.argv[i+1] == "S":
                alType = "s"
            elif sys.argv[i+1] == "G":
                alType = "g"
            elif sys.argv[i+1] == "L":
                alType = "l"
            elif sys.argv[i+1] != "G":
                print("invalid -atype")

    # call functions

    # Test
    # seq2 = "AACCTATAGCT"
    # seq1 = "GCGATATA"
    if protein:
        matrix = "BLOSUM45"
    else:
        matrix = "dnaMatrix"
        
    if alType == "g":
        alignedSeq1, alignedSeq2, start, score = globalSeqAlignment(MatrixRead(matrix), seq1, seq2)
        printOut(alignedSeq1, alignedSeq2, seq1, seq2, score, outFile, name1, name2)

    elif alType == "s":
        alignedSeq1, alignedSeq2, start, score = semiglobalSeqAlignment(MatrixRead(matrix), seq1, seq2)
        printOut(alignedSeq1, alignedSeq2, seq1, seq2, score, outFile, name1, name2)

    elif alType == "l":
        localSeqAlignments = localSeqAlignment(MatrixRead(matrix), seq1, seq2)
        # add all alignments with the same score
        for alignment in localSeqAlignments:
            out = printOutLocal(alignment[0], alignment[1], seq1, seq2, alignment[2], name1, name2)
            f = open(outFile, "w")
            f.write(out)

