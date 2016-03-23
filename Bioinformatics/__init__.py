__author__ = 'Megan Nilles'
#due Feb 12, 2015
#OYOP6
'''This program uses a semi-global alignment to align sequences
because it does not penalize gaps for unequal length sequences
It will then calculate the distance metric using the Jukes-Cantor,
Kimura or Tamura. It is user dependent.
'''

import math

#functions:

#function that builds alignment from directional string
def buildDirectionalString (s1, s2, n ,m ,dString):
    #step 3 Build alignment using directional strings
    s1Pos = n-1
    s2Pos = m-1
    dirPos = 0
    alS1 = ''
    alS2 = ''
    alignmentScore = 0
    matrixOut = ''

    while(dirPos < len(dString)):

        if(dString[dirPos] == "D"): #Align sequence for match
            alS1 = s1[s1Pos] + alS1
            alS2 = s2[s2Pos] + alS2
            s1Pos = s1Pos - 1
            s2Pos = s2Pos - 1
        elif(dString[dirPos] == "V"):#Place a space for a gap
            alS1 = s1[s1Pos] + alS1
            alS2 = '-' + alS2
            s1Pos = s1Pos - 1
        else:#H
            alS1 = '-' + alS1 #Place a space for a gap
            alS2 = s2[s2Pos] + alS2
            s2Pos = s2Pos -1

        dirPos = dirPos + 1

    alS1List = list(alS1)
    alS2List = list(alS2)
    return (alS1List,alS2List,alS1,alS2)


# Function that builds a string that provides direction for tracing
# a path back through the matrix to build the best match.
def buildDirectional(matrix, s1Len, s2Len,gap):
    dString = ''
    currentRow = s1Len
    currentCol = s2Len
    gapScore = gap
    while(currentRow != 0 or currentCol != 0):
        if(currentRow == 0):
            dString = dString + ('H')
            currentCol = currentCol - 1
        elif(currentCol == 0):
            dString = dString + ('V')
            currentRow = currentRow - 1
        elif(matrix[currentRow][currentCol - 1] + int(gapScore) == matrix[currentRow][currentCol]):
            dString = dString + ('H')
            currentCol = currentCol - 1
        elif(matrix[currentRow - 1][currentCol] + int(gapScore) == matrix[currentRow][currentCol]):
            dString = dString + ('V')
            currentRow = currentRow - 1
        else:
            dString = dString + ('D')
            currentRow = currentRow - 1
            currentCol = currentCol - 1
    return dString

#This is the alignment function. It will score and print the percentages if the user wants.
def semiGlobal(s1Len, s2Len,match, mismatch):
    #step 1:build matrix
    #variable declaration
    m = s2Len
    n = s1Len
    #build matrix
    matrix = [[0 for x in range(m+1)] for x in range(n+1)]
    for i in range(1, n+1):
        matrix[i][0] = matrix[0][0]
    for j in range(1, m+1):
        matrix[0][j] = matrix[0][0]

    for i in range(1, n+1):
        for j in range(1, m+1):
            if s1[i-1] == s2 [j-1]:
                score1 = matrix[i-1][j-1] + int(match)
            else:
                score1 = matrix[i-1][j-1] + int(mismatch)
            score2 = matrix[i][j-1]
            score3 = matrix[i-1][j]
            matrix[i][j] = max(score1,score2,score3)

    score = matrix[n][m]#score overall
    return matrix
    return score


#Jukes-Cantor function for distance calculation
def jukes(alS1, alS2):
    difCtr = 0
    lenCtr = 0

    for i in range(0,len(s1)):
        if(alS1[i] != '-' and alS2[i] != '-'):
            lenCtr += 1
        if(alS1[i] != alS2[2]):
            difCtr += 1
    difCtr = float(difCtr)/float(lenCtr)
    jukes = abs((-3/4) * math.log(abs(1-(4/3)*difCtr)))
    print jukes


#Kimura for distance calculation
def kimura(alS1, alS2):
    difCtr = 0
    lenCtr = 0
    s = 0 #transitions
    v = 0#transversions
    for i in range(0,len(s1)):
        if(alS1[i] != '-' and alS2[i] != '-'):
            lenCtr += 1
        if(alS1[i] != alS2[2]):
            difCtr = difCtr+1
            if alS1[i] == 'G':
                if alS2[i] == 'A':
                    s = s+1
                else:#transversion
                    v = v+1
            elif alS1[i]=='A':
                if alS2[i] == 'G':
                    s = s+1
                else:#transversion
                    v = v+1
            elif alS1[i] == 'T':
                if alS2[i] == 'C':
                    s = s+1
                else:#transversion
                    v = v+1
            else:#C
                if alS2[i] == 'T':
                    s = s+1
                else:#transversion
                    v = v+1
            difCtr +1
    s = float(s/difCtr)
    v = float(s/difCtr)
    kimura = ((1/2)*math.log(abs(1/(1-(2*s)-v))))+((1/4)*math.log(abs(1/(1-(2*v)))))
    print kimura

#Tamura for distance calculation
def tamura(alS1, alS2):
    difCtr = 0
    lenCtr = 0
    s = 0 #transitions
    v = 0#transversions
    #transitions and transversions
    for i in range(0,len(s1)):
        if(alS1[i] != '-' and alS2[i] != '-'):
            difCtr += 1
        if(alS1[i] != alS2[2]):
            difCtr=difCtr + 1
            if alS1[i] == 'G':
                if alS2[i] == 'A':
                    s = s+1
                else:#transversion
                    v = v+1
            elif alS1[i]=='A':
                if alS2[i] == 'G':
                    s = s+1
                else:#transversion
                    v = v+1
            elif alS1[i] == 'T':
                if alS2[i] == 'C':
                    s = s+1
                else:#transversion
                    v = v+1
            else:#C
                if alS2[i] == 'T':
                    s = s+1
                else:#transversion
                    v = v+1
            difCtr +1
    s = float(s/difCtr)
    v = float(s/difCtr)
    #determine C for tamura
    c1=s1.count('C')
    g1=s1.count('G')
    c2=s2.count('C')
    g2=s2.count('G')
    gc1 = float((c1 + g1)/len(s1))
    gc2 = float((c2 + g2)/len(s2))
    c = float((gc1 + gc2)-(2*gc1*gc2))

    tamura = ((-c)*math.log(abs(1-(s/c)-v)))-(1/2(1-c))*math.log(abs(1-(2*v)))
    print tamura

#this is main, but not called main
#user input for scoring data and sequences
print ('This program will align your sequences and output the distance scoring in whichever method you choose. \n')
data = raw_input ('Are your sequences in a file?\n')
#get names of files and read in lines or direct input
if data == 'Y' or data == 'y':
    infile1 = raw_input ('Please enter the name of your file including the .txt \n')
    print ('You entered ', infile1)
    infile1 = open(infile1, 'r')
    data = raw_input('Is the second sequence in another file?\n')
    if data == 'Y' or data == 'y':
        infile2 = raw_input ('Please enter the name of your second file including the .txt \n')
        print ('You entered ',infile2)
        infile2 = open(infile2, 'r')

    s1 = ""
    infile1.readline() #bypass > headers
    for line in infile1:
        line = line.replace('\n','') #takes out new line formatting
        s1 = s1 + line
    s1 = s1.upper() #make it uppercase

    s2 = ""
    infile2.readline() #bypass > header line
    for line in infile2:
        line = line.replace('\n','') #takes out new line formatting here too
        s2 = s2 + line
    s2 = s2.upper()


else:
    s1 = raw_input('Please enter your sequence. ')
    s2 = raw_input('Please enter your second sequence. ')
    s1 = s1.upper()
    s2 = s2.upper()

#align sequence
s1Len = len(s1)
s2Len = len(s2)

match = raw_input('Please enter your match score. ')
mismatch = raw_input('Please enter your mismatch score. ')
gap = raw_input('Please enter your gap score. ')

#step1: build matrix
matrix = semiGlobal(s1Len, s2Len, match, mismatch)


#step 2 create directional string
dString = ''
dString = buildDirectional(matrix, s1Len, s2Len,gap)

#step3: build directional string
alS1List,alS2List,alS1, alS2 = buildDirectionalString (s1, s2, s1Len,s2Len,dString)


#output for aligned sequences
connectors = ''
alignmentScore = 0
for char1, char2 in zip(alS1List, alS2List):
    if char1 == char2:
        connectors = connectors + "|" #Create list for relationship between sequences
        alignmentScore = alignmentScore + 1

    else:
        if char1 == "-" or char2 == "-":
            connectors = connectors + " "
        else:
            connectors = connectors + "."

line = 0
lineCount = 1
charPerLine = 50
# Prints the alignment with a limit of 50 characters per line.
for line in range(line,len(alS1),charPerLine):
    print(str(lineCount) + " " + alS1[line:line + charPerLine])
    print(len(str(lineCount)) * " " + " " + connectors[line:line + charPerLine])
    print(str(lineCount) + " " + alS2[line:line+charPerLine])
    lineCount = lineCount + 1

print("Match score is: " + str(alignmentScore/float(len(alS1))))
print("\nMatch percentage is: " + str(alignmentScore/float(len(alS1))) + '\n')
print("The alignment score is: " + str(matrix[s1Len][s2Len]))

answer = raw_input('would you like to know the distance? ')
while answer == 'Y' or answer == 'y':
    print ('What metric would you like to use? ')
    data = raw_input('(J)ukes-Cantor, (K)imura, or (T)amura?')
    if data == 'J' or data == 'j':
        jukes(alS1,alS2)
    elif data == 'K' or data == 'k':
        kimura(alS1,alS2)
    else:#T
        tamura(alS1,alS2)
    answer = raw_input('Would you like to see them in a different metric? ')
