#!/usr/bin/env python
from os import remove
from os import path  #Importing two methods from os module
from re import sub   #This imports regular expression usage
from optparse import OptionParser    #Imports the option parser module

###### OPTIONS and USAGE ######
parser = OptionParser(usage = """fastq_collapse.py -i INFILE -o OUTFILE -j JOB_ID

fastq_collapse.py - Discards exact fastq sequence duplicates, preserving the
    highest average quality score. It also makes a file containing all the 
    duplicated reads.

                                      by
                                Kevin Weitemier
                                 February 2013

Copyright (c) 2011,2012,2013  Kevin Weitemier.
Version 0.03
This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version. A copy of this license is available at <http://www.gnu.org/licenses/>.
Great effort has been taken to make this software perform its said
task, however, this software comes with ABSOLUTELY NO WARRANTY,
not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

Input - A fastq file where each sequence occupies four (4) lines.""")
parser.add_option("-i", action="store", type="string", dest="inname",
    help="Input filename", default="")
parser.add_option("-o", action="store", type="string", dest="outname",
    help="Output filename", default="")
parser.add_option("-j", action="store", type="string", dest="job_id",
    default="1", help="""Unique job ID for performing multiple runs at once in
the same directory. Don't use characters that might cause problems in filenames.
(Default = "1")""")
(options, args) = parser.parse_args()

# Makes sure all filenames are given
if options.inname == "":
    parser.error("Please include an input file using -i.")
if options.outname == "":
    parser.error("Please include an output file using -o.")

###### OPENING INPUT/OUTPUT FILES ######
if path.isfile(options.outname):
    parser.error("""The output filename exists. Please delete it first or \
choose a new name.""")
DupName = "%s_Duplicated_reads" % (options.outname)
if path.isfile(DupName):
    parser.error("""The duplicate reads filename exists. Please delete is first\
 or choose a new name.""")
InFile = open(options.inname, 'r')
OutFile = open(options.outname, 'w')
DupFile = open(DupName, 'w')
Line = InFile.readline()
Bases = ['A','T','C','G','N']
#Temp files
PathDict = {}
HandleDict = {}
for Base in Bases:
    PathDict[Base] = "%s.%s_reads.temp" % (options.job_id,Base)
    HandleDict[Base] = open(PathDict[Base], 'w')

###### ANALYSES ######
###### Parsing reads to temp files ##############
print ("Parsing reads to temp files.")
while Line: #Reads file line-by-line
    #Reads and stores the 4 fastq lines
    AtLine = Line.strip()
    Line = InFile.readline() #This should be the sequence line
    Seq = Line.strip()
    Line = InFile.readline() #This should be the qual ID line
    PlusLine = Line.strip()
    Line = InFile.readline() #This should be the quality line
    Qual = Line.strip()

    for Base, Handle in HandleDict.iteritems():
         if Seq[0] == Base or Seq[0] == Base.lower(): #write the reads to appropriate files.
             Handle.write("%s\n%s\n%s\n%s\n" % (AtLine,Seq,PlusLine,Qual))
             break

    Line = InFile.readline() #This goes to the next @ line to continue loop.

InFile.close()
# Closing the temp files so they can be opened to read later.
for Handle in HandleDict.values():
    Handle.close()

###### ANALYZE EACH TEMP FILE ##################
for Base in Bases:
    print ("Checking for duplicates starting with'%s'." % (Base))
    TempFile = open(PathDict[Base], 'r')
    Line = TempFile.readline()
    SeqDict = {}
    DupDict = {}
    while Line:
        Line = Line.strip()
        SeqID = sub('^@(.+)$',r"\1",Line)
        Line = TempFile.readline() #This should be Seq line
        Seq = Line.strip()
        QualID = TempFile.readline().strip() #This should be qual ID line
        Line = TempFile.readline() #This should be Qual line
        Qual = Line.strip()
        if Seq in SeqDict: #Tests if sequence has already been seen
            #Add to the duplicate dictionary
            DupDict[SeqID] = [Seq,Qual,QualID]
            DupDict[SeqDict[Seq][0]] = [Seq,SeqDict[Seq][1],SeqDict[Seq][2]]
            #Do comparison
            OldScore = 0
            NewScore = 0
            for Position in Qual:
                Score = ord(Position) - 64
                NewScore += Score
            for Place in SeqDict[Seq][1]:
                Score = ord(Place) - 64
                OldScore += Score
            if NewScore > OldScore:
                SeqDict[Seq] = [SeqID,Qual,QualID]
        else:
            SeqDict[Seq] = [SeqID,Qual,QualID]
        Line = TempFile.readline() #Advance to next @line to continue
    
    print ("Writing reads starting with '%s' to file." % (Base))
    for Key, Value in SeqDict.iteritems():
        OutFile.write('@%s\n%s\n%s\n%s\n' % (Value[0],Key,Value[2],Value[1]))
    for Key, Value in DupDict.iteritems():
        DupFile.write('@%s\n%s\n%s\n%s\n' % (Key,Value[0],Value[2],Value[1]))
    TempFile.close()
    SeqDict.clear()
    DupDict.clear()

###### CLEANUP ##########
for Path in PathDict.values():
    remove(Path)
OutFile.close()
DupFile.close()
###### END OF FILE ########
