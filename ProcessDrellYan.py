#!/usr/bin/env python
  
import sys
import subprocess

if (len(sys.argv)!=5) :
    print("ERROR: this script expects exactly 4 arguments")
    print("Variables: "+sys.argv[1]+", "+sys.argv[2]+", "+sys.argv[3]+", "+sys.argv[4])
    sys.exit(0)

runType     = sys.argv[1]
sampleType  = sys.argv[2]
lepType     = sys.argv[3]
fileName    = sys.argv[4]

print("Lepton type: "+lepType)
print("Sample type: "+sampleType)
print("Run type:    "+runType)
print("Input file:  "+fileName)

scriptSpecs = 'test.C+('+ str(runType) + ','+str(sampleType)+','+str(lepType)+','+str(fileName)+')'
print(scriptSpecs)
rootCommand = ['root']
rootCommand.append('-b')
rootCommand.append('-q')
rootCommand.append(scriptSpecs)

(out,err) = subprocess.Popen(rootCommand,stdout=subprocess.PIPE).communicate()
print(out)
print("Any errors?")
print(err)
