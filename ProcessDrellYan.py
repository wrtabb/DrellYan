#!/usr/bin/env python
  
import sys
import subprocess

if (len(sys.argv)!=4) :
    print("ERROR: this script expects exactly 3 arguments")
    print("Variables: "+sys.argv[1]+", "+sys.argv[2]+", "+sys.argv[3])
    sys.exit(0)

runType     = sys.argv[1]
sampleType  = sys.argv[2]
lepType     = sys.argv[3]

print("Lepton type: "+lepType)
print("Sample type: "+sampleType)
print("Run type: "+runType)

scriptSpecs = 'testC('+ str(runType) + ','+str(sampleType)+','+str(lepType)+')'
print(scriptSpecs)
rootCommand = ['root']
rootCommand.append('-b')
rootCommand.append('-q')
rootCommand.append(scriptSpecs)

(out,err) = subprocess.Popen(rootCommand,stdout=subprocess.PIPE).communicate()
print(out)
print("Any errors?")
print(err)
