import os
#This function will find the variable in config and look for it to the left of equals sign.
#If the string is defined by another variable then look this up and output the correct string
def findVar(var):
    parts = []
    with open('config/config.cfg') as f:
        lines = f.readlines()
        for i,line in enumerate(lines):
         #   print "line: " , line
            if line.find(var) >= 0 :
                print "Found String: ",line
                parts = line.split("=")
                if parts[0].find(var) >= 0:
                    print "Found equality"
                    parts[1] = replaceSymbols(parts[1])
                    return parts[0], parts[1]





def replaceSymbols(input):
    eutel = os.environ['EUTELESCOPE']
    print "EUTel location: ", eutel
    print "Before ", input
    input = input.replace("%(eutelescopepath)s",eutel)
    print "After: ", input 
    return input


