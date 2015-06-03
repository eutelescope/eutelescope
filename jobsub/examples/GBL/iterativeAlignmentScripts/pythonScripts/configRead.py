import fnmatch
def readConfig(base):
    start=0
    found=False
    parts=[]
    path =  base + '/config/config.cfg'
  #  print path 
    with open(path) as f:
         lines = f.readlines()
    for i,line in enumerate(lines):
        if line.find("[iterativeAlignment]") >= 0 :
            start=i
       #     print "Start line position", start, "  Line: ", line
            found=True
            continue #Continue on to look for end of this processor
            #[*] as wildcard does not work?
            #Look for next processor after iterativeAlignment
        if fnmatch.fnmatch(line, '[*') and  found:
            end=i
         #   print "Here is the line number we end on ", end, " Line: ", line
            break
    if not found:
        print "YOU HAVE NOT ADDED ITERATIVE ALIGNMENT TO THE CONFIG FILE CORRECTLY"
        exit(-1)

 #   print "Here are the variables between: ", start," ",end
    for i in range(start+1,end):
        if not lines[i].isspace():
            print lines[i]
            lines[i] = lines[i].replace('"',"")
            lines[i] = lines[i].replace('\n',"")
            parts.append(lines[i].split("="))
    return parts
