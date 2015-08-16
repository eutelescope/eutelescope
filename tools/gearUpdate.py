#!/usr/bin/env python2
'''
Created on Jul 8, 2015

@author: A.Morton
'''
import argparse
from xml.dom.minidom import parse 
import sys
import shutil
## \file
#This file will use the results of millepede to update the geometry file. The location of the new gear must also be given.
#It will update the gear file using the CMS convention.
#Assumed parameters labelling: label = 10*(sensor ID) + (parameterNumber) + 1
#X/Y/Z shifts => parameterNumber (0,1,2)
#X/Y/Z rotations => parameterNumber (3,4,5)


radianCutOff = 0.1 

## Create results vector
#
# This function will read in a results file and return a vector with alignment parameters,results and errors   
#
# @param[in] res Results file 
# @return resLink Output [alignment parameter,results value,error]   
#
def getVector(res):
    vecParRes = []
    with open(res) as result:
        resLines = result.readlines()
    for line in resLines:
        parRes = getParRes(line)
        if parRes is not None:
            vecParRes.append(parRes)
    return vecParRes
## Get alignment parameters. 
#
# The function will return the alignment parameter, result and error if there is one on this line in the results file   
# It will return nothing if the absolute value is zero
#
# @param[in] line  
# @return parRes list [alignment parameter,results value]   
#
def getParRes(line):        
    try:
        parts = line.split()
        ##Check the first string is an alignment parameter
        float(parts[0])
        if(abs(float(parts[1])) > 0):
            return [parts[0],parts[1],parts[4]]  #alignment ID, result, error.
        else:
           return None 
    except ValueError:
        return None
    except IndexError:
        print "No error give for millepede!!!!! Check Steering file has solution method added"
        return None

## Get information associated to this alignment parameter 
#
# Millepede results file is read in and transformed to list to be placed inside gear file   
# This will find the alignment mode (shift, rotation etc) and ID from the alignment parameters.
#
# @param[in] resVec Vector with alignment ID and results value
# @return resLink Output [ID,Gear alignment string,results value]   
#
##\todo We reverse the string ID at the end. This might be a bit convoluted. I think this can be improved but would have to refactor this function.

def getInfo(resVec):
    IDModeRes = []
    link = [ "positionX" , "positionY","positionZ","rotationZY","rotationZX","rotationXY"] 
    for res in resVec:
        par = list(res[0]);
        ID = ""
        mode = ""
        for count, num in enumerate(reversed(par)):
            if(count == 0 and len(par) == 1): #If we have a single number this must be sensor 0 and the mode is just the number
                mode = link[int(num) - 1]
                ID = str(0)
            elif(count == 0 and  len(par) != 1 ): #If the length is greater than 1 then the mode is just the first number in the reversed list. 
                mode = link[int(num) -1]
            else:
                ID = ID + num
        ID = ID[::-1]
        IDModeRes.append([ID,mode,res[1],res[2]]) 
    return IDModeRes

## Update the gear  
#
#   The gear is read using a xml parse and the elements updated. 
#   Copy the old gear then access the newgear which is a copied and updated
#
# @param[in] IDModeRes [ID,mode,result,error]
# @param[in] oldGear 
# @return newGear Gear with added results returned   
#

def updateGear(IDModeRes,oldGear,newGear):
    try:
        #This is a list of the strings we should change from radians to degrees. 
        trans = [ "rotationZY","rotationZX","rotationXY"] 

        shutil.copyfile(oldGear, newGear)
        g = parse(newGear)

        for imd in IDModeRes:
            layers = g.getElementsByTagName("layer")
            for i, layer in enumerate(layers):
                ladder = layer.getElementsByTagName("ladder")[0]
                id = int(ladder.getAttribute("ID")) 
                if(int(id) == int(imd[0])):
                    val =  ladder.getAttribute(imd[1])
                    ##If the result is a rotation, multiply by 57.3(180/pi)
                    if imd[1] in trans:
                        newVal = float(val) - float(imd[2])*57.3
                    else:
                        newVal = float(val) - float(imd[2])

                    ladder.setAttribute(str(imd[1]), str(newVal))
        g.writexml(open(newGear, 'w'))
    except  shutil.Error:
        print "ERROR: YOU ARE USING THE SAME FILE AS OUTPUT AS INPUT!!!!!!!!!!!!!!!!! GEAR FILE NOT UPDATED!!!!!!"

## Screen results for large rotations and small results compared to error. 
#  These will be removed and not passed back via return.
# @param[in] IDModeResErr [ID,mode,result,error]
# @return   IDModeResErr [ID,mode,result,error] updated to remove large rotations and small result values
#

def screen(IDModeRes):
    IDModeResUpdate = []
    #This is a list of the strings we should check for large rotations. 
    trans = [ "rotationZY","rotationZX","rotationXY"] 
    for imd in IDModeRes:
        if abs(float(imd[2])/float(imd[3])) > 1 : #Check the update is large than the error. Otherwise there is no point in saving.
            if imd[1] in trans: #If a rotation then check size. Otherwise just add.
                if abs(float(imd[2])) > radianCutOff:
                    print "WARNING large rotation. ID: " , imd[0], " mode: ", imd[1], " result ", imd[2], " error ", imd[3], "    Removed! "  
                else:
                    IDModeResUpdate.append(imd) 
            else:
                IDModeResUpdate.append(imd) 

        else:
            print "INFO change smaller than error. ID: " , imd[0], " mode: ", imd[1], " result ", imd[2], " error ", imd[3], "    Removed! "  

    return IDModeResUpdate

#We might want to import this at some point. Remove this then.
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-og', help="The old gear file", default="gear.xml")
    parser.add_argument('-ng', help="The new gear file location to save. (Give new name too)", default="gear.xml")
    parser.add_argument('-r', help="The results file",default="millepede.res" )
    args = parser.parse_args()
    resVec = getVector(args.r)
    IDModeRes =  getInfo(resVec)
    IDModeResCor =  screen(IDModeRes)

    print "The alignment parameters to update with these corrections:"
    for imr in IDModeResCor:
        print "The sensor ID: ", imr[0] , " Type: ", imr[1], " Value: ", imr[2] , " Error: ", imr[3]
    updateGear(IDModeResCor,args.og,args.ng)
