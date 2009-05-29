#!/usr/bin/env python

from ROOT import *
import lcio
from optparse import OptionParser

def main() :


    cvsVersion = "$Revision: 1.1 $"
    version = cvsVersion[10:len(cvsVersion) - 1 ]

    parser = OptionParser( version=version )

    parser.add_option( "-o", "--output-root-file",
                       action="store",
                       dest="output",
                       help="Use this option to specify the ROOT output file name")

    parser.add_option("-t", "--template",
                      action="store",
                      dest="template",
                      help="Use this option to specify the input file template in this format myfile-X.slcio where"
                      " XXX will be replaced with the correct value")

    parser.add_option("-f", "--full-statistics",
                      action="store",
                      dest="full_stats",
                      help="Use this option to specify the file contaning the full statistics" )

    parser.add_option("--min",
                      action="store",
                      dest="minvalue",
                      type="int",
                      help="Minimum value of the scan")

    parser.add_option("--max",
                      action="store",
                      dest="maxvalue",
                      type="int",
                      help="Maximum value of the scan")

    parser.add_option("--step",
                      action="store",
                      dest="step",
                      type="int",
                      help="Incremental step")

    options , args = parser.parse_args()


    if options.output == None:
        options.output = "file.root"

    outputFile = TFile(  options.output , "RECREATE" )

    lcioFactory = lcio.LCFactory.getInstance()
    dataReader  = lcioFactory.createLCReader()
    scaleFactor = 1000

    if options.full_stats != None:

        dataReader.open( options.full_stats )
        event = dataReader.readNextEvent()
        collection = event.getCollection( "alignment" )
        fullXShiftGraphList = []
        fullYShiftGraphList = []
        fullPhyGraphList    = []

        for pos in range( collection.getNumberOfElements() ):

            element = collection.getElementAt( pos )

            xShift     = scaleFactor * element.getDoubleVal(  0 )
            yShift     = scaleFactor * element.getDoubleVal(  1 )
            phy        = element.getDoubleVal(  5 )
            xShiftErr  = scaleFactor * element.getDoubleVal(  6 )
            yShiftErr  = scaleFactor * element.getDoubleVal(  7 )
            phyErr     = element.getDoubleVal( 11 )

            xShiftGraph = TGraph( 5 )
            xShiftGraph.SetPoint( 0 , options.minvalue, xShift - 0.5 * xShiftErr )
            xShiftGraph.SetPoint( 1 , options.maxvalue, xShift - 0.5 * xShiftErr )
            xShiftGraph.SetPoint( 2 , options.maxvalue, xShift + 0.5 * xShiftErr )
            xShiftGraph.SetPoint( 3 , options.minvalue, xShift + 0.5 * xShiftErr )
            xShiftGraph.SetPoint( 4 , options.minvalue, xShift - 0.5 * xShiftErr )

            xShiftGraph.SetFillColor ( kBlue )
            xShiftGraph.SetFillStyle( 3005 )

            fullXShiftGraphList.append( xShiftGraph )

            yShiftGraph = TGraph( 5 )
            yShiftGraph.SetPoint( 0 , options.minvalue, yShift - 0.5 * yShiftErr )
            yShiftGraph.SetPoint( 1 , options.maxvalue, yShift - 0.5 * yShiftErr )
            yShiftGraph.SetPoint( 2 , options.maxvalue, yShift + 0.5 * yShiftErr )
            yShiftGraph.SetPoint( 3 , options.minvalue, yShift + 0.5 * yShiftErr )
            yShiftGraph.SetPoint( 4 , options.minvalue, yShift - 0.5 * yShiftErr )

            yShiftGraph.SetFillColor ( kBlue )
            yShiftGraph.SetFillStyle( 3005 )

            fullYShiftGraphList.append( yShiftGraph )

            phyGraph = TGraph( 5 )
            phyGraph.SetPoint( 0 , options.minvalue, phy - 0.5 * phyErr )
            phyGraph.SetPoint( 1 , options.maxvalue, phy - 0.5 * phyErr )
            phyGraph.SetPoint( 2 , options.maxvalue, phy + 0.5 * phyErr )
            phyGraph.SetPoint( 3 , options.minvalue, phy + 0.5 * phyErr )
            phyGraph.SetPoint( 4 , options.minvalue, phy - 0.5 * phyErr )

            phyGraph.SetFillColor ( kBlue )
            phyGraph.SetFillStyle( 3005 )

            fullPhyGraphList.append( phyGraph )
            

    xShiftGraphList = []
    yShiftGraphList = []
    phyGraphList    = []
    for index in range( 6 ):

        xShiftGraph = TGraphErrors( len( range( options.minvalue, options.maxvalue, options.step ) ))
        name = "xShift_d%(i)d" % { "i": index }
        title = name
        xShiftGraph.SetName( name )
        xShiftGraph.SetTitle( title )
        xShiftGraph.GetXaxis().SetTitle( "time [kEvents]")
        xShiftGraph.GetYaxis().SetTitle( "shift [#mum]")

        yShiftGraph = TGraphErrors( len( range( options.minvalue, options.maxvalue, options.step ) ))
        name = "yShift_d%(i)d" % { "i": index }
        title = name
        yShiftGraph.SetName( name )
        yShiftGraph.SetTitle( title )
        yShiftGraph.GetXaxis().SetTitle( "time [kEvents]")
        yShiftGraph.GetYaxis().SetTitle( "shift [#mum]")


        phyGraph    = TGraphErrors( len( range( options.minvalue, options.maxvalue, options.step ) ))
        name = "phy_d%(i)d" % { "i": index }
        title = name
        phyGraph.SetName( name )
        phyGraph.SetTitle( title )
        phyGraph.GetXaxis().SetTitle( "time [kEvents]")
        phyGraph.GetYaxis().SetTitle( "rotation [mrad]")


        xShiftGraphList.append( xShiftGraph )
        yShiftGraphList.append( yShiftGraph )
        phyGraphList.append( phyGraph )



    template = options.template.replace( "X", "%(value)s" )

    index = 0
    for value in range( options.minvalue, options.maxvalue, options.step ) :
        file = template % { "value": value }

        dataReader.open( "%(file)s" % { "file": file } )
        event = dataReader.readNextEvent()
        collection = event.getCollection( "alignment" )

        for pos in range( collection.getNumberOfElements() ):
            if pos == 0 or pos == collection.getNumberOfElements() - 1:
                continue


            element = collection.getElementAt( pos )

            xShift     = scaleFactor * element.getDoubleVal(  0 )
            yShift     = scaleFactor * element.getDoubleVal(  1 )
            phy        = element.getDoubleVal(  5 )
            xShiftErr  = scaleFactor * element.getDoubleVal(  6 )
            yShiftErr  = scaleFactor * element.getDoubleVal(  7 )
            phyErr     = element.getDoubleVal( 11 )

            xShiftGraphList[pos].SetPoint(  index, value, xShift )
            xShiftGraphList[pos].SetPointError( index, 0, xShiftErr )

            yShiftGraphList[pos].SetPoint(  index, value, yShift )
            yShiftGraphList[pos].SetPointError( index, 0, yShiftErr )

            phyGraphList[pos].SetPoint( index, value, phy )
            phyGraphList[pos].SetPointError( index, 0, phyErr )

        dataReader.close()
        index = index + 1

    for index, value in enumerate( xShiftGraphList ) :
        xShiftGraphList[ index ].Write()
        yShiftGraphList[ index ].Write()
        phyGraphList[ index ].Write()



    # prepare the canvas
    xShiftCanvas = TCanvas("xShift", "xShift", 800, 800 )
    xShiftCanvas.Divide( 2 , 2 )
    padCounter = 0
    for index, value in enumerate( xShiftGraphList ) :
        if index == 0 or index == len( xShiftGraphList ) - 1 :
            continue
        padCounter = padCounter + 1
        xShiftCanvas.cd( padCounter )
        value.Draw( "ALP " )
        if options.full_stats != None:
            fullXShiftGraphList[ index ].Draw("F")
    xShiftCanvas.Write()

    # prepare the canvas
    yShiftCanvas = TCanvas("yShift", "yShift", 800, 800 )
    yShiftCanvas.Divide( 2 , 2 )
    padCounter = 0
    for index, value in enumerate( yShiftGraphList ) :
        if index == 0 or index == len( yShiftGraphList ) - 1 :
            continue
        padCounter = padCounter + 1
        yShiftCanvas.cd( padCounter )
        value.Draw( "ALP " )
        if options.full_stats != None:
            fullYShiftGraphList[ index ].Draw("F")
    yShiftCanvas.Write()


    # prepare the canvas
    phyCanvas = TCanvas("phy", "phy", 800, 800 )
    phyCanvas.Divide( 2 , 2 )
    padCounter = 0
    for index, value in enumerate( phyGraphList ) :
        if index == 0 or index == len( phyGraphList ) - 1 :
            continue
        padCounter = padCounter + 1
        phyCanvas.cd( padCounter )
        value.Draw( "ALP " )
        if options.full_stats != None:
            fullPhyGraphList[ index ].Draw("F")
    phyCanvas.Write()

    raw_input( "any key to continue... " )
    outputFile.Close()


if __name__ == "__main__" :
    main()
