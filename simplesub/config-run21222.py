def runConfig(params):
    params["GearFile"]  = "gear_ibl_phase1_rot.xml"
    params["HistoInfo"] = "histoinfo.xml"
    params["DUTPlanes"] = "10 11 12"
    params["BADPlanes"] = "13 14 15 16 17"
    params["TELPlanes"] = "0 1 2 3 4 5"
    params["TelAlignmentDb"] = "run21198-alignment-tel-db.slcio"
    params["DutAlignmentDb"] = "run21198-alignment-dut-db.slcio"
    params["NoiseDb"]        = "run21198-noise-db.slcio"
    #Trimmed dataset
    params["DataSet"] = "[21230, 21231, 21232, 21233, 21234, 21235, 21236, 21237, 21238, 21239, 21240, 21241, 21242, 21243, 21244, 21245, 21246, 21247, 21248, 21249, 21250, 21251, 21252, 21253, 21254, 21255, 21256, 21257, 21258, 21259, 21260, 21261, 21262, 21263, 21264, 21265, 21266, 21267, 21268, 21269, 21270, 21271, 21272, 21273, 21274, 21275, 21276, 21277, 21278, 21279, 21280, 21281, 21282, 21283, 21284]"
#Full dataset
#"[21222, 21223, 21224, 21225, 21226, 21227, 21228, 21229, 21230, 21231, 21232, 21233, 21234, 21235, 21236, 21237, 21238, 21239, 21240, 21241, 21242, 21243, 21244, 21245, 21246, 21247, 21248, 21249, 21250, 21251, 21252, 21253, 21254, 21255, 21256, 21257, 21258, 21259, 21260, 21261, 21262, 21263, 21264, 21265, 21266, 21267, 21268, 21269, 21270, 21271, 21272, 21273, 21274, 21275, 21276, 21277, 21278, 21279, 21280, 21281, 21282, 21283, 21284, 21285, 21286, 21287, 21288, 21289, 21290]"#, 21291, 21292, 21293, 21294, 21295, 21296, 21297, 21298, 21299, 21300, 21301, 21302, 21303, 21304, 21305, 21306, 21307]"

def noise(params):
    params["TemplateFile"] = "noise-tmp.xml"

def rawtohit(params):
    params["TemplateFile"] = "rawtohit-tmp.xml"
    params["TelOccupancyThresh"] = "0.001"
    params["DUTOccupancyThresh"] = "0.001"

def aligntel(params):
    params["TemplateFile"] = "kalman-align-tel-tmp.xml"
    params["RecordNumber"] = "1000000000"
    params["SkipNEvents"] = "0"
    params["RunPede"] = "True"
    params["UseResidualCuts"] = "True"
    params["ResidualXMin"] = " 130   -20  -360  -9999  -9999  -9999  -60   80   -40  -9999  -9999  -9999  -9999"
    params["ResidualXMax"] = " 230    80  -270   9999   9999   9999   60  120    40   9999   9999   9999   9999"
    params["ResidualYMin"] = "  60  -180  -160  -9999  -9999  -9999  -80  100  -110  -9999  -9999  -9999  -9999"
    params["ResidualYMax"] = " 200     0   -50   9999   9999   9999   80  150    00   9999   9999   9999   9999"
    params["TelescopeResolution"] = "10 10 10 10000 10000 10000 10 10 10 10000 10000 10000 10000"
    params["DistanceMax"] = "125.0"
    params["MaxChi2"] = "2550"
    params["MinDxDz"] = "-0.0009"
    params["MaxDxDz"] = "0.0001"
    params["MinDyDz"] = "-0.00035" 
    params["MaxDyDz"] = "0.0005"
    params["ExcludePlanes"] = params["DUTPlanes"] + " " + params["BADPlanes"]
    params["FixedPlanes"] = "0"
    params["FixedTranslations"] = "4"
    params["FixedScales"] = "2"
    params["FixedZRotations"] = ""

def aligndut(params):
    params["TemplateFile"] = "daf-align-dut.xml"
    params["RecordNumber"] = "3000000000000"
    params["SkipNEvents"] = "0"
    params["MakePlots"] = "True"
    params["AddToLCIO"] = "False"
    params["TelescopePlanes"] = params["TELPlanes"]
    params["DutPlanes"] = params["DUTPlanes"]
    params["BeamEnergy"] = "120.0"
    params["TelResolution"] = "4.3"
    params["DutResolutionX"] = "1002"
    params["DutResolutionY"] = "1165"
    params["FinderRadius"] = "300.0"
    params["Chi2Cutoff"] = "100.0"
    params["RequireNTelPlanes"] = "4.0"
    params["MaxChi2OverNdof"] = "10.0"
    params["NominalDxdz"] = "-0.0003"
    params["NominalDydz"] = "0.000"
    params["ScaleScatter"] = "1.00"
    params["ResidualXMin"] = "   450  -850  2875 "
    params["ResidualXMax"] = "   620  -500  3030 "
    params["ResidualYMin"] = " -1000  -220 -540 "
    params["ResidualYMax"] = "  -490   340  -60 " 
    params["Translate"] = params["DUTPlanes"]
    params["ZRotate"] = params["DUTPlanes"]
    params["Scale"] = ""#params["DUTPlanes"]
    params["ScaleY"] = ""#params["DUTPlanes"]
    params["ScaleX"] = ""#params["DUTPlanes"]
    params["NDutHits"] = "2"

def fitter(params):
    params["TemplateFile"] = "kalman-fitter-tmp.xml"
    params["RecordNumber"] = "100000000"
    params["SkipNEvents"] = "0"
    params["UseResidualCuts"] = "False"
    #params["InTimeCheck"] = params["DUTPlanes"]
    params["InTimeCheck"] = ""
    params["ResidualXMin"] = "-100  -100  -100   -300 -300  -300  -9999  -100  -100  -100"
    params["ResidualXMax"] = " 100   100   100    300  300   300   9999   100   100   100"
    params["ResidualYMin"] = "-100  -100  -100    -75  -75   -75  -9999  -100  -100  -100"
    params["ResidualYMax"] = " 100   100   100     75   75    75   9999   100   100   100"
    params["MaxChi2"] = "15"
    params["MinDxDz"] = "-.0008"
    params["MaxDxDz"] = ".0002"
    params["MinDyDz"] = "-0.0004" 
    params["MaxDyDz"] = "0.0005"
    params["DistanceMax"] = "40.0"
    params["ExcludePlanes"] = params["BADPlanes"] + " " + params["DUTPlanes"]
    params["AllowNSkippedPlanes"] = "3"
    params["BeamEnergy"] = "120"

def daffittersingle(params):
    daffitter(params)
    params["TemplateFile"] = "daf-fitter-single-tmp.xml"

def daffitter(params):
    params["TemplateFile"] = "daf-fitter-tmp.xml"
    params["RecordNumber"] = "100000"
    params["SkipNEvents"] = "0"
    params["MakePlots"] = "True"
    params["FitDuts"] = "False"
    params["AddToLCIO"] = "True"
    params["TelescopePlanes"] = params["TELPlanes"]
    params["DutPlanes"] = params["DUTPlanes"]
    params["BeamEnergy"] = "120.0"
    params["TelResolution"] = "4.5"
    params["DutResolutionX"] = "10.0"
    params["DutResolutionY"] = "110.0 "
    params["FinderRadius"] = "300.0"
    params["Chi2Cutoff"] = "15.0"
    params["RequireNTelPlanes"] = "6.0"
    params["MaxChi2OverNdof"] = "3.0"
    params["NominalDxdz"] = "-0.0002"
    params["NominalDydz"] = "0.0001"
    params["ScaleScatter"] = "1.0"
    params["NDutHits"]= "0"
    params["RadiationLengths"] = ""#"0.0073 0.0073 0.0073 0.0897 0.0871 0.0846 0.0073 0.0073 0.0073"

def testfitter(params):
    params["TemplateFile"] = "testfitter-tmp.xml"
    params["RecordNumber"] = "100000"
    params["SkipNEvents"] = "0"
    params["AllowedMissingHits"] = "0"
    params["AllowedSkipHits"] = "0"
    params["DistanceMax"] = "1000"
    params["Chi2Max"] = "40"
    params["BeamEnergy"] = "120"
    params["MissingHitPenalty"] = "1"
    params["SkipHitPenalty"] = "1"
    params["PassiveLayerIDs"] = "10 11 12"
    

# run Marlin
if __name__ == "__main__":
    from steeringGenerator import *
    functions = {"noise": noise,
                 "rawtohit": rawtohit,
                 "aligntel": aligntel,
                 "aligndut": aligndut,
                 "daffitter": daffitter,
                 "testfitter": testfitter,
                 "daffittersingle": daffittersingle}
    jobMaker(functions, runConfig)
