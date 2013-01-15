def runConfig(params):
    params["GearFile"]  = "gear_telescope_apix_cern_november2010.xml"
    params["HistoInfo"] = "histoinfo.xml"
    params["DUTPlanes"] = "10 11 12"
    params["BADPlanes"] = "13 14 15 16 17"
    params["TELPlanes"] = "0 1 2 3 4 5"
    params["TelAlignmentDb"] = "run21198-alignment-tel-db.slcio"
    params["DutAlignmentDb"] = "run21198-alignment-dut-db.slcio"
    params["NoiseDb"]        = "run21198-noise-db.slcio"
    params["DataSet"] = "[21198, 21199, 21200, 21201, 21202, 21203, 21204, 21205, 21206, 21207, 21208, 21209, 21210, 21211, 21212, 21213, 21214, 21215, 21216, 21217, 21218, 21219, 21220, 21221]"

def noise(params):
    params["TemplateFile"] = "noise-tmp.xml"

def rawtohit(params):
    params["TemplateFile"] = "rawtohit-tmp.xml"
    params["TelOccupancyThresh"] = "0.001"
    params["DUTOccupancyThresh"] = "0.001"

def aligntel(params):
    params["TemplateFile"] = "kalman-align-tel-tmp.xml"
    params["RecordNumber"] = "10000000"
    params["SkipNEvents"] = "0"
    params["RunPede"] = "True"
    params["UseResidualCuts"] = "True"
    params["ResidualXMin"] = " 120   -20  -360  -9999  -9999  -9999  -60   60   -50  -9999  -9999  -9999  -9999"
    params["ResidualXMax"] = " 230    80  -260   9999   9999   9999   60  130    60   9999   9999   9999   9999"
    params["ResidualYMin"] = "  60  -180  -160  -9999  -9999  -9999  -80   90  -120  -9999  -9999  -9999  -9999"
    params["ResidualYMax"] = " 220     0   -50   9999   9999   9999   80  160    10   9999   9999   9999   9999"
    params["TelescopeResolution"] = "10 10 10 10000 10000 10000 10 10 10 10000 10000 10000 10000"
    params["DistanceMax"] = "125.0"
    params["MaxChi2"] = "2700"
    params["MinDxDz"] = "-0.0009"
    params["MaxDxDz"] = "0.0001"
    params["MinDyDz"] = "-0.0004" 
    params["MaxDyDz"] = "0.0005"
    params["ExcludePlanes"] = params["DUTPlanes"] + " " + params["BADPlanes"]
    params["FixedPlanes"] = "0"
    params["FixedTranslations"] = "4"
    params["FixedScales"] = "2"
    params["FixedZRotations"] = ""

def aligndut(params):
    params["TemplateFile"] = "daf-align-dut.xml"
    params["RecordNumber"] = "1000000"
    params["SkipNEvents"] = "0"
    params["MakePlots"] = "True"
    params["FitDuts"] = "False"
    params["AddToLCIO"] = "True"
    params["TelescopePlanes"] = params["TELPlanes"]
    params["DutPlanes"] = params["DUTPlanes"]
    params["BeamEnergy"] = "120.0"
    params["TelResolution"] = "4.3"
    params["DutResolutionX"] = "144"
    params["DutResolutionY"] = "1165"
    params["FinderRadius"] = "300.0"
    params["Chi2Cutoff"] = "100.0"
    params["RequireNTelPlanes"] = "4.0"
    params["MaxChi2OverNdof"] = "12.0"
    params["NominalDxdz"] = "-0.0003"
    params["NominalDydz"] = "0.000"
    params["ScaleScatter"] = "1.00"
    params["ResidualXMin"] = "  1270   600   4825"
    params["ResidualXMax"] = "  1425   800   4960"
    params["ResidualYMin"] = " -1000  -280   -560"
    params["ResidualYMax"] = " -500   280    -90"
    params["Translate"] = params["DUTPlanes"]
    params["ZRotate"] = params["DUTPlanes"]
    params["Scale"] = params["DUTPlanes"]
    params["ScaleY"] = ""#params["DUTPlanes"]
    params["ScaleX"] = ""# params["DUTPlanes"]
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

def daffitter(params):
    params["TemplateFile"] = "daf-fitter-tmp.xml"
    params["RecordNumber"] = "10000000"
    params["SkipNEvents"] = "0"
    params["MakePlots"] = "True"
    params["FitDuts"] = "False"
    params["AddToLCIO"] = "True"
    params["TelescopePlanes"] = params["TELPlanes"]
    params["DutPlanes"] = params["DUTPlanes"]
    params["BeamEnergy"] = "120.0"
    params["TelResolution"] = "4.3"
    params["DutResolutionX"] = "7.2"
    params["DutResolutionY"] = "58.0 "
    params["FinderRadius"] = "300.0"
    params["Chi2Cutoff"] = "40.0"
    params["RequireNTelPlanes"] = "6.0"
    params["MaxChi2OverNdof"] = "7.0"
    params["NominalDxdz"] = "-0.0002"
    params["NominalDydz"] = "0.0001"
    params["ScaleScatter"] = ".5"
    params["NDutHits"]= "1"
    params["RadiationLengths"] = "0.0073 0.0073 0.0073 0.0897 0.0871 0.0846 0.0073 0.0073 0.0073"

# run Marlin
if __name__ == "__main__":
    from steeringGenerator import *
    functions = {"noise": noise,
                 "rawtohit": rawtohit,
                 "aligntel": aligntel,
                 "aligndut": aligndut,
                 "daffitter": daffitter}
    jobMaker(functions, runConfig)
