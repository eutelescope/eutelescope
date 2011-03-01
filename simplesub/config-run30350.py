def runConfig(params):
    params["GearFile"]  = " gear_ibl_desy2011_02_16.xml"
    params["HistoInfo"] = "histoinfo.xml"
    params["DUTPlanes"] = "10 20"
    params["BADPlanes"] = ""
    params["TELPlanes"] = "1 2 3 4 5"
    params["TelAlignmentDb"] = "run30350-alignment-tel-db.slcio"
    params["DutAlignmentDb"] = "run30350-alignment-dut-db.slcio"
    params["NoiseDb"]        = "run30350-noise-db.slcio"
    params["DataSet"] = "[30350, 30351, 30352, 30353, 30354, 30355, 30356, 30357, 30358, 30359]"
    params["BeamEnergy"] = "4.0"

def noise(params):
    params["TemplateFile"] = "noise-tmp.xml"

def rawtohit(params):
    params["TemplateFile"] = "rawtohit-tmp.xml"
    params["TelOccupancyThresh"] = "0.001"
    params["DUTOccupancyThresh"] = "0.001"

def aligntel(params):
    params["TemplateFile"] = "kalman-align-tel-tmp.xml"
    params["RecordNumber"] = "100000000"
    params["SkipNEvents"] = "0"
    params["RunPede"] = "True"
    params["UseResidualCuts"] = "True"
    params["ResidualXMin"] = " -150   20  -75 -9999 -9999  -60   20 -200"
    params["ResidualXMax"] = "    0  120   75  9999  9999  100   80  100"
    params["ResidualYMin"] = " -200  -75   50 -9999 -9999 -150  100 -300"
    params["ResidualYMax"] = "   50   75  350  9999  9999  150  200   50"
    params["TelescopeResolution"] = "10 10 10 10000 10000 10 10 10"
    params["DistanceMax"] = "200.0"
    params["MaxChi2"] = "1700"
    params["MinDxDz"] = "-.002"
    params["MaxDxDz"] = ".002"
    params["MinDyDz"] = "-.003" 
    params["MaxDyDz"] = ".001"
    params["ExcludePlanes"] = "10 20"
    params["FixedPlanes"] = "0"
    params["FixedTranslations"] = "4"
    params["FixedScales"] = "2"
    params["FixedZRotations"] = ""

def daffitter(params):
    params["TemplateFile"] = "daf-fitter-tmp.xml"
    params["RecordNumber"] = "10000000"
    params["SkipNEvents"] = "0"
    params["MakePlots"] = "True"
    params["FitDuts"] = "False"
    params["AddToLCIO"] = "True"
    params["TelescopePlanes"] = params["TELPlanes"]
    params["DutPlanes"] = params["DUTPlanes"]
    params["TelResolution"] = "4.3"
    params["DutResolutionX"] = "100"
    params["DutResolutionY"] = "200"
    params["FinderRadius"] = "500.0"
    params["Chi2Cutoff"] = "15.0"
    params["RequireNTelPlanes"] = "5.0"
    params["MaxChi2OverNdof"] = "2.0"
    params["NominalDxdz"] = "-0.0006"
    params["NominalDydz"] = "-0.0016"
    params["ScaleScatter"] = "1.0"
    params["RadiationLengths"] = ""#"0.0007 0.0007 0.0007 0.0007 0.0007 0.0007"
    params["NDutHits"] = "1"

def aligndut(params):
    params["TemplateFile"] = "daf-align-dut.xml"
    params["RecordNumber"] = "10000000"
    params["SkipNEvents"] = "0"
    params["ColMin"] = " 1 1 "
    params["ColMax"] = " 16 78 "
    params["RowMin"] = " 0 10 "
    params["RowMax"] = " 144 325 "
    params["MakePlots"] = "True"
    params["AddToLCIO"] = "True"
    params["TelescopePlanes"] = params["TELPlanes"]
    params["DutPlanes"] = params["DUTPlanes"]
    params["BeamEnergy"] = "2.0"
    params["TelResolution"] = "4.3"
    params["DutResolutionX"] = "7200.0"
    params["DutResolutionY"] = "1000.0"
    params["FinderRadius"] = "700.0"
    params["Chi2Cutoff"] = "12.0"
    params["RequireNTelPlanes"] = "5.0"
    params["MaxChi2OverNdof"] = "1.4"
    params["NominalDxdz"] = "-0.0006"
    params["NominalDydz"] = "-0.0016"
    params["ScaleScatter"] = "1.00"
    params["ResidualXMin"] = " 1800 -5200"
    params["ResidualXMax"] = " 2600 -4650 "
    params["ResidualYMin"] = " 840   6650"
    params["ResidualYMax"] = " 1450  7050"
    params["Translate"] = params["DUTPlanes"]
    params["ZRotate"] = params["DUTPlanes"]
    params["Scale"] = params["DUTPlanes"]
    params["ScaleY"] = ""#params["DUTPlanes"]
    params["ScaleX"] = ""#params["DUTPlanes"]
    params["NDutHits"] = "2"

# run Marlin
if __name__ == "__main__":
    from steeringGenerator import *
    functions = {"noise": noise,
                 "rawtohit": rawtohit,
                 "aligntel": aligntel,
                 "aligndut": aligndut,
                 "daffitter": daffitter}
    jobMaker(functions, runConfig)
