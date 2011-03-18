def runConfig(params):
    params["GearFile"]  = " gear_ibl_desy_0_0.xml"
    params["HistoInfo"] = "histoinfo.xml"
    params["DUTPlanes"] = "10 20"
    params["BADPlanes"] = ""
    params["TELPlanes"] = "0 1 2 3 4 5"
    params["TelAlignmentDb"] = "run30377-alignment-tel-db.slcio"
    params["MilleAlignmentDb"] = "run30377-alignment-mille-db.slcio"
    params["DutAlignmentDb"] = "run30377-alignment-dut-db.slcio"
    params["NoiseDb"]        = "run30377-noise-db.slcio"
    params["PreAlignDb"]     = "run30377-prealign-db.slcio"
    params["DataSet"] = "[30377]"
    params["BeamEnergy"] = "2.0"

def noise(params):
    params["TemplateFile"] = "noise-pysub-tmp.xml"
    params["TelOccupancyThresh"] = "0.001"

def rawtohit(params):
    params["TemplateFile"] = "rawtohit-pysub-tmp.xml"

def prealign(params):
    params["TemplateFile"] = "prealigner-tmp.xml"
    params["RecordNumber"] = "100000"
    params["SkipNEvents"] = "0"

def aligntel(params):
    params["TemplateFile"] = "kalman-align-tel-pre-tmp.xml"
    params["RecordNumber"] = "100000000"
    params["SkipNEvents"] = "0"
    params["RunPede"] = "True"
    params["UseResidualCuts"] = "True"
    params["ResidualXMin"] = " -120  -120  -120  -120  -120  -120  -120  -120  -120 "
    params["ResidualXMax"] = "  120   120   120   120   120   120   120   120   120 "
    params["ResidualYMin"] = " -120  -120  -120  -120  -120  -120  -120  -120  -120 "
    params["ResidualYMax"] = "  120   120   120   120   120   120   120   120   120 "
    params["TelescopeResolution"] = "10 10 10 10000 10000 10 10 10"
    params["DistanceMax"] = "200.0"
    params["MaxChi2"] = "250"
    params["MinDxDz"] = "-.003"
    params["MaxDxDz"] = ".003"
    params["MinDyDz"] = "-.003" 
    params["MaxDyDz"] = ".003"
    params["ExcludePlanes"] = params["DUTPlanes"]
    params["FixedPlanes"] = "0"
    params["FixedTranslations"] = "4"
    params["FixedScales"] = "2"
    params["FixedZRotations"] = ""

def aligndut(params):
    params["TemplateFile"] = "daf-align-dut-pre.xml"
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
    params["NominalDxdz"] = "-0.0"
    params["NominalDydz"] = "-0.0004"
    params["ScaleScatter"] = "1.00"
    params["ResidualXMin"] = "-200  -250"
    params["ResidualXMax"] = " 450   270 "
    params["ResidualYMin"] = "-120  -250"
    params["ResidualYMax"] = " 170   300"
    params["Translate"] = params["DUTPlanes"]
    params["ZRotate"] = params["DUTPlanes"]
    params["Scale"] = params["DUTPlanes"]
    params["ScaleY"] = ""#params["DUTPlanes"]
    params["ScaleX"] = ""#params["DUTPlanes"]
    params["NDutHits"] = "2"

def millealign(params):
    params["TemplateFile"] = "eutelmille-align-tmp.xml"
    params["RunPede"] = "1"
    params["RecordNumber"] = "100000000"
    params["SkipNEvents"] = "0"
    params["UseResidualCuts"] = "1"
    params["ResidualXMin"] = " -80  -50  -50  -200  -250  -250  -400  -400 "
    params["ResidualXMax"] = "  80   50   50   500   300   450   600   600 "
    params["ResidualYMin"] = " -50  -50  -50  -150  -200  -400  -600  -500 "
    params["ResidualYMax"] = "  50   50   50   150   200   500   600   700 "
    params["DistanceMax"] = "1000"
    params["ExcludePlanes"] = ""
    params["FixedPlanes"] = "0 5"

def dafcommon(params):
    params["RecordNumber"] = "1000000"
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
    params["Chi2Cutoff"] = "20.0"
    params["RequireNTelPlanes"] = "4.0"
    params["MaxChi2OverNdof"] = "2.0"
    params["NominalDxdz"] = "0.0001"
    params["NominalDydz"] = "-0.0004"
    params["ScaleScatter"] = "1.0"
    params["RadiationLengths"] = ""#"0.0007 0.0007 0.0007 0.0007 0.0007 0.0007"
    params["NDutHits"] = "1"

def daffitter(params):
    params["TemplateFile"] = "daf-fitter-pre-tmp.xml"
    dafcommon(params)

def dafmille(params):
    params["TemplateFile"] = "daf-fitter-mille-tmp.xml"
    dafcommon(params)

if __name__ == "__main__":
    from steeringGenerator import *
    functions = {"noise": noise,
                 "rawtohit": rawtohit,
                 "prealign": prealign,
                 "aligntel": aligntel,
                 "aligndut": aligndut,
                 "millealign": millealign,
                 "dafmille": dafmille,
                 "daffitter": daffitter}
    jobMaker(functions, runConfig)

