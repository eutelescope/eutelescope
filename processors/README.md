Description of the Processors and parameter options.

For May 2018: Focus on the common processors using GBL:
* EUTelAlignGBL (for alignment)
* EUTelGBLFitter (for track fitting)

# EUTelAlignGBL

Parameters (required):
* "Ebeam": "Beam energy [GeV]", _eBeam, static_cast<double>(4.0)
* "IsFirstAlignStep": "Bool: 1/0 (yes/no)", _IsFirstAlignStep, static_cast<int>(0)
* "kappa", "Global factor to Highland formula", _kappa, static_cast <double>(1.0), 1.0 means HL as is, 1.2 means 20% additional scattering

Parameters (optional):
* "ExcludePlanes", "Exclude planes from fit according to their sensor ids.", _excludePlanes_sensorIDs ,std::vector<int>());
* "FixedPlanes","Fix sensor planes in the fit according to their sensor ids.",_FixedPlanes_sensorIDs ,std::vector<int>()
* "MaxTrackCandidatesTotal","Maximal number of track candidates (Total).",_maxTrackCandidatesTotal, static_cast <int> (10000000)
* "MaxTrackCandidates","Maximal number of track candidates.",_maxTrackCandidates, static_cast <int> (2000)
* "BinaryFilename","Name of the Millepede binary file.",_binaryFilename, string ("mille.bin")
* "TelescopeResolution","Resolution of the telescope for Millepede (sigma_x=sigma_y.",_telescopeResolution, static_cast <float> (0.010)
* "AlignMode","Number of alignment constants used. Available mode are: "
                              "\nXYZShifts - shifts in X and Y"
                              "\nXYShiftsRotZ - shifts in X and Y and rotation around the Z axis,"
                              "\nXYZShiftsRotZ - shifts in X,Y and Z and rotation around the Z axis",
                              _alignModeString, std::string("XYShiftsRotZ")
* "triCut", "Upstream triplet residual cut [um]", _triCut, 0.30
* "driCut", "Downstream triplet residual cut [um]", _driCut, 0.40
* "sixCut", "Upstream-Downstream Track matching cut [um]", _sixCut, 0.60
* "slopeCut", "t(d)riplet slope cut [radian]", _slopeCut, 0.01
* "GeneratePedeSteerfile","Generate a steering file for the pede program.",_generatePedeSteerfile, static_cast <int> (0)
* "PedeSteerfileName","Name of the steering file for the pede program.",_pedeSteerfileName, string("steer_mille.txt")

# EUTelGBLFitter

tbd
