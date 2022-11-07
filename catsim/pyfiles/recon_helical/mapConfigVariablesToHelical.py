from math import ceil
from copy import deepcopy

def mapConfigVariablesToHelical(cfg):

    #sid = 10  # The distance from source to object (cm)
    #sdd = 15  # The distance from detector to object center(cm)
    #YL = 350  # Detector cell number along the horizontal direction of detector array
    #ZL = 20  # Detector cell number along the vertical direction of detector array
    sid = cfg.scanner.sid
    sdd = cfg.scanner.sdd
    YL = int(cfg.scanner.detectorColCount)
    ZL = int(cfg.scanner.detectorRowCount)

    #N_Turn = 4  # The number of turns for the whole helical scan
    #N_2pi = 360  # The projections/views number for each turn of scan
    #h = 0.1  # Helical pitch related to detector height added by Jiayong: this definition is different than what is xcist defined
    N_2pi = cfg.protocol.viewsPerRotation # total numbers of view per rotation
    N_Turn = cfg.protocol.viewCount/cfg.protocol.viewsPerRotation
    #h = cfg.protocol.tableSpeed
    h = 16.

    #dectorYoffset = 0
    #dectorZoffset = 0
    dectorYoffset = -cfg.scanner.detectorColOffset
    dectorZoffset = cfg.scanner.detectorRowOffset

    # question: what values should I choose for these values?
    # The following lines are used to define the reconstruction paramters
    k1 = 5  # The order to define the 3D weighting function
    delta = 60  # The range to define smoothness of 2D weigthing function
    HSCoef = 0.6  # This is used to define the half-scan range

    #nMod=1
    #rowSize=1
    #modWidth=1
    nMod = ceil(cfg.scanner.detectorColCount/cfg.scanner.detectorColsPerMod)
    rowSize = cfg.scanner.detectorRowSize
    modWidth = cfg.scanner.detectorColsPerMod*cfg.scanner.detectorColSize

    #imageSize = 128  # Define the size of the reconstructed image
    #sliceCount = 1
    #sliceThickness = 1
    #fov  = 1.1  # Diameter of imaginh object
    #kernelType = 'SL'
    #centerOffset = [0,0,0]

    imageSize = cfg.recon.imageSize
    sliceCount = cfg.recon.sliceCount
    sliceThickness = cfg.recon.sliceThickness
    fov = cfg.recon.fov
    kernelType = cfg.recon.kernelType
    centerOffset = deepcopy(cfg.recon.centerOffset)
    # Pass desired X as Y
    centerOffset[1] = deepcopy(cfg.recon.centerOffset[0])
    # Pass desired Y as -X
    centerOffset[0] = -deepcopy(cfg.recon.centerOffset[1])

    return  sid, sdd, YL, ZL, N_Turn, N_2pi, h, dectorYoffset, dectorZoffset, \
            k1, delta, HSCoef,nMod, rowSize, modWidth, imageSize,sliceCount, sliceThickness, centerOffset, fov, kernelType
