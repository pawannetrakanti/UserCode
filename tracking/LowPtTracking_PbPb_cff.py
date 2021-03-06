import FWCore.ParameterSet.Config as cms

from RecoHI.HiTracking.HighPtTracking_PbPb_cff import *

#Commented by Pawan  (Default)
#hiPixel3PrimTracks.RegionFactoryPSet.RegionPSet.ptMin = 0.9
#hiPixel3PrimTracks.RegionFactoryPSet.RegionPSet.originRadius = 0.1
#hiPixel3PrimTracks.FilterPSet.ptMin = 0.9
#ckfBaseTrajectoryFilter.filterPset.minPt = 0.9

#Lowered pT threshold
hiPixel3PrimTracks.RegionFactoryPSet.RegionPSet.ptMin = 0.4
hiPixel3PrimTracks.RegionFactoryPSet.RegionPSet.originRadius = 0.1
hiPixel3PrimTracks.FilterPSet.ptMin = 0.4
ckfBaseTrajectoryFilter.filterPset.minPt = 0.4


