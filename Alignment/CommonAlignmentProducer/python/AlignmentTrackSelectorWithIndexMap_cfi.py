import FWCore.ParameterSet.Config as cms
from Alignment.CommonAlignmentProducer.AlignmentTrackSelector_cfi import AlignmentTrackSelector

# Same parameters as AlignmentTrackSelector_cfi.AlignmentTrackSelector but
# bound to the AlignmentTrackSelectorWithIndexMapModule C++ class, which
# additionally emits an edm::ValueMap<unsigned int> (instance label
# "originalIndex") of source-track indices alongside the cloned tracks.
# Use this in any chain where downstream consumers (DeDxValueMapProjector,
# VertexCompositeCandidateRemapper) need to navigate from a cloned track
# back to its source TrackCollection entry.
AlignmentTrackSelectorWithIndexMap = cms.EDFilter(
    'AlignmentTrackSelectorWithIndexMapModule',
    **{name: getattr(AlignmentTrackSelector, name) for name in AlignmentTrackSelector.parameterNames_()}
)
