import FWCore.ParameterSet.Config as cms

# WMass: LHE weight tables in the grouped layout of the 10_6 custom NanoAOD
# (LHEScaleWeight[AltSetN], LHEPdfWeight[AltSetN], MEParamWeight[AltSetN],
# UnknownWeight[AltSetN]); see plugins/LHEWeightGroupsTableProducer.cc. The
# pdfSetsIndex fixes the AltSetN numbering to the one of the 10_6 production
# (LHAPDF 6.2.1 index), which the analysis code addresses by number.
lheWeightGroupsTable = cms.EDProducer("LHEWeightGroupsTableProducer",
    lheInfo = cms.VInputTag(cms.InputTag("externalLHEProducer"), cms.InputTag("source")),
    lheWeightPrecision = cms.int32(14),
    # Only the first character is read, and the capitalisation matters
    # ('scale', 'PDF', 'matrix element', 'unknown', 'parton shower'); table order.
    weightgroups = cms.vstring('scale', 'PDF', 'matrix element', 'unknown'),
    # Max number of groups to store per type above, -1 = all found
    maxGroupsPerType = cms.vint32(-1, -1, -1, -1),
    # Unknown groups only if one of these types has no group at all
    unknownOnlyIfEmpty = cms.vstring('scale', 'PDF'),
    pdfSetsIndex = cms.FileInPath("PhysicsTools/NanoAOD/data/pdfsets_lhapdf-6.2.1-pafccj3.index"),
)

lheWeightGroupsTableTask = cms.Task(lheWeightGroupsTable)


def nanoAOD_addLHEWeightGroups(process):
    """Schedule the grouped LHE weight tables with the gen weight tables and
    drop the stock LHEScale/LHEPdf products (same content, other layout) from
    the NanoAOD output. The vector-of-tables product is kept through
    NanoAODEDMEventContent_cff ("keep nanoaodFlatTables_*Table_*_*")."""
    process.lheWeightGroupsTable = lheWeightGroupsTable
    process.genWeightsTableTask.add(process.lheWeightGroupsTable)
    for name, out in process.outputModules.items():
        if out.type_() == "NanoAODOutputModule":
            out.outputCommands.extend([
                "drop nanoaodFlatTable_genWeightsTable_LHEScale_*",
                "drop nanoaodFlatTable_genWeightsTable_LHEPdf_*",
            ])
    return process
