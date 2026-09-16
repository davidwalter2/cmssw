// WMass: LHE weight tables in the grouped layout of the 10_6 custom NanoAOD.
//
// One FlatTable per weight group of the LHE header, named <Type>[AltSetN], with
// N counting the groups of that type in header order:
//   LHEScaleWeight, LHEScaleWeightAltSet1, ...
//   LHEPdfWeight,   LHEPdfWeightAltSet1,   ...
//   MEParamWeight,  MEParamWeightAltSet1,  ...
//   UnknownWeight,  UnknownWeightAltSet1,  ...
// The stock GenWeightsTableProducer keeps genWeight, PSWeight, the named and the
// LHEReweightingWeight tables; its own LHEScale/LHEPdf products overlap with the
// tables here and are dropped from the output by the WMass customise.
//
// The classification follows the 10_6 GeneratorInterface/Core WeightHelper
// (never merged upstream), decided on the group name:
//   scale    : contains "scale_variation" or "Central scale variation"
//   PDF      : contains "PDF_variation", or is an LHAPDF set name
//   ME param : contains "mg_reweighting" or "variation"
//              (mass_variation, width_variation, sthw2_variation, q0_variation, ...)
//   PS       : contains isr / fsr / nominal / baseline / emission
//   unknown  : anything else, including weights outside any <weightgroup>
//              (each of those is its own group, as in 10_6)
// A PDF group is split into one table per LHAPDF *set*: the set of a weight is the
// index entry with the largest LHAID <= the weight's LHAID (LHAPDF::lookupPDF), and
// a new table starts whenever that set changes along the header. The index is a
// parameter (pdfSetsIndex) rather than the release's LHAPDF: the 10_6 production
// ran with LHAPDF 6.2.1, whose index predates NNPDF4.0, CT18, MSHT20, PDF4LHC21...,
// so those blocks were lumped into the preceding known set, and the analysis code
// addresses the AltSetN tables by number. Reproducing the numbering means
// reproducing that index; the shipped default is the 6.2.1 index of CMSSW_10_6_26.
//
// Weights are w / originalXWGTUP, the scale group is reordered to the 10_6 / stock
// (muR, muF) order when all nine variations are present.

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/transform.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <regex>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

  enum class WType { Scale = 0, Pdf, MEParam, PS, Unknown };
  const std::vector<std::string> kTypeTable = {"LHEScaleWeight", "LHEPdfWeight", "MEParamWeight", "PSWeight", "UnknownWeight"};

  // the 10_6 weightgroups strings ("scale", "PDF", "matrix element", "unknown", "parton shower"):
  // only the first character is read, and the capitalisation matters
  WType typeFromString(const std::string& s) {
    if (s.empty())
      throw cms::Exception("Configuration", "empty weight group type");
    switch (s[0]) {
      case 's':
        return WType::Scale;
      case 'P':
        return WType::Pdf;
      case 'm':
        return WType::MEParam;
      case 'p':
        return WType::PS;
      case 'u':
        return WType::Unknown;
    }
    throw cms::Exception("Configuration", "unknown weight group type '" + s + "'");
  }

  struct ParsedWeight {
    std::string id;
    std::string content;                                   // text between <weight ...> and </weight>
    std::unordered_map<std::string, std::string> attrs;   // tag attributes
  };

  struct Group {
    WType type = WType::Unknown;
    std::string name;
    std::string doc;
    std::vector<std::string> ids;  // weight ids, in the order the values are written
    // PDF
    int parentLhaid = -1;
    // scale: (muR, muF) -> id, for the reordering
    std::map<std::pair<float, float>, std::string> scaleIds;
    bool scaleHasDyn = false;
  };

  struct RunGroups {
    std::vector<Group> groups;  // the groups to store, in table order
    std::vector<std::string> names;
  };

  // LHAPDF's pdfsets.index: "LHAID SetName Version"
  class PdfSetsIndex {
  public:
    explicit PdfSetsIndex(const std::string& path) {
      std::ifstream in(path);
      if (!in)
        throw cms::Exception("Configuration", "cannot open pdfSetsIndex " + path);
      std::string line;
      while (std::getline(in, line)) {
        std::istringstream ss(line);
        int id;
        std::string name;
        if (ss >> id >> name) {
          byId_[id] = name;
          byName_.emplace(name, id);
        }
      }
      if (byId_.empty())
        throw cms::Exception("Configuration", "no entries in pdfSetsIndex " + path);
    }
    // LHAPDF::lookupPDF(int): (set name, member) of the index entry with the largest LHAID <= lhaid
    std::pair<std::string, int> lookupPDF(int lhaid) const {
      auto it = byId_.upper_bound(lhaid);
      if (it == byId_.begin())
        return {"", -1};
      --it;
      return {it->second, lhaid - it->first};
    }
    // LHAPDF::lookupLHAPDFID(name)
    int lookupLHAPDFID(const std::string& name) const {
      auto it = byName_.find(name);
      return it == byName_.end() ? -1 : it->second;
    }
    int parentLhaid(int lhaid) const { return lhaid - lookupPDF(lhaid).second; }

  private:
    std::map<int, std::string> byId_;
    std::unordered_map<std::string, int> byName_;
  };

  // the 10_6 WeightHelper attribute aliases
  const std::unordered_map<std::string, std::vector<std::string>> kAttributeNames = {
      {"muf", {"muF", "MUF", "muf", "facscfact"}},
      {"mur", {"muR", "MUR", "mur", "renscfact"}},
      {"pdf", {"PDF", "PDF set", "lhapdf", "pdf", "pdf set", "pdfset"}},
      {"dyn", {"DYN_SCALE"}},
  };

  std::string trim(const std::string& s) {
    auto b = s.find_first_not_of(" \t\r\n\"");
    if (b == std::string::npos)
      return "";
    auto e = s.find_last_not_of(" \t\r\n\"");
    return s.substr(b, e - b + 1);
  }

  // WeightHelper::searchAttributes: by tag attribute first, then by "label = value" in the content
  std::string searchAttribute(const std::string& label, const ParsedWeight& w) {
    const auto& names = kAttributeNames.at(label);
    for (const auto& n : names) {
      auto it = w.attrs.find(n);
      if (it != w.attrs.end())
        return trim(it->second);
    }
    std::smatch m;
    for (const auto& n : names) {
      std::regex floatExpr(n + "\\s*=\\s*([0-9.]+(?:[eE][+-]?[0-9]+)?)");
      std::regex strExpr(n + "\\s*=\\s*([^=]+)");
      if (std::regex_search(w.content, m, floatExpr))
        return trim(m.str(1));
      if (std::regex_search(w.content, m, strExpr))
        return trim(m.str(1));
    }
    return "";
  }

  // "1d0", "0.5d0" (fortran) and plain floats; nan when not a number
  float toFloat(const std::string& s) {
    if (s.empty())
      return std::nanf("");
    std::string t = s;
    auto d = t.find_first_of("dD");
    if (d != std::string::npos)
      t[d] = 'e';
    try {
      size_t used = 0;
      float v = std::stof(t, &used);
      return used > 0 ? v : std::nanf("");
    } catch (...) {
      return std::nanf("");
    }
  }

  bool containsCI(const std::string& s, const std::string& what) {
    std::string l = s;
    std::transform(l.begin(), l.end(), l.begin(), ::tolower);
    return l.find(what) != std::string::npos;
  }

  const std::vector<std::pair<float, float>> kScaleOrder = {
      {0.5, 0.5}, {0.5, 1}, {0.5, 2}, {1, 0.5}, {1, 1}, {1, 2}, {2, 0.5}, {2, 1}, {2, 2}};

}  // namespace

class LHEWeightGroupsTableProducer : public edm::global::EDProducer<edm::RunCache<RunGroups>> {
public:
  LHEWeightGroupsTableProducer(edm::ParameterSet const& params)
      : lheLabel_(params.getParameter<std::vector<edm::InputTag>>("lheInfo")),
        lheTag_(edm::vector_transform(lheLabel_,
                                      [this](const edm::InputTag& tag) { return mayConsume<LHEEventProduct>(tag); })),
        lheRunTag_(edm::vector_transform(
            lheLabel_, [this](const edm::InputTag& tag) { return mayConsume<LHERunInfoProduct, edm::InRun>(tag); })),
        precision_(params.getParameter<int32_t>("lheWeightPrecision")),
        weightgroups_(edm::vector_transform(params.getParameter<std::vector<std::string>>("weightgroups"),
                                            [](const std::string& s) { return typeFromString(s); })),
        maxGroupsPerType_(params.getParameter<std::vector<int32_t>>("maxGroupsPerType")),
        unknownOnlyIfEmpty_(edm::vector_transform(params.getParameter<std::vector<std::string>>("unknownOnlyIfEmpty"),
                                                  [](const std::string& s) { return typeFromString(s); })),
        pdfIndex_(params.getParameter<edm::FileInPath>("pdfSetsIndex").fullPath()) {
    if (weightgroups_.size() != maxGroupsPerType_.size())
      throw cms::Exception("Configuration", "'weightgroups' and 'maxGroupsPerType' must have equal size");
    produces<std::vector<nanoaod::FlatTable>>("LHEWeightTableVec");
  }

  ~LHEWeightGroupsTableProducer() override {}

  // ---------------------------------------------------------------- run: parse the header

  std::shared_ptr<RunGroups> globalBeginRun(edm::Run const& iRun, edm::EventSetup const&) const override {
    auto out = std::make_shared<RunGroups>();
    // LHERunInfoProduct is put at endRun by externalLHEProducer, so a getByToken here
    // throws ("read a Run product before endRun()", cms-sw#18499); the legacy getByLabel
    // still serves it from the input file (with a LogicError warning) -- the same
    // access the stock GenWeightsTableProducer and the 10_6 LHEWeightProductProducer use.
    edm::Handle<LHERunInfoProduct> lheInfo;
    for (const auto& label : lheLabel_) {
      iRun.getByLabel(label, lheInfo);
      if (lheInfo.isValid())
        break;
    }
    if (!lheInfo.isValid())
      return out;

    std::string header;
    for (auto it = lheInfo->headers_begin(), end = lheInfo->headers_end(); it != end; ++it) {
      if (it->tag() != "initrwgt")
        continue;
      for (const auto& line : it->lines())
        header += line + "\n";
      break;
    }
    if (header.empty())
      return out;
    // some samples write the header HTML-escaped
    const std::vector<std::pair<std::string, std::string>> reps = {{"&lt;", "<"}, {"&gt;", ">"}};
    for (const auto& rep : reps) {
      size_t pos = 0;
      while ((pos = header.find(rep.first, pos)) != std::string::npos) {
        header.replace(pos, rep.first.size(), rep.second);
        pos += rep.second.size();
      }
    }

    std::vector<Group> all = buildGroups(header);

    // groupsToStore: per requested type, the groups of that type in header order (up to the max);
    // unknown groups only if none of the unknownOnlyIfEmpty types was found
    bool storeUnknown = unknownOnlyIfEmpty_.empty();
    for (auto t : unknownOnlyIfEmpty_) {
      if (std::none_of(all.begin(), all.end(), [t](const Group& g) { return g.type == t; }))
        storeUnknown = true;
    }
    for (size_t i = 0; i < weightgroups_.size(); ++i) {
      WType t = weightgroups_[i];
      int maxStore = maxGroupsPerType_[i];
      if (maxStore == 0 || (t == WType::Unknown && !storeUnknown))
        continue;
      int n = 0;
      for (const auto& g : all) {
        if (g.type != t)
          continue;
        if (maxStore > 0 && n >= maxStore)
          break;
        std::string name = kTypeTable[static_cast<int>(t)];
        if (n > 0)
          name += "AltSet" + std::to_string(n);
        out->groups.push_back(g);
        out->names.push_back(name);
        ++n;
      }
    }

    std::ostringstream summary;
    for (size_t i = 0; i < out->groups.size(); ++i)
      summary << "\n  " << out->names[i] << " <- '" << out->groups[i].name << "' (" << out->groups[i].ids.size()
              << " weights)";
    edm::LogInfo("LHEWeightGroupsTableProducer") << "LHE weight groups of run " << iRun.run() << ":" << summary.str();
    return out;
  }

  void globalEndRun(edm::Run const&, edm::EventSetup const&) const override {}

  // ---------------------------------------------------------------- event: fill the tables

  void produce(edm::StreamID id, edm::Event& iEvent, const edm::EventSetup& iSetup) const override {
    auto tables = std::make_unique<std::vector<nanoaod::FlatTable>>();
    const RunGroups* rg = runCache(iEvent.getRun().index());

    edm::Handle<LHEEventProduct> lheInfo;
    for (const auto& tag : lheTag_) {
      iEvent.getByToken(tag, lheInfo);
      if (lheInfo.isValid())
        break;
    }
    if (rg && lheInfo.isValid() && !rg->groups.empty()) {
      std::unordered_map<std::string, double> byId;
      byId.reserve(lheInfo->weights().size());
      for (const auto& w : lheInfo->weights())
        byId[w.id] = w.wgt;
      const double w0 = lheInfo->originalXWGTUP();

      tables->reserve(rg->groups.size());
      for (size_t ig = 0; ig < rg->groups.size(); ++ig) {
        const Group& g = rg->groups[ig];
        std::vector<float> vals;
        vals.reserve(g.ids.size());
        for (const auto& wid : g.ids) {
          auto it = byId.find(wid);
          if (it == byId.end())
            throw cms::Exception("LHEWeightGroupsTableProducer")
                << "weight id '" << wid << "' of group '" << g.name << "' (table " << rg->names[ig]
                << ") listed in the LHE header is missing from the event weights";
          vals.push_back(it->second / w0);
        }
        tables->emplace_back(vals.size(), rg->names[ig], false, false);
        tables->back().addColumn<float>("", vals, g.doc, precision_);
      }
    }
    iEvent.put(std::move(tables), "LHEWeightTableVec");
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<std::vector<edm::InputTag>>("lheInfo", std::vector<edm::InputTag>{{"externalLHEProducer"}, {"source"}})
        ->setComment("tag(s) for the LHE information (LHEEventProduct and LHERunInfoProduct)");
    desc.add<int32_t>("lheWeightPrecision", 14)->setComment("mantissa bits of the stored weights");
    desc.add<std::vector<std::string>>("weightgroups", {"scale", "PDF", "matrix element", "unknown"})
        ->setComment(
            "weight group types to store, in table order; only the first character is read: "
            "'scale', 'PDF', 'matrix element', 'unknown', 'parton shower'");
    desc.add<std::vector<int32_t>>("maxGroupsPerType", {-1, -1, -1, -1})
        ->setComment("max number of groups per type above, -1 = all found");
    desc.add<std::vector<std::string>>("unknownOnlyIfEmpty", {"scale", "PDF"})
        ->setComment("store unknown groups only if one of these types has no group");
    desc.add<edm::FileInPath>("pdfSetsIndex",
                              edm::FileInPath("PhysicsTools/NanoAOD/data/pdfsets_lhapdf-6.2.1-pafccj3.index"))
        ->setComment(
            "LHAPDF pdfsets.index used to split PDF groups into sets (defines the AltSetN numbering); "
            "default = the LHAPDF 6.2.1 index of CMSSW_10_6_26, for parity with the 10_6 custom NanoAOD");
    descriptions.add("lheWeightGroupsTable", desc);
  }

private:
  // ---------------------------------------------------------------- header parsing

  static std::string groupNameOf(const std::string& tag) {
    // <weightgroup name="..."> or <weightgroup type="...">; text after a '.' is dropped (10_6 parseGroupName)
    std::smatch m;
    for (const char* att : {"name", "type"}) {
      std::regex r(std::string(att) + "\\s*=\\s*\"([^\"]*)\"");
      if (std::regex_search(tag, m, r)) {
        std::string name = m.str(1);
        auto dot = name.find('.');
        if (dot != std::string::npos)
          name.erase(dot);
        return name;
      }
    }
    return "";
  }

  static ParsedWeight parseWeight(const std::string& tag, const std::string& content) {
    ParsedWeight w;
    w.content = content;
    std::regex attr("([A-Za-z_][A-Za-z0-9_]*)\\s*=\\s*\"([^\"]*)\"");
    for (auto it = std::sregex_iterator(tag.begin(), tag.end(), attr); it != std::sregex_iterator(); ++it)
      w.attrs[(*it).str(1)] = (*it).str(2);
    // also id='...' (single quotes, some MG headers)
    std::smatch m;
    if (w.attrs.find("id") == w.attrs.end() && std::regex_search(tag, m, std::regex("id\\s*=\\s*'([^']*)'")))
      w.attrs["id"] = m.str(1);
    auto it = w.attrs.find("id");
    w.id = it == w.attrs.end() ? "" : it->second;
    return w;
  }

  int lhapdfId(const ParsedWeight& w, const std::string& groupname) const {
    std::string text = searchAttribute("pdf", w);
    if (!text.empty()) {
      try {
        return std::stoi(text);
      } catch (...) {
      }
    }
    return pdfIndex_.lookupLHAPDFID(groupname);
  }

  // WeightHelper::buildGroup: the type is decided on the group name (and, for an orphan
  // PDF weight, on the pdf attribute of the first weight)
  WType classify(const std::string& name, const ParsedWeight& first, std::string& nameOut) const {
    nameOut = name;
    if (name.find("scale_variation") != std::string::npos || name.find("Central scale variation") != std::string::npos)
      return WType::Scale;
    if (name.find("PDF_variation") != std::string::npos || pdfIndex_.lookupLHAPDFID(name) != -1)
      return WType::Pdf;
    if (name.find("mg_reweighting") != std::string::npos || name.find("variation") != std::string::npos)
      return WType::MEParam;
    for (const char* ps : {"isr", "fsr", "nominal", "baseline", "emission"})
      if (containsCI(name, ps))
        return WType::PS;
    // orphan PDF weight: a pdf attribute pointing at the first member of a set
    std::string text = searchAttribute("pdf", first);
    if (!text.empty()) {
      try {
        auto set = pdfIndex_.lookupPDF(std::stoi(text));
        if (!set.first.empty() && set.second == 0) {
          nameOut = set.first;
          return WType::Pdf;
        }
      } catch (...) {
      }
    }
    return WType::Unknown;
  }

  void addToGroup(Group& g, const ParsedWeight& w) const {
    g.ids.push_back(w.id);
    if (g.type == WType::Scale) {
      float muR = toFloat(searchAttribute("mur", w));
      float muF = toFloat(searchAttribute("muf", w));
      if (!searchAttribute("dyn", w).empty())
        g.scaleHasDyn = true;
      else if (!std::isnan(muR) && !std::isnan(muF))
        g.scaleIds[{muR, muF}] = w.id;
    }
  }

  void startPdfGroup(Group& g, int lhaid) const {
    g.parentLhaid = pdfIndex_.parentLhaid(lhaid);
    g.doc = "Uncertainty sets for LHAPDF set " + pdfIndex_.lookupPDF(g.parentLhaid).first +
            " with LHAID = " + std::to_string(g.parentLhaid) + "; ";
  }

  void finishGroup(Group& g) const {
    std::ostringstream doc;
    if (g.type == WType::Scale) {
      // 10_6 orderedScaleWeights / stock: (muR, muF) ascending when all nine are there
      doc << "LHE scale variation weights (w_var / w_nominal); ";
      bool wellFormed = !g.scaleHasDyn;
      for (const auto& s : kScaleOrder)
        wellFormed = wellFormed && g.scaleIds.count(s);
      if (wellFormed) {
        std::vector<std::string> ordered;
        for (size_t i = 0; i < kScaleOrder.size(); ++i) {
          ordered.push_back(g.scaleIds.at(kScaleOrder[i]));
          doc << "[" << i << "] is muR=" << kScaleOrder[i].first << " muF=" << kScaleOrder[i].second
              << (i + 1 < kScaleOrder.size() ? "; " : "");
        }
        g.ids = ordered;
      } else {
        size_t nstore = std::min<size_t>(kScaleOrder.size(), g.ids.size());
        g.ids.resize(nstore);
        doc << "WARNING: Unexpected format found. Contains first " << nstore << " elements of weights vector, unordered";
      }
      g.doc = doc.str();
    } else if (g.type == WType::Pdf) {
      g.doc += std::to_string(g.ids.size()) + " weights";
    } else {
      // ME parameter / unknown: the group name and the content of every weight
      // (e.g. "mass=91.0876"), so the variation grid is documented in the file
      doc << (g.name.empty() ? std::string("(no weight group)") : g.name) << "; " << g.ids.size() << " weights";
      g.doc = doc.str();
    }
  }

  std::vector<Group> buildGroups(const std::string& header) const {
    std::vector<Group> groups;
    // the docs of ME-param / unknown groups carry the per-weight content
    std::vector<std::vector<std::string>> contents;

    std::regex token("<weightgroup\\b[^>]*>|</weightgroup>|<weight\\b[^>]*>[\\s\\S]*?</weight>");
    bool inGroup = false;
    std::string groupName;
    bool groupOpen = false;  // a Group has been created for the current <weightgroup>
    auto newGroup = [&](const std::string& name, const ParsedWeight& first) {
      Group g;
      g.type = classify(name, first, g.name);
      groups.push_back(g);
      contents.emplace_back();
    };
    for (auto it = std::sregex_iterator(header.begin(), header.end(), token); it != std::sregex_iterator(); ++it) {
      const std::string tok = (*it).str();
      if (tok.compare(0, 12, "<weightgroup") == 0) {
        inGroup = true;
        groupOpen = false;
        groupName = groupNameOf(tok);
      } else if (tok.compare(0, 13, "</weightgroup") == 0) {
        inGroup = false;
        groupOpen = false;
      } else {
        auto gt = tok.find('>');
        auto end = tok.rfind("</weight>");
        ParsedWeight w = parseWeight(tok.substr(0, gt + 1), tok.substr(gt + 1, end - gt - 1));
        if (w.id.empty())
          continue;
        if (!inGroup) {
          // a weight outside any group is a group of its own (10_6 parseWeights)
          newGroup("", w);
          groupOpen = false;
        } else if (!groupOpen) {
          newGroup(groupName, w);
          groupOpen = true;
        }
        Group& g = groups.back();
        if (g.type == WType::Pdf) {
          // WeightHelper::splitPdfWeight + updatePdfInfo: a new set starts a new group
          int lhaid = lhapdfId(w, g.name);
          if (g.parentLhaid < 0) {
            if (lhaid > 0)
              startPdfGroup(g, lhaid);
          } else if (lhaid > 0 && pdfIndex_.parentLhaid(lhaid) != g.parentLhaid) {
            Group split;
            split.type = WType::Pdf;
            split.name = g.name;
            groups.push_back(split);
            contents.emplace_back();
            startPdfGroup(groups.back(), lhaid);
          }
        }
        addToGroup(groups.back(), w);
        contents.back().push_back(trim(w.content));
      }
    }
    for (size_t i = 0; i < groups.size(); ++i) {
      finishGroup(groups[i]);
      if (groups[i].type == WType::MEParam || groups[i].type == WType::Unknown) {
        std::ostringstream doc;
        doc << groups[i].doc;
        size_t shown = 0;
        for (size_t k = 0; k < contents[i].size() && doc.tellp() < 1800; ++k, ++shown)
          doc << "; [" << k << "] " << contents[i][k].substr(0, 40);
        if (shown < contents[i].size())
          doc << "; ...";
        groups[i].doc = doc.str();
      }
    }
    return groups;
  }

  const std::vector<edm::InputTag> lheLabel_;
  const std::vector<edm::EDGetTokenT<LHEEventProduct>> lheTag_;
  const std::vector<edm::EDGetTokenT<LHERunInfoProduct>> lheRunTag_;
  const int32_t precision_;
  const std::vector<WType> weightgroups_;
  const std::vector<int32_t> maxGroupsPerType_;
  const std::vector<WType> unknownOnlyIfEmpty_;
  const PdfSetsIndex pdfIndex_;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(LHEWeightGroupsTableProducer);
