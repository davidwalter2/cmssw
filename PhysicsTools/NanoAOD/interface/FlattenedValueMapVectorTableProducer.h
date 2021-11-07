#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include <vector>
#include <iostream>
#include <boost/ptr_container/ptr_vector.hpp>

template<typename T>
class FlattenedValueMapVectorTableProducer : public edm::stream::EDProducer<> {
    public:

        FlattenedValueMapVectorTableProducer( edm::ParameterSet const & params ):
            name_( params.getParameter<std::string>("name") ),
            doc_(params.existsAs<std::string>("doc") ? params.getParameter<std::string>("doc") : ""),
            src_(consumes<edm::View<T>>(params.getParameter<edm::InputTag>("src"))),
            cut_(params.existsAs<std::string>("cut") ? params.getParameter<std::string>("cut") : "", true),
            countPrecision_(params.existsAs<int>("countPrecision") ? params.getParameter<int>("countPrecision") : -1) 
        {
            edm::ParameterSet const & varsPSet = params.getParameter<edm::ParameterSet>("variables");
            for (const std::string & vname : varsPSet.getParameterNamesForType<edm::ParameterSet>()) {
                const auto & varPSet = varsPSet.getParameter<edm::ParameterSet>(vname);
                const std::string & type = varPSet.getParameter<std::string>("type");
                if (type == "std::vector<int>") {
                    intVecMaps_.emplace_back(consumes<edm::ValueMap<std::vector<int>>>(varPSet.getParameter<edm::InputTag>("src")));
                    intNames_.emplace_back(vname);
                    intDocs_.emplace_back(varPSet.getParameter<std::string>("doc"));
                    intPrecisions_.emplace_back(varPSet.getParameter<int>("precision"));
                }
                else if (type == "std::vector<float>") {
                    floatVecMaps_.emplace_back(consumes<edm::ValueMap<std::vector<float>>>(varPSet.getParameter<edm::InputTag>("src")));
                    floatNames_.emplace_back(vname);
                    floatDocs_.emplace_back(varPSet.getParameter<std::string>("doc"));
                    floatPrecisions_.emplace_back(varPSet.getParameter<int>("precision"));
                }
                else throw cms::Exception("Configuration", "unsupported type "+type+" for variable "+vname);
            }

            produces<nanoaod::FlatTable>("floatcounts");
            produces<nanoaod::FlatTable>("intcounts");
            produces<nanoaod::FlatTable>("intvecs");
            produces<nanoaod::FlatTable>("floatvecs");
        }

        ~FlattenedValueMapVectorTableProducer() override {}

        template<typename P>
        std::vector<P> readVals(const edm::ValueMap<std::vector<P>>& vmap, edm::PtrVector<T>& objs, std::vector<int>& sizes) {
            std::vector<P> allvals;
            // Will be at least this long
            allvals.reserve(objs.size());
            for (size_t i = 0; i < objs.size(); i++) {
                auto& vals = vmap[objs[i]];
                if (sizes[i] == 0)
                    sizes[i] = vals.size();
                else if (sizes[i] != static_cast<int>(vals.size())) {
                    throw cms::Exception("Unmatched values", std::string("All vector value maps in the table need to have the same length. ") +
                        std::string("If you want to add vectors of a different length, make a new table call\n ") +
                        std::string("Size of initial collection is ") + std::to_string(sizes[i]) + 
                        " tried to add collection of size " + std::to_string(vals.size()) + std::string(" at index ") + std::to_string(i) +
                        std::string(". Number of objects is ") + std::to_string(objs.size()));
                }
                allvals.insert(std::end(allvals), std::begin(vals), std::end(vals));
            }

            return allvals;
        }

        void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override {
            edm::Handle<edm::View<T>> src;
            iEvent.getByToken(src_, src);

            edm::PtrVector<T> objs;
            for (size_t i = 0; i < src->size(); ++i) {
                edm::Ptr<T> obj = src->ptrAt(i);
                if (cut_(*obj)) { 
                    objs.push_back(obj);
                }
            }

            std::vector<int> sizes(objs.size(), 0);
            auto intsizetab = std::make_unique<nanoaod::FlatTable>(objs.size(), this->name_, false, false);
            auto floatsizetab = std::make_unique<nanoaod::FlatTable>(objs.size(), this->name_, false, true);
            std::unique_ptr<nanoaod::FlatTable> floatvectab; ;
            std::unique_ptr<nanoaod::FlatTable> intvectab; ;

            for (size_t i = 0; i < intVecMaps_.size(); i++) {
                edm::Handle<edm::ValueMap<std::vector<int>>> vmap;
                iEvent.getByToken(intVecMaps_[i], vmap);
                const auto& results = readVals(*vmap, objs, sizes);
                if (i == 0)
                    intvectab = std::make_unique<nanoaod::FlatTable>(results.size(), this->name_+"_IntVals", false, false);
                intvectab->addColumn<int>(intNames_[i], results, intDocs_[i], nanoaod::FlatTable::IntColumn, intPrecisions_[i]);
            }
            intsizetab->template addColumn<int>("IntCounts", sizes, "Number of entries per object", nanoaod::FlatTable::IntColumn, countPrecision_);
            std::fill(sizes.begin(), sizes.end(), 0);
            for (size_t i = 0; i < floatVecMaps_.size(); i++) {
                edm::Handle<edm::ValueMap<std::vector<float>>> vmap;
                iEvent.getByToken(floatVecMaps_[i], vmap);
                const auto& results = readVals(*vmap, objs, sizes);
                if (i == 0)
                    floatvectab = std::make_unique<nanoaod::FlatTable>(results.size(), this->name_+"_FloatVals", false, false);
                floatvectab->addColumn<float>(floatNames_[i], results, floatDocs_[i], nanoaod::FlatTable::FloatColumn, floatPrecisions_[i]);
            }
            floatsizetab->template addColumn<int>("FloatCounts", sizes, "Number of entries per object", nanoaod::FlatTable::IntColumn, countPrecision_);


            intsizetab->setDoc(doc_);
            floatsizetab->setDoc(doc_);
            floatvectab->setDoc(doc_);
            intvectab->setDoc(doc_);

            iEvent.put(std::move(intsizetab), "intcounts");
            iEvent.put(std::move(floatsizetab), "floatcounts");
            iEvent.put(std::move(floatvectab), "floatvecs");
            iEvent.put(std::move(intvectab), "intvecs");
        }

    protected:
        const std::string name_; 
        const std::string doc_;
        const edm::EDGetTokenT<edm::View<T>> src_;
        const StringCutObjectSelector<T> cut_;
        int countPrecision_;
        std::vector<edm::EDGetTokenT<edm::ValueMap<std::vector<int>>>> intVecMaps_;
        std::vector<edm::EDGetTokenT<edm::ValueMap<std::vector<float>>>> floatVecMaps_;
        std::vector<std::string> intNames_;
        std::vector<std::string> floatNames_;
        std::vector<std::string> intDocs_;
        std::vector<std::string> floatDocs_;
        std::vector<int> intPrecisions_;
        std::vector<int> floatPrecisions_;
};
