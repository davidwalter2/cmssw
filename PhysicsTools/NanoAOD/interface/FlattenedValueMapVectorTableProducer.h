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

            // One (Counts + Vals) table pair per variable, each its own named
            // collection. This is the 10_6 behaviour: distinct-length vectors
            // (e.g. globalIdxs, jacRef, momCov) each get their own table, so
            // there is no "same length" constraint and no collision with the
            // main object table.
            for (size_t i = 0; i < intVecMaps_.size(); i++) {
                produces<nanoaod::FlatTable>("intcounts"+std::to_string(i));
                produces<nanoaod::FlatTable>("intvec"+std::to_string(i));
            }
            for (size_t i = 0; i < floatVecMaps_.size(); i++) {
                produces<nanoaod::FlatTable>("floatcounts"+std::to_string(i));
                produces<nanoaod::FlatTable>("floatvec"+std::to_string(i));
            }
        }

        ~FlattenedValueMapVectorTableProducer() override {}

        template<typename P>
        std::vector<P> readVals(const edm::ValueMap<std::vector<P>>& vmap, edm::PtrVector<T>& objs, std::vector<int>& sizes) {
            std::vector<P> allvals;
            // Will be at least this long
            allvals.reserve(objs.size());
            for (size_t i = 0; i < objs.size(); i++) {
                auto& vals = vmap[objs[i]];
                sizes[i] = vals.size();
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

            // One (Counts + Vals) table pair per variable, each with its own
            // named collection (name_ + "_" + var[ + "_Counts"]). No same-length
            // constraint across variables, and no collision with the main
            // object table (name_). Matches the 10_6 behaviour.
            std::vector<int> sizes(objs.size(), 0);
            for (size_t i = 0; i < intVecMaps_.size(); i++) {
                edm::Handle<edm::ValueMap<std::vector<int>>> vmap;
                iEvent.getByToken(intVecMaps_[i], vmap);
                const auto& results = readVals(*vmap, objs, sizes);

                auto countstab = std::make_unique<nanoaod::FlatTable>(objs.size(), this->name_ + "_" + intNames_[i] + "_Counts", false, false);
                countstab->template addColumn<int>("", sizes, "Number of entries per object", countPrecision_);
                countstab->setDoc(doc_);
                auto vectab = std::make_unique<nanoaod::FlatTable>(results.size(), this->name_ + "_" + intNames_[i], false, false);
                vectab->template addColumn<int>("Vals", results, intDocs_[i], intPrecisions_[i]);
                vectab->setDoc(doc_);

                iEvent.put(std::move(countstab), "intcounts"+std::to_string(i));
                iEvent.put(std::move(vectab), "intvec"+std::to_string(i));
            }
            for (size_t i = 0; i < floatVecMaps_.size(); i++) {
                edm::Handle<edm::ValueMap<std::vector<float>>> vmap;
                iEvent.getByToken(floatVecMaps_[i], vmap);
                const auto& results = readVals(*vmap, objs, sizes);

                auto countstab = std::make_unique<nanoaod::FlatTable>(objs.size(), this->name_ + "_" + floatNames_[i] + "_Counts", false, false);
                countstab->template addColumn<int>("", sizes, "Number of entries per object", countPrecision_);
                countstab->setDoc(doc_);
                auto vectab = std::make_unique<nanoaod::FlatTable>(results.size(), this->name_ + "_" + floatNames_[i], false, false);
                vectab->template addColumn<float>("Vals", results, floatDocs_[i], floatPrecisions_[i]);
                vectab->setDoc(doc_);

                iEvent.put(std::move(countstab), "floatcounts"+std::to_string(i));
                iEvent.put(std::move(vectab), "floatvec"+std::to_string(i));
            }
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
