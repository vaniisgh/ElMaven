#ifndef PEAKGROUP_H
#define PEAKGROUP_H

#include "datastructures/mzSlice.h"
#include "Fragment.h"
#include "Peak.h"
#include "standardincludes.h"

class mzSample;
class Isotope;
class MassCalculator;
class Compound;
class Peak;
class Scan;
class EIC;
class MassCutoff;

using namespace std;

class PeakGroup{

    public:
        enum GroupType {None=0, C13=1, Adduct=2, Covariant=4, Isotope=5 };     //group types
        enum QType	   {AreaTop=0,
                        Area=1,
                        Height=2,
                        AreaNotCorrected=3,
                        RetentionTime=4,
                        Quality=5,
                        SNRatio=6,
                        AreaTopNotCorrected=7};

        enum class ClassifiedLabel {
            None,                // the group has not been classified using ML yet
            Noise,               // the group is too noisy to be a metbolite
            Signal,              // the group shows a clear signal but may not be an interesting metabolite
            Correlation,         // the group is a signal and is correlated to one or more signals
            Pattern,             // the group is a signal and shows interesting inter-cohort intensity pattern
            CorrelationAndPattern // the group is a signal which shows correlation as well as intensity pattern
        };

        PeakGroup();
        PeakGroup(const PeakGroup& o);
        PeakGroup& operator=(const PeakGroup& o);

        bool operator==(const PeakGroup* o);
        /**
         * [copyObj ]
         * @method copyObj
         * @param  o       []
         */
        void copyObj(const PeakGroup& o);

        /**
         * [copy ]
         * @method copy
         * @param  o    []
         */
        void copy(const PeakGroup* o);

        ~PeakGroup();

        PeakGroup* parent;

        // have to do this since `GroupType` enum also has an Adduct.
        // In future use "enum class" instead. Also from MAVEN (upstream).
        class Adduct* adduct;

        vector<Peak> peaks;
        vector<PeakGroup> children;
        vector<PeakGroup> childrenBarPlot;
        vector<PeakGroup> childrenIsoWidget;
        vector<mzSample*> samples;  //this varibale will hold only those sample which has been
                                    //used for peak detection
        string srmId;
        string tagString;

        // Stores the name of Peak Table this group belongs to.
        string searchTableName;

        /**
         * @brief returns group name
         * @method getName
         * @detail returns compound name, tagString, srmID, meanMz@meanRt or groupId in this order of preference
         * @return string
         */
        string getName();

        bool isFocused;

        QType quantitationType;

        int groupId;
        int metaGroupId;
        int clusterId;

        bool deletedFlag;

        float maxIntensity;
        float maxAreaTopIntensity;
        float maxAreaIntensity;
        float maxHeightIntensity;
        float maxAreaNotCorrectedIntensity;
        float maxAreaTopNotCorrectedIntensity;
        float currentIntensity;
        float meanRt;
        float meanMz;
        float expectedMz;
        int totalSampleCount;

        int  ms2EventCount;

        FragmentationMatchScore fragMatchScore;
        Fragment fragmentationPattern;

        //isotopic information
        float expectedAbundance;
        int   isotopeC13count;

        double minIntensity;

        //int quantileIntensityPeaks;
        //int quantileQualityPeaks;

        float minRt;
        float maxRt;
        float minMz;
        float maxMz;

        float blankMax;
        float blankMean;
        unsigned int blankSampleCount;

        int sampleCount;
        float sampleMean;
        float sampleMax;

        unsigned int maxNoNoiseObs;
        unsigned int  maxPeakOverlap;
        float maxQuality;
        float avgPeakQuality;

        double minQuality;
        float maxPeakFracionalArea;
        float maxSignalBaseRatio;
        float maxSignalBaselineRatio;
        int goodPeakCount;
        float expectedRtDiff;
        float groupRank;

        //for sample contrasts  ratio and pvalue
        float changeFoldRatio;
        float changePValue;

        bool isMS1();

        /**
         * [hasSrmId ]
         * @method hasSrmId
         * @return []
         */
        bool  	hasSrmId() const   { return srmId.empty(); }

        /**
         * [setSrmId ]
         * @method setSrmId
         * @param  id       []
         * @return []
         */
        void  	setSrmId(string id)	  { srmId=id; }

        /**
         * [getSrmId ]
         * @method getSrmId
         * @return []
         */
        inline  string getSrmId() const { return srmId; }
    
        /**
         * [hasCompoundLink ]
         * @method hasCompoundLink
         * @return []
         */
        bool hasCompoundLink() const;

        /**
         * [isEmpty ]
         * @method isEmpty
         * @return []
         */
        inline bool isEmpty() const   { if(peaks.size() == 0) return true; return false; }

        /**
         * [peakCount ]
         * @method peakCount
         * @return []
         */
        inline unsigned int peakCount()  const { return peaks.size(); 	  }

        /**
         * [childCount ]
         * @method childCount
         * @return []
         */
        inline unsigned int childCount() const { return children.size(); }

        inline unsigned int childCountBarPlot() const { return childrenBarPlot.size(); }

        inline unsigned int childCountIsoWidget() const { return childrenIsoWidget.size(); }

        Compound* getCompound();

        void setCompound(Compound* compound);

        void setSlice(const mzSlice& slice);

        const mzSlice& getSlice() const;

        /**
         * @brief Check whether a slice has previosuly been set for this group.
         * @return true if a slice has been set, false otherwise.
         */
        bool hasSlice() const;

        /**
         * @brief Check whether both bounds of the group's slice are close to
         * zero in either m/z or rt dimensions.
         * @return `true` if either slice's `mzmin` and `mzmax` are both close
         * to zero or slice's `rtmin` and `rtmax` are both close to zero, false
         * otherwise.
         */
        bool sliceIsZero() const;

        /**
         * [getParent ]
         * @method getParent
         * @return []
         */
        inline PeakGroup* getParent() { return parent; }


        inline vector<Peak>& getPeaks() { return peaks; }


        inline vector<PeakGroup>& getChildren()  { return children; }

        vector<Scan*> getRepresentativeFullScans(); //TODO: Sahil - Kiran, Added while merging mainwindow

        /**
         * @brief find all MS2 scans for this group
         * @return vector of all MS2 scans for this group
         */
        vector<Scan*> getFragmentationEvents();

        /**
         * @brief build a consensus fragment spectra for this group
         * @param productPpmTolr ppm tolerance for fragment m/z
         */
        void computeFragPattern(float productPpmTolr);

        Scan* getAverageFragmentationScan(float productPpmTolr);

        void matchFragmentation(float ppmTolerance, string scoringAlgo);
        
        double getExpectedMz(int charge);

        /**
         * [setParent ]
         * @method setParent
         * @param  p         []
         */
        inline void setParent(PeakGroup* p) {parent=p;}

        /**
         * @brief Assign a label to this peak group.
         * @details If this method is called with one of the valid user labels,
         * then user decision overrides any PeakML classified predictions.
         * @param label A character that is either 'g' (good), 'b' (bad) or
         * '\0' (unset).
         */
        void setUserLabel(const char label);

        /**
         * @brief Get the current good/bad label assigned to this peak group.
         * @return A character, denoting a good peak if 'g', a bad peak if 'b'
         * or '\0' if unassigned.
         */
        char userLabel() const { return _userLabel; }

        /**
         * [ppmDist ]
         * @method ppmDist
         * @param  cmass   []
         * @return []
         */
        float massCutoffDist(float cmass,MassCutoff *massCutoff);

        /**
         * [addPeak ]
         * @method addPeak
         * @param  peak    []
         */
        void addPeak(const Peak& peak); 

        /**
         * [addChild ]
         * @method addChild
         * @param  child    []
         */
        inline void addChild(const PeakGroup& child) { children.push_back(child); children.back().parent = this;   }

        inline void addChildBarPlot(const PeakGroup& child) { childrenBarPlot.push_back(child); childrenBarPlot.back().parent = this;   }

        inline void addChildIsoWidget(const PeakGroup& child) { childrenIsoWidget.push_back(child); childrenIsoWidget.back().parent = this;   }

        inline void setGroupIdForChildren() { for (auto& child : children) child.groupId = groupId; }

        /**
         * [getPeak ]
         * @method getPeak
         * @param  sample  []
         * @return []
         */

        Peak* getPeak(mzSample* sample);

        GroupType _type;

        /**
         * [type ]
         * @method type
         * @return []
         */
        inline GroupType type() const { return _type; }
        /**
         * [setType ]
         * @method setType
         * @param  t       []
         */
        inline void setType(GroupType t)  { _type = t; }

        void setQuantitationType(QType type) {quantitationType = type;}

        /**
         * [isIsotope ]
         * @method isIsotope
         * @return []
         */
        inline bool isIsotope() const { return _type == Isotope; }

        /**
         * [isAdduct ]
         * @method isAdduct
         * @return []
         */
        inline bool isAdduct() const {  return _type == Adduct; }

        /**
         * [summary ]
         * @method summary
         */
        void summary();

        /**
         * [groupStatistics ]
         * @method groupStatistics
         */
        void groupStatistics();

        /**
         * [updateQuality ]
         * @method updateQuality
         */
        void updateQuality();

        /**
         * [medianRt ]
         * @method medianRt
         * @return []
         */
        float medianRt();

        /**
         * [meanRtW ]
         * @method meanRtW
         * @return []
         */
        float meanRtW();

        /**
         * [reduce ]
         * @method reduce
         */
        void reduce();

        /**
         * [fillInPeaks ]
         * @method fillInPeaks
         * @param  eics        []
         */
        void fillInPeaks(const vector<EIC*>& eics);

        /**
         * [computeAvgBlankArea ]
         * @method computeAvgBlankArea
         * @param  eics                []
         */
        void computeAvgBlankArea(const vector<EIC*>& eics);

        /**
         * [groupOverlapMatrix ]
         * @method groupOverlapMatrix
         */
        void groupOverlapMatrix();

        /**
         * [getSamplePeak ]
         * @method getSamplePeak
         * @param  sample        []
         * @return []
         */
        Peak* getSamplePeak(mzSample* sample);

        /**
         * [deletePeaks ]
         * @method deletePeaks
         */
        void deletePeaks();

        /**
         * [deletePeak ]
         * @method deletePeak
         * @param  index      []
         * @return []
         */
        bool deletePeak(unsigned int index);

        /**
         * [clear ]
         * @method clear
         */
        void clear();

        /**
         * [deleteChildren ]
         * @method deleteChildren
         */
        void deleteChildren();

        /**
         * [deleteChild ]
         * @method deleteChild
         * @param  index       []
         * @return []
         */
        bool deleteChild(unsigned int index);

        /**
         * [deleteChild ]
         * @method deleteChild
         * @param  child       []
         * @return []
         */
        bool deleteChild(PeakGroup* child);

        /**
         * [copyChildren ]
         * @method copyChildren
         * @param  other        []
         */
        void copyChildren(const PeakGroup& other);

        vector<float> getOrderedIntensityVector(vector<mzSample*>& samples, QType type);

        /**
         * [reorderSamples ]
         * @method reorderSamples
         */
        void reorderSamples();

        /**
         * [compRt ]
         * @method compRt
         * @param  a      []
         * @param  b      []
         * @return []
         */
        static bool compRt(const PeakGroup& a, const PeakGroup& b ) { return(a.meanRt < b.meanRt); }

        /**
         * [compMz ]
         * @method compMz
         * @param  a      []
         * @param  b      []
         * @return []
         */
        static bool compMz(const PeakGroup& a, const PeakGroup& b ) { return(a.meanMz > b.meanMz); }

        /**
         * [compIntensity ]
         * @method compIntensity
         * @param  a             []
         * @param  b             []
         * @return []
         */
        static bool compIntensity(const PeakGroup& a, const PeakGroup& b ) { return(a.maxIntensity > b.maxIntensity); }

        /**
         * [compArea ]
         * @method compArea
         * @param  a        []
         * @param  b        []
         * @return []
         */
        static bool compArea(const PeakGroup& a, const PeakGroup& b ) { return(a.maxPeakFracionalArea > b.maxPeakFracionalArea); }

        /**
         * [compQuality ]
         * @method compQuality
         * @param  a           []
         * @param  b           []
         * @return []
         */
        static bool compQuality(const PeakGroup& a, const PeakGroup& b ) { return(a.maxQuality > b.maxQuality); }
        //static bool compInfoScore(const PeakGroup& a, const PeakGroup& b ) { return(a.informationScore > b.informationScore); }

        /**
         * [compRank ]
         * @method compRank
         * @param  a        []
         * @param  b        []
         * @return []
         */
        static bool compRank(const PeakGroup& a, const PeakGroup& b ) { return(a.groupRank < b.groupRank); }

        /**
         * [compRatio ]
         * @method compRatio
         * @param  a         []
         * @param  b         []
         * @return []
         */
        static bool compRatio(const PeakGroup& a, const PeakGroup& b ) { return(a.changeFoldRatio < b.changeFoldRatio); }

        /**
         * [compPvalue ]
         * @method compPvalue
         * @param  a          []
         * @param  b          []
         * @return []
         */
        static bool compPvalue(const PeakGroup* a, const PeakGroup* b ) { return(a->changePValue< b->changePValue); }

        /**
         * [compC13 ]
         * @method compC13
         * @param  a       []
         * @param  b       []
         * @return []
         */
        static bool compC13(const PeakGroup* a, const PeakGroup* b) { return(a->isotopeC13count < b->isotopeC13count); }

        /**
         * [compMetaGroup ]
         * @method compMetaGroup
         * @param  a             []
         * @param  b             []
         * @return []
         */
        static bool compMetaGroup(const PeakGroup& a, const PeakGroup& b) { return(a.metaGroupId < b.metaGroupId); }
        bool operator< (const PeakGroup* b) const { return this->maxIntensity < b->maxIntensity; }

        void calGroupRank(bool deltaRtCheckFlag,
                            int qualityWeight,
                            int intensityWeight,
                            int deltaRTWeight);
        /**
         * @brief take list of sample and filter those which are marked as selected
         * @details this method used for assigning samples to this group based on whether that samples
         * are marked as selected.
        */
        void setSelectedSamples(vector<mzSample*> vsamples);

        /**
         * @brief Set the classification label for this PeakGroup.
         * @details Anything other than `PeakGroup::ClassifiedLabel::Noise` or
         * `PeakGroup::ClassifiedLabel::None` will be marked as a good ('g')
         * group.
         * @param label The PeakML class to assign.
         * @param probability The probability for this classification.
         */
        void setPredictedLabel(const ClassifiedLabel label,
                               const float probability);

        /**
         * @brief Set the inference values for assigned classification.
         * @param inference A multi-map of computed values mapping to each peak
         * attribute. A `multimap` is used so that "views" can later on use it
         * to display top N features responsible for classification, which would
         * be a tedious process if we used a map of attributes:values instead.
         */
        void setPredictionInference(const multimap<float, string>& inference);

        /**
         * @brief Get the predicted label for this peak group.
         * @return A `PeakGroup::ClassifiedLabel` denoting the predicted class.
         */
        ClassifiedLabel predictedLabel() const;

        /**
         * @brief The proability with which this group has been assigned its
         * current PeakML class.
         * @return A floating point value representing confidence of
         * prediction.
         */
        float predictionProbability() const;

        /**
         * @brief Obtain inferences that can help illustrate why the peak-group
         * was assigned its current probability.
         * @return A `std::multimap` storing values mapping to every known
         * attribute.
         */
        multimap<float, string> predictionInference() const;

        /**
         * @brief Converts an integer prediction-class identifier to its
         * corresponding `PeakGroup::ClassifiedLabel`.
         * @details This function should be regarded as complementary to
         * `PeakGroup::integralValueForClass`, assisting in conversion from a
         * primitive representation (integer).
         * @param An integer which has to be converted to a `ClassifiedLabel`.
         * @return A `PeakGroup::ClassifiedLabel`.
         */
        static ClassifiedLabel classificationLabelForValue(int value)
        {
            if (value == 0)
                return PeakGroup::ClassifiedLabel::Noise;
            if (value == 1)
                return PeakGroup::ClassifiedLabel::Signal;
            if (value == 2)
                return PeakGroup::ClassifiedLabel::Correlation;
            if (value == 3)
                return PeakGroup::ClassifiedLabel::Pattern;
            if (value == 4)
                return PeakGroup::ClassifiedLabel::CorrelationAndPattern;
            return PeakGroup::ClassifiedLabel::None;
        };

        /**
         * @brief Converts a `PeakGroup::ClassifiedLabel` identifier to an
         * integer value, guaranteed to be unique for all possible classes.
         * @details This function should be regarded as complementary to
         * `PeakGroup::classificationLabelForValue`, assisting in conversion to
         * a primitive representation (integer).
         * @param label A `PeakGroup::ClassifiedLabel` that needs to be
         * converted to an integer.
         * @return An integer, unique for each enum value.
         */
        static int integralValueForLabel(ClassifiedLabel label)
        {
            if (label == PeakGroup::ClassifiedLabel::Noise)
                return 0;
            if (label == PeakGroup::ClassifiedLabel::Signal)
                return 1;
            if (label == PeakGroup::ClassifiedLabel::Correlation)
                return 2;
            if (label == PeakGroup::ClassifiedLabel::Pattern)
                return 3;
            if (label == PeakGroup::ClassifiedLabel::CorrelationAndPattern)
                return 4;
            return -1;
        }

    private:
        mzSlice _slice;
        bool _sliceSet;

        // user classification label
        char _userLabel;

        // properties for PeakML classification
        ClassifiedLabel _predictedLabel;
        float _predictionProbability;
        multimap<float, string> _predictionInference;
};
#endif
