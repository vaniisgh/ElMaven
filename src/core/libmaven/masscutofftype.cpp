#include "masscutofftype.h"

double MassCutoff::massCutoffValue(double mz){
    
    if(_massCutoffType=="ppm"){
        //cerr<<"mass cutoff type:  "<<_massCutoffType<<"  value: "<<_massCutoff<<endl;
        return _massCutoff*mz/1e6;
    }
    else if(_massCutoffType=="mDa"){
        //cerr<<"mass cutoff type:  "<<_massCutoffType<<"  value: "<<_massCutoff<<endl;
        return _massCutoff;
    }
    else {
        //cerr<<"mass cutoff type:  "<<"unknown"<<"  value: "<<0<<endl;
        return 0;
    }
}

void MassCutoff::setMassCutoffAndType(double massCutoff, string massCutoffType){
    _massCutoffType=massCutoffType;
    _massCutoff=massCutoff;
}