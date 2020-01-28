#include "numeric_treewidgetitem.h"
#include "PeakGroup.h"

bool NumericTreeWidgetItem::operator<( const QTreeWidgetItem & other ) const{
    int sortCol = treeWidget()->sortColumn();

    // takes care of sorting based on PeakML class labels
    QString thisLabelString = this->data(sortCol, Qt::UserRole).value<QString>();
    QString otherLabelString = other.data(sortCol, Qt::UserRole).value<QString>();
    auto thisLabel = PeakGroup::labelForString(thisLabelString.toStdString());
    auto otherLabel = PeakGroup::labelForString(otherLabelString.toStdString());
    if (thisLabel != PeakGroup::ClassifiedLabel::None
        || otherLabel != PeakGroup::ClassifiedLabel::None) {
        return thisLabel < otherLabel;
    }

    QString thisText = text(sortCol);
    QString otherText = other.text(sortCol);

    QCollator collator;
    collator.setNumericMode(true);

    QRegExp exp("^[-+]?([0-9]*\\.?[0-9]+)[eE][-+]?([0-9]+)$");
    if (thisText.contains(exp)) {
        double mantissa = exp.cap(1).toDouble();
        double exponent = exp.cap(2).toDouble();
        double trueValue = mantissa * (pow(10, exponent));
        thisText = QString::fromStdString(to_string(trueValue).c_str());
    }
    if (otherText.contains(exp)) {
        double mantissa = exp.cap(1).toDouble();
        double exponent = exp.cap(2).toDouble();
        double trueValue = mantissa * (pow(10, exponent));
        otherText = QString::fromStdString(to_string(trueValue).c_str());
    }

    return collator.compare(thisText , otherText) < 0;
}
