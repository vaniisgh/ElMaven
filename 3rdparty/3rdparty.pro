include($$mac_compiler)

TEMPLATE = subdirs
CONFIG += ordered qt thread

SUBDIRS += obiwarp       \
           libneural     \
           libcdfread    \
           libcsvparser  \
           pugixml/src   \
           libpillow     \
           libpls        \
           Logger        \
           ErrorHandling \
           NimbleDSP     \
           libmgf

macx {
  SUBDIRS -= libcdfread
}
