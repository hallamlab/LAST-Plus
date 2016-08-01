// Copyright 2015 Martin C. Frith

#include "LastEvaluer.hh"

#include "GeneticCode.hh"

#include <algorithm>
#include <cstring>

#define COUNTOF(a) (sizeof (a) / sizeof *(a))

namespace cbrc {

  struct FrameshiftEvalueParameters {
    const char *matrixName;
    int gapOpen;
    int gapEpen;
    int frameshiftCost;
    Sls_P::AlignmentEvaluerParameters parameters;
  };

  const FrameshiftEvalueParameters frameshiftEvalueParameters[] = {
    {"BL62", 11, 1, 15, {0.31457181182385774, 0.077024909125411836,
                          3.5355057419386005, -39.014329998056937,
                          1.1739847579695837, -12.896364921780187,
                          160.29749789587885, -3200.2716722761552,
                          17.539096792506459, -349.06406999516196,
                          51.79266617536797, -1027.401229133852}},

    // this one is included only to make lest-test.sh faster:
    {"BL62", 11, 2, 12, {0.32704828292493776, 0.10543687626821494,
                          2.7166240938670798, -19.361170444340445,
                          0.89899140520122844, -6.29652445533966,
                          74.744763449386667, -1147.0060455603425,
                          8.3059555754127565, -127.46868078491309,
                          24.78179014970371, -379.14020451790986}},

    {"BL62", 11, 2, 15, {0.33388770870821022, 0.11532516007803961,
                          2.4678058049483518, -14.50532580281522,
                          0.82160372991753583, -4.8091552692419572,
                          55.103072639059761, -731.90592162187147,
                          6.1010321043683131, -80.763060603166991,
                          18.237750047538203, -240.59017890476582}},

    {"BL80", 11, 1, 15, {0.35400649542314511, 0.13270256108942211,
                          1.8960749679829285, -14.061923223904673,
                          0.62940827123451903, -4.6245065070665863,
                          35.18772909081801, -583.19242886423649,
                          3.8260214679558033, -62.789729751450864,
                          11.072568656496113, -178.63729131744145}},

    {"BL80", 11, 2, 15, {0.3652492341706855, 0.16850422398182774,
                          1.5316138005575799, -5.7577598061709985,
                          0.5101720323233776, -1.9097398376324572,
                          17.427364219333899, -170.0223112776693,
                          1.9259860816444827, -18.621287186644096,
                          5.7329546520583801, -54.693768145180513}},
  };

  static bool isEqual(const char *x, const char *y) {
    return std::strcmp(x, y) == 0;
  }

  static bool isHit(const FrameshiftEvalueParameters &p,
      const char *n, int a, int b, int f) {
    return isEqual(p.matrixName, n) &&
      p.gapOpen == a && p.gapEpen == b && p.frameshiftCost == f;
  }

  // Hardcoded to make it the first one. May refactor this later
  void LastEvaluer::init_LASTP(){
    const FrameshiftEvalueParameters &p = frameshiftEvalueParameters[0];
    return frameshiftEvaluer.initParameters(p.parameters);
  }
}
