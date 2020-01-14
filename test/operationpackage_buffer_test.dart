import "package:test/test.dart";
import 'package:dart_jts/dart_jts.dart';
import "dart:math" as math;
import 'testing_utilities.dart';

double TOLERANCE = 1E-5;

PrecisionModel precisionModel = PrecisionModel.fixedPrecision(1);

GeometryFactory geometryFactory = GeometryFactory.withPrecisionModelSrid(precisionModel, 0);

WKTReader reader = WKTReader.withFactory(geometryFactory);

void main() {
  group("Buffer tests - ", () {
    test("test1", () {
      BufferValidator(0, "POINT (100 100)").setEmptyBufferExpected(true).test();
    });
    test("test2", () {
      BufferValidator(0, "LINESTRING (10 10, 100 100)").setEmptyBufferExpected(true).test();
    });
    test("test1a", () {
      BufferValidator(-1, "POINT (100 100)").setEmptyBufferExpected(true).test();
    });
    test("test2a", () {
      BufferValidator(-1, "LINESTRING (10 10, 100 100)").setEmptyBufferExpected(true).test();
    });
    test("test3", () {
      BufferValidator(10, "LINESTRING (100 100, 200 100, 200 200, 100 200, 100 100)").test();
    });
    test("test4", () {
      BufferValidator(50, "LINESTRING (40 40, 160 40, 100 180, 40 80)").test();
    });
    test("test5", () {
      BufferValidator(0, "POLYGON ((80 300, 280 300, 280 300, 280 300, 280 80, 80 80, 80 300))").test();
    });
    test("test6", () {
      BufferValidator(10, "POLYGON ((60 300, 60 160, 240 160, 240 300, 60 300))").test();
    });
    test("test7", () {
      BufferValidator(10, "POLYGON ((80 300, 280 300, 280 80, 80 80, 80 300), (260 280, 180 200, 100 280, 100 100, 260 100, 260 280))").test();
    });
    test("test8", () {
      BufferValidator(200, "POLYGON ((80 300, 280 300, 280 80, 80 80, 80 300), (260 280, 180 200, 100 280, 100 100, 260 100, 260 280))").test();
    });
    test("test9", () {
      BufferValidator(-10, "POLYGON ((80 300, 280 300, 280 80, 80 80, 80 300))").test();
    });
    test("test10", () {
      BufferValidator(10, "POLYGON ((100 300, 300 300, 300 100, 100 100, 100 300), (220 220, 180 220, 180 180, 220 180, 220 220))").test();
    });
    test("test11", () {
      BufferValidator(5,
              "POLYGON ((260 400, 220 300, 80 300, 180 220, 40 200, 180 160, 60 20, 200 80, 280 20, 260 140, 440 20, 340 180, 520 160, 280 220, 460 340, 300 300, 260 400), (260 320, 240 260, 220 220, 160 180, 220 160, 200 100, 260 160, 300 140, 320 180, 260 200, 260 320))")
          .test();
    });
    test("test12", () {
      BufferValidator(-17, "POLYGON ((260 320, 240 260, 220 220, 160 180, 220 160, 260 160, 260 200, 260 320))").test();
    });
    test("test13", () {
      BufferValidator(-17, "POLYGON ((260 320, 240 260, 220 220, 260 160, 260 320))").test();
    });
    test("test14", () {
      BufferValidator(-14, "POLYGON ((260 320, 240 260, 220 220, 260 160, 260 320))").test();
    });
    test("test15", () {
      BufferValidator(26, "LINESTRING (260 160, 260 200, 260 320, 240 260, 220 220)").test();
    });
    test("test16", () {
      BufferValidator(-7,
              "POLYGON ((260 400, 220 300, 80 300, 180 220, 40 200, 180 160, 60 20, 200 80, 280 20, 260 140, 440 20, 340 180, 520 160, 280 220, 460 340, 300 300, 260 400), (260 320, 240 260, 220 220, 160 180, 220 160, 200 100, 260 160, 300 140, 320 180, 260 200, 260 320))")
          .test();
    });
    test("test17", () {
      BufferValidator(-183,
              "POLYGON ((32 136, 27 163, 30 236, 34 252, 49 291, 72 326, 83 339, 116 369, 155 391, 176 400, 219 414, 264 417, 279 416, 339 401, 353 395, 380 381, 394 372, 441 328, 458 303, 463 294, 480 251, 486 205, 486 183, 473 115, 469 105, 460 85, 454 74, 423 33, 382 2, 373 -3, 336 -19, 319 -24, 275 -31, 252 -32, 203 -27, 190 -24, 149 -10, 139 -5, 84 37, 76 46, 52 81, 36 121, 32 136))")
          .test();
    });
    test("test18", () {
      BufferValidator(20, "POLYGON((-4 225, -17 221, -16 223, -15 224, -13 227, -4 225))").test();
    });
    test("test19", () {
      BufferValidator(21, "POLYGON ((184 369, 181 368, 180 368, 179 367, 176 366, 185 357, 184 369 ))").test();
    });
    test("test20", () {
      BufferValidator(1000, "POLYGON ((13841 1031, 13851 903, 13853 885, 13853 875, 13856 862, 13859 831, 13670 900, 13841 1031))").test();
    });
    test("test21", () {
      BufferValidator(18, "POLYGON ((164 84, 185 91, 190 75, 187 76, 182 77, 179 79, 178 79, 174 81, 173 81, 172 82, 169 83,  164 84 ))").test();
    });
    test("test22", () {
      BufferValidator(15, "POLYGON ((224 271, 225 261, 214 258, 210 266, 212 267, 214 267, 217 268, 218 268, 219 268, 221 269, 222 270,  224 271 ))").test();
    });
    test("test23", () {
      BufferValidator(25,
              "POLYGON ((484 76, 474 79, 492 122, 502 119, 501 117, 500 112, 499 111, 498 107, 497 104, 496 103, 494 98, 493 96, 491 92, 490 90, 489 86, 487 81, 486 79, 485 77, 484 76 ))")
          .test();
    });
    test("test24", () {
      BufferValidator(160, "POLYGON ((20 60, 20 20, 240 20, 40 21, 240 22, 40 22, 240 23, 240 60, 20 60))").test();
    });
    test("test25", () {
      BufferValidator(-3,
              "POLYGON ((233 195, 232 195, 231 194, 222 188, 226 187, 227 187, 229 187, 230 186, 232 186, 234 185, 236 184, 237 183, 238 182, 237 184, 236 185, 236 186, 235 187, 235 188, 234 189, 234 191, 234 192, 233 193, 233 195 ))")
          .test();
    });
    test("test26", () {
      BufferValidator(6,
              "LINESTRING (233 195, 232 195, 231 194, 222 188, 226 187, 227 187, 229 187, 230 186, 232 186, 234 185, 236 184, 237 183, 238 182, 237 184, 236 185, 236 186, 235 187, 235 188, 234 189, 234 191, 234 192, 233 193, 233 195 )")
          .test();
    });
    test("test27", () {
      BufferValidator(-30,
              "POLYGON ((2330 1950, 2320 1950, 2310 1940, 2220 1880, 2260 1870, 2270 1870, 2290 1870, 2300 1860, 2320 1860, 2340 1850, 2360 1840, 2370 1830, 2380 1820, 2370 1840, 2360 1850, 2360 1860, 2350 1870, 2350 1880, 2340 1890, 2340 1910, 2340 1920, 2330 1930, 2330 1950 ))")
          .test();
    });
    test("test28", () {
      BufferValidator(30,
              "LINESTRING (2330 1950, 2320 1950, 2310 1940, 2220 1880, 2260 1870, 2270 1870, 2290 1870, 2300 1860, 2320 1860, 2340 1850, 2360 1840, 2370 1830, 2380 1820, 2370 1840, 2360 1850, 2360 1860, 2350 1870, 2350 1880, 2340 1890, 2340 1910, 2340 1920, 2330 1930, 2330 1950 )")
          .test();
    });
    test("test29", () {
      BufferValidator(26,
              "POLYGON ((440 -93, 440 -67, 475 -67, 471 -71, 469 -72, 468 -73, 467 -74, 463 -78, 459 -81, 458 -82, 454 -84, 453 -85, 452 -86, 450 -86, 449 -87, 448 -88, 444 -90, 443 -91, 441 -92, 440 -93 ))")
          .test();
    });
    test("test30", () {
      BufferValidator(260,
              "POLYGON ((4400 -930, 4400 -670, 4750 -670, 4710 -710, 4690 -720, 4680 -730, 4670 -740, 4630 -780, 4590 -810, 4580 -820, 4540 -840, 4530 -850, 4520 -860, 4500 -860, 4490 -870, 4480 -880, 4440 -900, 4430 -910, 4410 -920, 4400 -930 ))")
          .test();
    });
    test("test31", () {
      BufferValidator(0.1,
              "POLYGON ((635074.6769928858 6184832.427381967, 635075.6723193424 6184799.950949265, 634717.5983159657 6184655.107092909, 634701.0176852546 6184648.498845058, 634697.7188197445 6184647.20632975, 634694.416887708 6184645.922033237, 634691.1138635761 6184644.642692243, 634687.8077729489 6184643.371570057, 634684.498667351 6184642.107006015, 634681.1875340013 6184640.847368483, 634677.8742698929 6184639.595978798, 634674.5570551592 6184638.351118257, 634671.2386969016 6184637.112873929, 634667.9173237421 6184635.881187774, 634664.5938713895 6184634.656088823, 634661.2674041622 6184633.437548058, 634657.9388577675 6184632.2255945075, 634654.6082322216 6184631.02022817, 634651.2745403448 6184629.823080709, 634647.9388208436 6184628.630859804, 634644.6000865338 6184627.4451971175, 634641.2592216335 6184626.267782336, 634637.9163291481 6184625.095294129, 634634.5713061031 6184623.931053837, 634631.2232683088 6184622.773371783, 634636.1918816608 6184608.365992378, 634633.2495506873 6184607.353869728, 634630.3051410569 6184606.348333739, 634627.3587557608 6184605.346063082, 634624.4102918282 6184604.3503790945, 634621.4607364619 6184603.359650123, 634618.5091539674 6184602.37384716, 634615.5564800596 6184601.392999219, 634612.6017790422 6184600.417077295, 634609.6450509242 6184599.446081388, 634606.6862442375 6184598.481672177, 634603.7263976521 6184597.52055733, 634600.7654082242 6184596.566058185, 634597.80145603 6184595.61645607, 634594.8364124894 6184594.671808995, 634591.8702261405 6184593.733777636, 634588.9020642313 6184592.799011653, 634585.9318238292 6184591.870832384, 634582.960543591 6184590.945947501, 634579.9871848791 6184590.027649342, 634577.0127348808 6184589.11430624, 634574.0362578988 6184588.205889201, 634571.0586381858 6184587.304087893, 634568.0790429743 6184586.405551985, 634565.0983050519 6184585.513631817, 634562.115540186 6184584.626637723, 634559.1316840936 6184583.744598703, 634556.1458010782 6184582.867485766, 634553.1587753976 6184581.996988578, 634550.1697742692 6184581.12975681, 634547.1796305005 6184580.269140797, 634544.1874598458 6184579.41345088, 634541.194198027 6184578.562716054, 634538.1998450506 6184577.716936321, 634535.2034137691 6184576.877743362, 634532.2059428038 6184576.041844833, 634529.2063935531 6184575.212533085, 634526.2057531879 6184574.388176442, 634523.2040217179 6184573.568774906, 634520.2002119992 6184572.755960159, 634517.1953626422 6184571.946439856, 634514.1893707667 6184571.143535337, 634510.267712847 6184585.871039091, 634281.9449709259 6184525.076957544, 633860.4859191478 6184412.861324424, 633664.3557212166 6184360.639468017, 633645.5884675509 6184355.641948889, 633486.222 6184313.208, 633485.7474265156 6184328.852301474, 633485.2749953512 6184344.496113185, 633650.4562371405 6184388.478170839, 633669.5206846121 6184393.553017912, 633852.6461183216 6184442.312440121, 634280.9949861752 6184556.364455, 634502.4254528129 6184615.324425217, 634505.716566367 6184616.204307566, 634509.0065372197 6184617.090806118, 634512.2953653594 6184617.983920872, 634515.5812308139 6184618.88193318, 634518.8659020835 6184619.788222348, 634522.1484948951 6184620.7010987215, 634525.4299963829 6184621.61893061, 634528.7093679372 6184622.545010364, 634531.9867124417 6184623.476016638, 634535.2619784358 6184624.413610098, 634538.5360501477 6184625.3594804, 634541.807159074 6184626.310248219, 634545.0771251463 6184627.267632202, 634548.3459483465 6184628.231632348, 634551.611808724 6184629.200529998, 634554.8764747474 6184630.177704472, 634558.138126462 6184631.161437107, 634561.3986867118 6184632.150125222, 634564.6571168633 6184633.147061154, 634567.9135198268 6184634.14892356, 634571.1678441261 6184635.157373106, 634574.421025444 6184636.172438785, 634577.6711923747 6184637.194062597, 634580.9192805979 6184638.222273537, 634584.1662257991 6184639.257100599, 634587.410208043 6184640.296825108, 634590.6529957657 6184641.3448264, 634593.8928205046 6184642.397725132, 634597.131502168 6184643.457239971, 634600.3671178726 6184644.524973575, 634603.6016934284 6184645.596001944, 634605.6958877691 6184646.2958191885, 634606.946276627 6184646.713661825, 634608.6177147877 6184647.275847967, 634610.2887808911 6184647.8379089665, 634613.6292576884 6184648.967082693, 634616.9666683461 6184650.104475331, 634620.3029357173 6184651.248484221, 634623.6361369188 6184652.4007120095, 634626.9663749042 6184653.557837355, 634630.2954695637 6184654.721578941, 634633.6214980055 6184655.893539405, 634636.9454988878 6184657.070426424, 634640.2664335228 6184658.255532312, 634643.5852890809 6184659.44722541, 634646.9020655481 6184660.6455057105, 634650.2167629072 6184661.850373212, 634653.5284454564 6184663.061798892, 634656.8370616816 6184664.2814434115, 634660.1436502574 6184665.506014449, 634663.4481081718 6184666.7388333315, 634666.7496027217 6184667.976549702, 634670.0489665787 6184669.222513908, 634673.3453155413 6184670.475036258, 634676.6396367943 6184671.732485097, 634679.9308916207 6184672.9981527375, 634683.2200672035 6184674.270407524, 634686.5062278372 6184675.54922043, 634689.789322 6184676.8362521175, 634693.0703883899 6184678.128210268, 634696.3493754807 6184679.426755545, 634699.6253475712 6184680.731858918, 634702.8982531315 6184682.04518105, 634706.1681951834 6184683.363400595, 635074.6769928858 6184832.427381967))")
          .test();
    });
    test("test32", () {
      BufferValidator(30,
              "MULTILINESTRING ((80 285, 85.5939933259177 234.65406006674084 ), (85.5939933259177 234.65406006674084, 98 123, 294 92, 344.3694502052736 126.0884157954882 ), (344.3694502052736 126.0884157954882, 393 159 ), (51 235, 85.5939933259177 234.65406006674084 ), (85.5939933259177 234.65406006674084, 251 233, 344.3694502052736 126.0884157954882 ), (344.3694502052736 126.0884157954882, 382 83 ))")
          .test();
    });
    test("test34", () {
      BufferValidator(1, "GEOMETRYCOLLECTION (POLYGON ((0 10, 10 0, 10 10, 0 10),  (4 8, 8 4, 8 8, 4 8)),   LINESTRING (6 6, 20 20))").test();
    });
    test("test35", () {
      BufferValidator(20,
              "GEOMETRYCOLLECTION (POINT (100 100), POLYGON ((400 260, 280 380, 240 220, 120 300, 120 100, 260 40, 200 160, 400 260)), LINESTRING (260 400, 220 280, 120 400, 20 280, 160 160, 60 40, 160 20, 360 140))")
          .test();
    });
    test("test36", () {
      BufferValidator(20, "GEOMETRYCOLLECTION (POINT (100 100), POLYGON ((400 260, 120 300, 120 100, 400 260)), LINESTRING (20 280, 160 160, 60 40))").test();
    });
    test("test37", () {
      BufferValidator(300, "POLYGON ((-140 700, 880 1120, 1280 -120, 300 -600, -480 -480, -140 700),   (0 360, 780 500, 240 -220, 0 360))").test();
    });
    test("test38", () {
      BufferValidator(300, "POLYGON ((-140 700, 880 1120, 1280 -120, 300 -600, -480 -480, -140 700),   (0 360, 240 -220, 780 500, 0 360))").test();
    });
    test("test39", () {
      BufferValidator(30, "MULTIPOLYGON (((0 400, 440 400, 440 0, 0 0, 0 400),(380 360, 160 120, 260 80, 380 360)), ((360 320, 200 120, 240 100, 360 320)))")
          .test();
    });
    test("test", () {});
    test("test", () {});
    test("test", () {});
    test("test", () {});
    test("test", () {});
  });
}

class Test implements Comparable {
  int priority;
  String name;
  Function testFunction;

  Test(String name, Function testFunction, {int priority}) {
    this.name = name;
    this.priority = priority;
    this.testFunction = testFunction;
  }

  String getName() {
    return name;
  }

  String toString() {
    return getName();
  }

  void test() {
    testFunction();
  }

  int compareTo(Object o) {
    return priority - (o as Test).priority;
  }
}

/**
 * @version 1.7
 */
class BufferValidator {
//   static void main(String[] args) throws Exception {
//  Geometry g =
//  new WKTReader().read(
//      "MULTILINESTRING (( 635074.5418406526 6184832.4888257105, 635074.5681951842 6184832.571842485, 635074.6472587794 6184832.575795664 ), ( 635074.6657069515 6184832.53889932, 635074.6933792098 6184832.451929366, 635074.5642420045 6184832.474330718 ))");
////System.out.println(g);
////System.out.println(g.buffer(0.01, 100));
////System.out.println("END");
//}

  Geometry original;
  double bufferDistance;
  Map nameToTestMap = {};
  Geometry buffer;
  static final int QUADRANT_SEGMENTS_1 = 100;
  static final int QUADRANT_SEGMENTS_2 = 50;
  String wkt;
  GeometryFactory geomFact = new GeometryFactory.defaultPrecision();
  WKTWriter wktWriter = new WKTWriter();
  WKTReader wktReader;

  BufferValidator(double bufferDistance, String wkt) : this.withBool(bufferDistance, wkt, true);

  BufferValidator.withBool(double bufferDistance, String wkt, bool _addContainsTest) {
// SRID = 888 is to test that SRID is preserved in computed buffers
    setFactory(new PrecisionModel(), 888);
    this.bufferDistance = bufferDistance;
    this.wkt = wkt;
    if (_addContainsTest) addContainsTest();
//addBufferResultValidatorTest();
  }

  void test() {
    try {
      for (Test test in nameToTestMap.values) {
        test.test();
      }
    } catch (e) {
      throw Exception(supplement(e.toString()) + "  " + e.toString());
    }
  }

  String supplement(String message) {
    String newMessage = "\n" + message + "\n";
    newMessage += "Original: " + wktWriter.writeFormatted(getOriginal()) + "\n";
    newMessage += "Buffer Distance: $bufferDistance\n";
    newMessage += "Buffer: " + wktWriter.writeFormatted(getBuffer()) + "\n";
    return newMessage.substring(0, newMessage.length - 1);
  }

  BufferValidator addTest(Test test) {
    nameToTestMap[test.getName()] = test;
    return this;
  }

  BufferValidator setExpectedArea(final double expectedArea) {
    var name = "Area Test";
    return addTest(new Test(name, () {
      double tolerance = (getBuffer().getArea() - getOriginal().buffer2(bufferDistance, QUADRANT_SEGMENTS_1 - QUADRANT_SEGMENTS_2).getArea()).abs();
      assertEqualsD(expectedArea, getBuffer().getArea(), tolerance, name);
    }));
  }

  BufferValidator setEmptyBufferExpected(final bool emptyBufferExpected) {
    return addTest(new Test("Empty Buffer Test", () {
      assertTrueMsg(supplement("Expected buffer " + (emptyBufferExpected ? "" : "not ") + "to be empty"), emptyBufferExpected == getBuffer().isEmpty());
    }, priority: 1));
  }

  BufferValidator setBufferHolesExpected(final bool bufferHolesExpected) {
    return addTest(new Test("Buffer Holes Test", () {
      assertTrueMsg(supplement("Expected buffer " + (bufferHolesExpected ? "" : "not ") + "to have holes"), hasHoles(getBuffer()) == bufferHolesExpected);
    }));
  }

  bool hasHoles(Geometry buffer) {
    if (buffer.isEmpty()) {
      return false;
    }
    if (buffer is Polygon) {
      return (buffer as Polygon).getNumInteriorRing() > 0;
    }
    MultiPolygon multiPolygon = buffer as MultiPolygon;
    for (int i = 0; i < multiPolygon.getNumGeometries(); i++) {
      if (hasHoles(multiPolygon.getGeometryN(i))) {
        return true;
      }
    }
    return false;
  }

  Geometry getOriginal() {
    if (original == null) {
      original = wktReader.read(wkt);
    }
    return original;
  }

  BufferValidator setPrecisionModel(PrecisionModel precisionModel) {
    wktReader = new WKTReader.withFactory(new GeometryFactory.withPrecisionModel(precisionModel));
    return this;
  }

  BufferValidator setFactory(PrecisionModel precisionModel, int srid) {
    wktReader = new WKTReader.withFactory(new GeometryFactory.withPrecisionModelSrid(precisionModel, srid));
    return this;
  }

  Geometry getBuffer() {
    if (buffer == null) {
      buffer = getOriginal().buffer2(bufferDistance, QUADRANT_SEGMENTS_1);
      if (getBuffer().runtimeType.toString() == "GeometryCollection" && getBuffer().isEmpty()) {
        try {
          //#contains doesn't work with GeometryLists [Jon Aquino
          // 10/29/2003]
          buffer = wktReader.read("POINT EMPTY");
        } catch (e) {
          Assert.shouldNeverReachHere();
        }
      }
    }
    return buffer;
  }

  void addContainsTest() {
    addTest(new Test("Contains Test", () {
      if (getOriginal().runtimeType.toString() == "GeometryCollection") {
        return;
      }
      Assert.isTrue(getOriginal().isValid());
      if (bufferDistance > 0) {
        assertTrueMsg(supplement("Expected buffer to contain original"), contains(getBuffer(), getOriginal()));
      } else {
        assertTrueMsg(supplement("Expected original to contain buffer"), contains(getOriginal(), getBuffer()));
      }
    }));
  }

  bool contains(Geometry a, Geometry b) {
    //JTS doesn't currently handle empty geometries correctly [Jon Aquino
    // 10/29/2003]
    if (b.isEmpty()) {
      return true;
    }
    bool isContained = a.contains(b);
    return isContained;
  }

  void addBufferResultValidatorTest() {
    addTest(new Test("BufferResultValidator Test", () {
      if (getOriginal().runtimeType.toString() == "GeometryCollection") {
        return;
      }
      assertTrueMsg(supplement("BufferResultValidator failure"), BufferResultValidator.isValidGDG(getOriginal(), bufferDistance, getBuffer()));
    }));
  }
}
