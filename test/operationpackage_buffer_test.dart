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
    test("tes1", () {});
    test("tes1", () {});
    test("tes1", () {});
    test("tes1", () {});
    test("tes1", () {});
    test("tes1", () {});
    test("tes1", () {});
    test("tes1", () {});
    test("tes1", () {});
    test("tes1", () {});
    test("tes1", () {});
    test("tes1", () {});
    test("tes1", () {});
    test("tes1", () {});
    test("tes1", () {});
    test("tes1", () {});
    test("tes1", () {});
    test("tes1", () {});
    test("tes1", () {});
    test("tes1", () {});
    test("tes1", () {});
    test("tes1", () {});
    test("tes1", () {});
    test("tes1", () {});
    test("tes1", () {});
    test("tes1", () {});
    test("tes1", () {});
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
