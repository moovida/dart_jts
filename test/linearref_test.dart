import 'package:dart_jts/dart_jts.dart';
import "package:test/test.dart";

import 'testing_utilities.dart';

/*
 * Copyright (c) 2016 Vivid Solutions.
 *
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse  License v1.0
 * and Eclipse Distribution License v. 1.0 which accompanies this distribution.
 * The Eclipse  License is available at http://www.eclipse.org/legal/epl-v10.html
 * and the Eclipse Distribution License is available at
 *
 * http://www.eclipse.org/org/documents/edl-v10.php.
 */

void main() {
  group("LocationIndexedLineTest - ", () {
    LocationIndexedLineTest testClass = LocationIndexedLineTest();

    test("testMultiLineStringSimple",
        () => testClass.testMultiLineStringSimple());
    test("testMultiLineString2", () => testClass.testMultiLineString2());
  });
  group("LinearLocationTest - ", () {
    LinearLocationTest testClass = LinearLocationTest();

    test(
        "testZeroLengthLineString", () => testClass.testZeroLengthLineString());
    test("testRepeatedCoordsLineString",
        () => testClass.testRepeatedCoordsLineString());
    test("testEndLocation", () => testClass.testEndLocation());
    test("testIsEndPoint", () => testClass.testIsEndPoint());
    test("testEndPointLowest", () => testClass.testEndPointLowest());
    test("testSameSegmentLineString",
        () => testClass.testSameSegmentLineString());
    test("testSameSegmentMultiLineString",
        () => testClass.testSameSegmentMultiLineString());
    test("testGetSegmentMultiLineString",
        () => testClass.testGetSegmentMultiLineString());
  });
  group("LengthIndexedLineTest - ", () {
    LengthIndexedLineTest testClass = LengthIndexedLineTest();

    test("testExtractLineBothIndicesAtEndpointXXX",
        () => testClass.testExtractLineBothIndicesAtEndpointXXX());
    test("testExtractLineBeyondRange",
        () => testClass.testExtractLineBeyondRange());
    test("testExtractLineReverse", () => testClass.testExtractLineReverse());

    // TODO check why this doesn work
    // test("testExtractLineReverseMulti",
    //     () => testClass.testExtractLineReverseMulti());
    test("testExtractLineNegative", () => testClass.testExtractLineNegative());
    test("testExtractLineNegativeReverse",
        () => testClass.testExtractLineNegativeReverse());
    test("testExtractLineIndexAtEndpoint",
        () => testClass.testExtractLineIndexAtEndpoint());
    test("testExtractLineIndexAtEndpointWithZeroLenComponents",
        () => testClass.testExtractLineIndexAtEndpointWithZeroLenComponents());
    test("testExtractLineBothIndicesAtEndpoint",
        () => testClass.testExtractLineBothIndicesAtEndpoint());
    test("testExtractLineBothIndicesAtEndpointNegative",
        () => testClass.testExtractLineBothIndicesAtEndpointNegative());
    test("testProjectExtractPoint", () => testClass.testProjectExtractPoint());
    test("testExtractPointBeyondRange",
        () => testClass.testExtractPointBeyondRange());
    test("testProjectPointWithDuplicateCoords",
        () => testClass.testProjectPointWithDuplicateCoords());
    test("testOffsetStartPointRepeatedPoint",
        () => testClass.testOffsetStartPointRepeatedPoint());
    test("testComputeZ", () => testClass.testComputeZ());

    // WKTReader is not good with 3d
    // test("testComputeZNaN", () => testClass.testComputeZNaN());
  });
}

/**
 * Tests the {@link LengthIndexedLine} class
 */
class LengthIndexedLineTest extends AbstractIndexedLineTest {
  void testExtractLineBothIndicesAtEndpointXXX() {
    checkExtractLine("MULTILINESTRING ((0 0, 10 0), (20 0, 25 0, 30 0))", -10,
        10, "LINESTRING (10 0, 10 0)");
  }

  void testExtractLineBeyondRange() {
    checkExtractLine(
        "LINESTRING (0 0, 10 10)", -100, 100, "LINESTRING (0 0, 10 10)");
  }

  void testExtractLineReverse() {
    checkExtractLine("LINESTRING (0 0, 10 0)", 9, 1, "LINESTRING (9 0, 1 0)");
  }

  void testExtractLineReverseMulti() {
    checkExtractLine("MULTILINESTRING ((0 0, 10 0), (20 0, 25 0, 30 0))", 19, 1,
        "MULTILINESTRING ((10 0, 1 0), (29 0, 25 0, 20 0))");
  }

  void testExtractLineNegative() {
    checkExtractLine("LINESTRING (0 0, 10 0)", -9, -1, "LINESTRING (1 0, 9 0)");
  }

  void testExtractLineNegativeReverse() {
    checkExtractLine("LINESTRING (0 0, 10 0)", -1, -9, "LINESTRING (9 0, 1 0)");
  }

  void testExtractLineIndexAtEndpoint() {
    checkExtractLine("MULTILINESTRING ((0 0, 10 0), (20 0, 25 0, 30 0))", 10,
        -1, "LINESTRING (20 0, 25 0, 29 0)");
  }

  /**
   * Tests that leading and trailing zero-length sublines are trimmed in the computed result,
   * and that zero-length extracts return the lowest extracted zero-length line
   */
  void testExtractLineIndexAtEndpointWithZeroLenComponents() {
    checkExtractLine(
        "MULTILINESTRING ((0 0, 10 0), (10 0, 10 0), (20 0, 25 0, 30 0))",
        10,
        -1,
        "LINESTRING (20 0, 25 0, 29 0)");
    checkExtractLine(
        "MULTILINESTRING ((0 0, 10 0), (10 0, 10 0), (20 0, 25 0, 30 0))",
        5,
        10,
        "LINESTRING (5 0, 10 0)");
    checkExtractLine(
        "MULTILINESTRING ((0 0, 10 0), (10 0, 10 0), (10 0, 10 0), (20 0, 25 0, 30 0))",
        10,
        10,
        "LINESTRING (10 0, 10 0)");
    checkExtractLine(
        "MULTILINESTRING ((0 0, 10 0), (10 0, 10 0), (10 0, 10 0), (10 0, 10 0), (20 0, 25 0, 30 0))",
        10,
        -10,
        "LINESTRING (10 0, 10 0)");
  }

  void testExtractLineBothIndicesAtEndpoint() {
    checkExtractLine("MULTILINESTRING ((0 0, 10 0), (20 0, 25 0, 30 0))", 10,
        10, "LINESTRING (10 0, 10 0)");
  }

  void testExtractLineBothIndicesAtEndpointNegative() {
    checkExtractLine("MULTILINESTRING ((0 0, 10 0), (20 0, 25 0, 30 0))", -10,
        10, "LINESTRING (10 0, 10 0)");
  }

  /**
   * From GEOS Ticket #323
   */
  void testProjectExtractPoint() {
    Geometry linearGeom = read("MULTILINESTRING ((0 2, 0 0), (-1 1, 1 1))");
    LengthIndexedLine indexedLine = new LengthIndexedLine(linearGeom);
    double index = indexedLine.project(new Coordinate(1, 0));
    Coordinate pt = indexedLine.extractPoint(index);
    assertTrue(pt.equals(new Coordinate(0, 0)));
  }

  void testExtractPointBeyondRange() {
    Geometry linearGeom = read("LINESTRING (0 0, 10 10)");
    LengthIndexedLine indexedLine = new LengthIndexedLine(linearGeom);
    Coordinate pt = indexedLine.extractPoint(100);
    assertTrue(pt.equals(new Coordinate(10, 10)));

    Coordinate pt2 = indexedLine.extractPoint(0);
    assertTrue(pt2.equals(new Coordinate(0, 0)));
  }

  void testProjectPointWithDuplicateCoords() {
    Geometry linearGeom = read("LINESTRING (0 0, 10 0, 10 0, 20 0)");
    LengthIndexedLine indexedLine = new LengthIndexedLine(linearGeom);
    double projIndex = indexedLine.project(new Coordinate(10, 1));
    assertTrue(projIndex == 10.0);
  }

  /**
   * These tests work for LengthIndexedLine, but not LocationIndexedLine
   *
   */
  void testOffsetStartPointRepeatedPoint() {
    runOffsetTest("LINESTRING (0 0, 10 10, 10 10, 20 20)", "POINT(0 0)", 1.0,
        "POINT (-0.7071067811865475 0.7071067811865475)");
    runOffsetTest("LINESTRING (0 0, 10 10, 10 10, 20 20)", "POINT(0 0)", -1.0,
        "POINT (0.7071067811865475 -0.7071067811865475)");
    runOffsetTest("LINESTRING (0 0, 10 10, 10 10, 20 20)", "POINT(10 10)", 5.0,
        "POINT (6.464466094067262 13.535533905932738)");
    runOffsetTest("LINESTRING (0 0, 10 10, 10 10, 20 20)", "POINT(10 10)", -5.0,
        "POINT (13.535533905932738 6.464466094067262)");
  }

  /**
   * Tests that z values are interpolated
   *
   */
  void testComputeZ() {
    Geometry linearGeom = read("LINESTRING Z(0 0 0, 10 10 10)");
    LengthIndexedLine indexedLine = new LengthIndexedLine(linearGeom);
    double projIndex = indexedLine.project(new Coordinate(5, 5));
    Coordinate projPt = indexedLine.extractPoint(projIndex);
//    System.out.println(projPt);
    assertTrue(projPt.equals3D(new Coordinate.fromXYZ(5, 5, 5)));
  }

  /**
   * Tests that if the input does not have Z ordinates, neither does the output.
   *
   */
  void testComputeZNaN() {
    Geometry linearGeom = read("LINESTRING Z(0 0, 10 10 10)");
    LengthIndexedLine indexedLine = new LengthIndexedLine(linearGeom);
    double projIndex = indexedLine.project(new Coordinate(5, 5));
    Coordinate projPt = indexedLine.extractPoint(projIndex);
    assertTrue(projPt.getZ().isNaN);
  }

  void checkExtractLine(String wkt, double start, double end, String expected) {
    Geometry linearGeom = read(wkt);
    LengthIndexedLine indexedLine = new LengthIndexedLine(linearGeom);
    Geometry result = indexedLine.extractLine(start, end);
    checkExpected(result, expected);
  }

  Geometry indicesOfThenExtract(Geometry linearGeom, Geometry subLine) {
    LengthIndexedLine indexedLine = new LengthIndexedLine(linearGeom);
    List<double> loc = indexedLine.indicesOf(subLine);
    Geometry result = indexedLine.extractLine(loc[0], loc[1]);
    return result;
  }

  bool indexOfAfterCheck(Geometry linearGeom, Coordinate testPt) {
    LengthIndexedLine indexedLine = new LengthIndexedLine(linearGeom);

    // check locations are consecutive
    double loc1 = indexedLine.indexOf(testPt);
    double loc2 = indexedLine.indexOfAfter(testPt, loc1);
    if (loc2 <= loc1) return false;

    // check extracted points are the same as the input
    Coordinate pt1 = indexedLine.extractPoint(loc1);
    Coordinate pt2 = indexedLine.extractPoint(loc2);
    if (!pt1.equals2D(testPt)) return false;
    if (!pt2.equals2D(testPt)) return false;

    return true;
  }

  bool indexOfAfterCheckWithAfterCoord(
      Geometry linearGeom, Coordinate testPt, Coordinate checkPt) {
    LengthIndexedLine indexedLine = new LengthIndexedLine(linearGeom);

    // check that computed location is after check location
    double checkLoc = indexedLine.indexOf(checkPt);
    double testLoc = indexedLine.indexOfAfter(testPt, checkLoc);
    if (testLoc < checkLoc) return false;

    return true;
  }

  Coordinate extractOffsetAt(
      Geometry linearGeom, Coordinate testPt, double offsetDistance) {
    LengthIndexedLine indexedLine = new LengthIndexedLine(linearGeom);
    double index = indexedLine.indexOf(testPt);
    return indexedLine.extractPointWithOffset(index, offsetDistance);
  }
}

/**
 * Tests methods involving only {@link LinearLocation}s
 * 
 * @author Martin Davis
 *
 */
class LinearLocationTest {
  WKTReader reader = new WKTReader();

  void testZeroLengthLineString() {
    Geometry line = reader.read("LINESTRING (10 0, 10 0)");
    LocationIndexedLine indexedLine = new LocationIndexedLine(line);
    LinearLocation loc0 = indexedLine.indexOf(new Coordinate(11, 0));
    assertTrue(
        loc0.compareTo(new LinearLocation.fromSegmentIndexFraction(1, 0.0)) ==
            0);
  }

  void testRepeatedCoordsLineString() {
    Geometry line = reader.read("LINESTRING (10 0, 10 0, 20 0)");
    LocationIndexedLine indexedLine = new LocationIndexedLine(line);
    LinearLocation loc0 = indexedLine.indexOf(new Coordinate(11, 0));
    assertTrue(
        loc0.compareTo(new LinearLocation.fromSegmentIndexFraction(1, 0.1)) ==
            0);
  }

  void testEndLocation() {
    Geometry line = reader.read("LINESTRING (10 0, 20 0)");
    LinearLocation loc0 = LinearLocation.getEndLocation(line);
    assertTrue(0 == loc0.getSegmentFraction());
    assertTrue(1 == loc0.getSegmentIndex());

    LocationIndexedLine indexedLine = new LocationIndexedLine(line);
    LinearLocation endLoc = indexedLine.getEndIndex();
    LinearLocation normLoc =
        new LinearLocation.fromComponentSegmentIndexFraction(
            endLoc.getComponentIndex(),
            endLoc.getSegmentIndex(),
            endLoc.getSegmentFraction());
    assertTrue(normLoc.getComponentIndex() == endLoc.getComponentIndex());
    assertTrue(normLoc.getSegmentIndex() == endLoc.getSegmentIndex());
    assertTrue(normLoc.getSegmentFraction() == endLoc.getSegmentFraction());
  }

  void testIsEndPoint() {
    Geometry line = reader.read("LINESTRING (10 0, 20 0)");

    assertTrue(
        !(new LinearLocation.fromSegmentIndexFraction(0, 0)).isEndpoint(line));
    assertTrue(!(new LinearLocation.fromSegmentIndexFraction(0, 0.5))
        .isEndpoint(line));
    assertTrue(!(new LinearLocation.fromSegmentIndexFraction(0, 0.9999))
        .isEndpoint(line));

    assertTrue(
        (new LinearLocation.fromSegmentIndexFraction(0, 1.0)).isEndpoint(line));

    assertTrue(
        (new LinearLocation.fromSegmentIndexFraction(1, 0.0)).isEndpoint(line));
    assertTrue(
        (new LinearLocation.fromSegmentIndexFraction(1, 0.5)).isEndpoint(line));
    assertTrue(
        (new LinearLocation.fromSegmentIndexFraction(1, 1.0)).isEndpoint(line));
    assertTrue(
        (new LinearLocation.fromSegmentIndexFraction(1, 1.5)).isEndpoint(line));

    assertTrue(
        (new LinearLocation.fromSegmentIndexFraction(2, 0.5)).isEndpoint(line));

    LinearLocation loc = new LinearLocation.fromSegmentIndexFraction(0, 0.0);
    loc.setToEnd(line);
    assertTrue(loc.isEndpoint(line));

    LinearLocation locLow = loc.toLowest(line);
    assertTrue(locLow.isEndpoint(line));
  }

  void testEndPointLowest() {
    Geometry line = reader.read("LINESTRING (10 0, 20 0, 30 10)");

    assertTrue(
        (new LinearLocation.fromSegmentIndexFraction(1, 1.0)).isEndpoint(line));
    assertTrue(
        (new LinearLocation.fromSegmentIndexFraction(2, 0.0)).isEndpoint(line));
    assertTrue(
        (new LinearLocation.fromSegmentIndexFraction(2, 0.5)).isEndpoint(line));

    LinearLocation loc = new LinearLocation.fromSegmentIndexFraction(0, 0.0);
    loc.setToEnd(line);
    assertTrue(loc.isEndpoint(line));
    assertEquals(2, loc.getSegmentIndex());
    assertEquals(0.0, loc.getSegmentFraction());

    LinearLocation locLow = loc.toLowest(line);
    assertTrue(locLow.isEndpoint(line));
    assertEquals(1, locLow.getSegmentIndex());
    assertEquals(1.0, locLow.getSegmentFraction());
  }

  void testSameSegmentLineString() {
    Geometry line = reader.read("LINESTRING (0 0, 10 0, 20 0, 30 0)");
    LocationIndexedLine indexedLine = new LocationIndexedLine(line);

    LinearLocation loc0 = indexedLine.indexOf(new Coordinate(0, 0));
    LinearLocation loc0_5 = indexedLine.indexOf(new Coordinate(5, 0));
    LinearLocation loc1 = indexedLine.indexOf(new Coordinate(10, 0));
    LinearLocation loc2 = indexedLine.indexOf(new Coordinate(20, 0));
    LinearLocation loc2_5 = indexedLine.indexOf(new Coordinate(25, 0));
    LinearLocation loc3 = indexedLine.indexOf(new Coordinate(30, 0));

    assertTrue(loc0.isOnSameSegment(loc0));
    assertTrue(loc0.isOnSameSegment(loc0_5));
    assertTrue(loc0.isOnSameSegment(loc1));
    assertTrue(!loc0.isOnSameSegment(loc2));
    assertTrue(!loc0.isOnSameSegment(loc2_5));
    assertTrue(!loc0.isOnSameSegment(loc3));

    assertTrue(loc0_5.isOnSameSegment(loc0));
    assertTrue(loc0_5.isOnSameSegment(loc1));
    assertTrue(!loc0_5.isOnSameSegment(loc2));
    assertTrue(!loc0_5.isOnSameSegment(loc3));

    assertTrue(!loc2.isOnSameSegment(loc0));
    assertTrue(loc2.isOnSameSegment(loc1));
    assertTrue(loc2.isOnSameSegment(loc2));
    assertTrue(loc2.isOnSameSegment(loc3));

    assertTrue(loc2_5.isOnSameSegment(loc3));

    assertTrue(!loc3.isOnSameSegment(loc0));
    assertTrue(loc3.isOnSameSegment(loc2));
    assertTrue(loc3.isOnSameSegment(loc2_5));
    assertTrue(loc3.isOnSameSegment(loc3));
  }

  void testSameSegmentMultiLineString() {
    Geometry line =
        reader.read("MULTILINESTRING ((0 0, 10 0, 20 0), (20 0, 30 0))");
    LocationIndexedLine indexedLine = new LocationIndexedLine(line);

    LinearLocation loc0 = indexedLine.indexOf(new Coordinate(0, 0));
    LinearLocation loc0_5 = indexedLine.indexOf(new Coordinate(5, 0));
    LinearLocation loc1 = indexedLine.indexOf(new Coordinate(10, 0));
    LinearLocation loc2 = indexedLine.indexOf(new Coordinate(20, 0));
    LinearLocation loc2B =
        new LinearLocation.fromComponentSegmentIndexFraction(1, 0, 0.0);

    LinearLocation loc2_5 = indexedLine.indexOf(new Coordinate(25, 0));
    LinearLocation loc3 = indexedLine.indexOf(new Coordinate(30, 0));

    assertTrue(loc0.isOnSameSegment(loc0));
    assertTrue(loc0.isOnSameSegment(loc0_5));
    assertTrue(loc0.isOnSameSegment(loc1));
    assertTrue(!loc0.isOnSameSegment(loc2));
    assertTrue(!loc0.isOnSameSegment(loc2_5));
    assertTrue(!loc0.isOnSameSegment(loc3));

    assertTrue(loc0_5.isOnSameSegment(loc0));
    assertTrue(loc0_5.isOnSameSegment(loc1));
    assertTrue(!loc0_5.isOnSameSegment(loc2));
    assertTrue(!loc0_5.isOnSameSegment(loc3));

    assertTrue(!loc2.isOnSameSegment(loc0));
    assertTrue(loc2.isOnSameSegment(loc1));
    assertTrue(loc2.isOnSameSegment(loc2));
    assertTrue(!loc2.isOnSameSegment(loc3));
    assertTrue(loc2B.isOnSameSegment(loc3));

    assertTrue(loc2_5.isOnSameSegment(loc3));

    assertTrue(!loc3.isOnSameSegment(loc0));
    assertTrue(!loc3.isOnSameSegment(loc2));
    assertTrue(loc3.isOnSameSegment(loc2B));
    assertTrue(loc3.isOnSameSegment(loc2_5));
    assertTrue(loc3.isOnSameSegment(loc3));
  }

  void testGetSegmentMultiLineString() {
    Geometry line =
        reader.read("MULTILINESTRING ((0 0, 10 0, 20 0), (20 0, 30 0))");
    LocationIndexedLine indexedLine = new LocationIndexedLine(line);

    LinearLocation loc0 = indexedLine.indexOf(new Coordinate(0, 0));
    LinearLocation loc0_5 = indexedLine.indexOf(new Coordinate(5, 0));
    LinearLocation loc1 = indexedLine.indexOf(new Coordinate(10, 0));
    LinearLocation loc2 = indexedLine.indexOf(new Coordinate(20, 0));
    LinearLocation loc2B =
        new LinearLocation.fromComponentSegmentIndexFraction(1, 0, 0.0);

    LinearLocation loc2_5 = indexedLine.indexOf(new Coordinate(25, 0));
    LinearLocation loc3 = indexedLine.indexOf(new Coordinate(30, 0));

    LineSegment seg0 = new LineSegment.fromCoordinates(
        new Coordinate(0, 0), new Coordinate(10, 0));
    LineSegment seg1 = new LineSegment.fromCoordinates(
        new Coordinate(10, 0), new Coordinate(20, 0));
    LineSegment seg2 = new LineSegment.fromCoordinates(
        new Coordinate(20, 0), new Coordinate(30, 0));

    assertTrue(loc0.getSegment(line).equals(seg0));
    assertTrue(loc0_5.getSegment(line).equals(seg0));

    assertTrue(loc1.getSegment(line).equals(seg1));
    assertTrue(loc2.getSegment(line).equals(seg1));

    assertTrue(loc2_5.getSegment(line).equals(seg2));
    assertTrue(loc3.getSegment(line).equals(seg2));
  }
}

class LocationIndexedLineTest extends AbstractIndexedLineTest {
  void testMultiLineStringSimple() {
    runExtractLine(
        "MULTILINESTRING ((0 0, 10 10), (20 20, 30 30))",
        new LinearLocation.fromComponentSegmentIndexFraction(0, 0, .5),
        new LinearLocation.fromComponentSegmentIndexFraction(1, 0, .5),
        "MULTILINESTRING ((5 5, 10 10), (20 20, 25 25))");
  }

  void testMultiLineString2() {
    runExtractLine(
        "MULTILINESTRING ((0 0, 10 10), (20 20, 30 30))",
        new LinearLocation.fromComponentSegmentIndexFraction(0, 0, 1.0),
        new LinearLocation.fromComponentSegmentIndexFraction(1, 0, .5),
        "MULTILINESTRING ((10 10, 10 10), (20 20, 25 25))");
  }

  void runExtractLine(
      String wkt, LinearLocation start, LinearLocation end, String expected) {
    Geometry geom = read(wkt);
    LocationIndexedLine lil = new LocationIndexedLine(geom);
    Geometry result = lil.extractLine(start, end);
    //System.out.println(result);
    checkExpected(result, expected);
  }

  Geometry indicesOfThenExtract(Geometry input, Geometry subLine) {
    LocationIndexedLine indexedLine = new LocationIndexedLine(input);
    List<LinearLocation> loc = indexedLine.indicesOf(subLine);
    Geometry result = indexedLine.extractLine(loc[0], loc[1]);
    return result;
  }

  bool indexOfAfterCheck(Geometry linearGeom, Coordinate testPt) {
    LocationIndexedLine indexedLine = new LocationIndexedLine(linearGeom);

    // check locations are consecutive
    LinearLocation loc1 = indexedLine.indexOf(testPt);
    LinearLocation loc2 = indexedLine.indexOfAfter(testPt, loc1);
    if (loc2.compareTo(loc1) <= 0) return false;

    // check extracted points are the same as the input
    Coordinate pt1 = indexedLine.extractPoint(loc1);
    Coordinate pt2 = indexedLine.extractPoint(loc2);
    if (!pt1.equals2D(testPt)) return false;
    if (!pt2.equals2D(testPt)) return false;

    return true;
  }

  bool indexOfAfterCheckWithAfterCoord(
      Geometry linearGeom, Coordinate testPt, Coordinate afterPt) {
    LocationIndexedLine indexedLine = new LocationIndexedLine(linearGeom);

    // check that computed location is after check location
    LinearLocation afterLoc = indexedLine.indexOf(afterPt);
    LinearLocation testLoc = indexedLine.indexOfAfter(testPt, afterLoc);
    if (testLoc.compareTo(afterLoc) < 0) return false;

    return true;
  }

  Coordinate extractOffsetAt(
      Geometry linearGeom, Coordinate testPt, double offsetDistance) {
    LocationIndexedLine indexedLine = new LocationIndexedLine(linearGeom);
    LinearLocation index = indexedLine.indexOf(testPt);
    return indexedLine.extractPointWithOffset(index, offsetDistance);
  }
}

/**
 * Base class for linear referencing class unit tests.
 */
abstract class AbstractIndexedLineTest {
  WKTReader reader = new WKTReader();

  void testFirst() {
    runOffsetTest(
        "LINESTRING (0 0, 20 20)", "POINT(20 20)", 0.0, "POINT (20 20)");
  }

  void testML() {
    runIndicesOfThenExtract("MULTILINESTRING ((0 0, 10 10), (20 20, 30 30))",
        "MULTILINESTRING ((1 1, 10 10), (20 20, 25 25))");
  }

  void testPartOfSegmentNoVertex() {
    runIndicesOfThenExtract(
        "LINESTRING (0 0, 10 10, 20 20)", "LINESTRING (1 1, 9 9)");
  }

  void testPartOfSegmentContainingVertex() {
    runIndicesOfThenExtract(
        "LINESTRING (0 0, 10 10, 20 20)", "LINESTRING (5 5, 10 10, 15 15)");
  }

  /**
   * Tests that duplicate coordinates are handled correctly.
   */
  void testPartOfSegmentContainingDuplicateCoords() {
    runIndicesOfThenExtract("LINESTRING (0 0, 10 10, 10 10, 20 20)",
        "LINESTRING (5 5, 10 10, 10 10, 15 15)");
  }

  /**
   * Following tests check that correct portion of loop is identified.
   * This requires that the correct vertex for (0,0) is selected.
   */

  void testLoopWithStartSubLine() {
    runIndicesOfThenExtract("LINESTRING (0 0, 0 10, 10 10, 10 0, 0 0)",
        "LINESTRING (0 0, 0 10, 10 10)");
  }

  void testLoopWithEndingSubLine() {
    runIndicesOfThenExtract("LINESTRING (0 0, 0 10, 10 10, 10 0, 0 0)",
        "LINESTRING (10 10, 10 0, 0 0)");
  }

  // test a subline equal to the parent loop
  void testLoopWithIdenticalSubLine() {
    runIndicesOfThenExtract("LINESTRING (0 0, 0 10, 10 10, 10 0, 0 0)",
        "LINESTRING (0 0, 0 10, 10 10, 10 0, 0 0)");
  }

  // test a zero-length subline equal to the start point
  void testZeroLenSubLineAtStart() {
    runIndicesOfThenExtract(
        "LINESTRING (0 0, 0 10, 10 10, 10 0, 0 0)", "LINESTRING (0 0, 0 0)");
  }

  // test a zero-length subline equal to a mid point
  void testZeroLenSubLineAtMidVertex() {
    runIndicesOfThenExtract("LINESTRING (0 0, 0 10, 10 10, 10 0, 0 0)",
        "LINESTRING (10 10, 10 10)");
  }

  void testIndexOfAfterSquare() {
    runIndexOfAfterTest(
        "LINESTRING (0 0, 0 10, 10 10, 10 0, 0 0)", "POINT (0 0)");
  }

  void testIndexOfAfterRibbon() {
    runIndexOfAfterTest(
        "LINESTRING (0 0, 0 60, 50 60, 50 20, -20 20)", "POINT (0 20)");
    runIndexOfAfterTestWithAfterWkt(
        "LINESTRING (0 0, 0 60, 50 60, 50 20, -20 20)",
        "POINT (0 20)",
        "POINT (30 60)");
  }

  void testIndexOfAfterBeyondEndRibbon() {
    runIndexOfAfterTestWithAfterWkt(
        "LINESTRING (0 0, 0 60, 50 60, 50 20, -20 20)",
        "POINT (-30 20)",
        "POINT (-20 20)");
  }

  void testOffsetStartPoint() {
    runOffsetTest("LINESTRING (0 0, 10 10, 20 20)", "POINT(0 0)", 1.0,
        "POINT (-0.7071067811865475 0.7071067811865475)");
    runOffsetTest("LINESTRING (0 0, 10 10, 20 20)", "POINT(0 0)", -1.0,
        "POINT (0.7071067811865475 -0.7071067811865475)");
    runOffsetTest("LINESTRING (0 0, 10 10, 20 20)", "POINT(10 10)", 5.0,
        "POINT (6.464466094067262 13.535533905932738)");
    runOffsetTest("LINESTRING (0 0, 10 10, 20 20)", "POINT(10 10)", -5.0,
        "POINT (13.535533905932738 6.464466094067262)");
  }

  void testOffsetStartPointRepeatedPoint() {
    runOffsetTest("LINESTRING (0 0, 10 10, 10 10, 20 20)", "POINT(0 0)", 1.0,
        "POINT (-0.7071067811865475 0.7071067811865475)");
    runOffsetTest("LINESTRING (0 0, 10 10, 10 10, 20 20)", "POINT(0 0)", -1.0,
        "POINT (0.7071067811865475 -0.7071067811865475)");
    // These tests work for LengthIndexedLine, but not LocationIndexedLine
    //runOffsetTest("LINESTRING (0 0, 10 10, 10 10, 20 20)", "POINT(10 10)", 5.0, "POINT (6.464466094067262 13.535533905932738)");
    //runOffsetTest("LINESTRING (0 0, 10 10, 10 10, 20 20)", "POINT(10 10)", -5.0, "POINT (13.535533905932738 6.464466094067262)");
  }

  void testOffsetEndPoint() {
    runOffsetTest(
        "LINESTRING (0 0, 20 20)", "POINT(20 20)", 0.0, "POINT (20 20)");
    runOffsetTest(
        "LINESTRING (0 0, 13 13, 20 20)", "POINT(20 20)", 0.0, "POINT (20 20)");
    runOffsetTest(
        "LINESTRING (0 0, 10 0, 20 0)", "POINT(20 0)", 1.0, "POINT (20 1)");
    runOffsetTest("LINESTRING (0 0, 20 0)", "POINT(10 0)", 1.0,
        "POINT (10 1)"); // point on last segment
    runOffsetTest("MULTILINESTRING ((0 0, 10 0), (10 0, 20 0))", "POINT(10 0)",
        -1.0, "POINT (10 -1)");
    runOffsetTest("MULTILINESTRING ((0 0, 10 0), (10 0, 20 0))", "POINT(20 0)",
        1.0, "POINT (20 1)");
  }

  Geometry read(String wkt) {
    try {
      return reader.read(wkt);
    } catch (ex) {
      throw new RuntimeException(ex);
    }
  }

  void runIndicesOfThenExtract(String inputStr, String subLineStr)
//
  {
    Geometry input = read(inputStr);
    Geometry subLine = read(subLineStr);
    Geometry result = indicesOfThenExtract(input, subLine);
    checkExpected(result, subLineStr);
  }

  void checkExpected(Geometry result, String expected) {
    Geometry subLine = read(expected);
    bool isEqual = result.equalsExactWithTol(subLine, 1.0e-5);
    if (!isEqual) {
      print("Computed result is: $result");
    }
    assertTrue(isEqual);
  }

  Geometry indicesOfThenExtract(Geometry input, Geometry subLine);

/*
  // example of indicesOfThenLocate method
   Geometry indicesOfThenLocate(LineString input, LineString subLine)
  {
    LocationIndexedLine indexedLine = new LocationIndexedLine(input);
    LineStringLocation[] loc = indexedLine.indicesOf(subLine);
    Geometry result = indexedLine.locate(loc[0], loc[1]);
    return result;
  }
*/

  void runIndexOfAfterTest(String inputStr, String testPtWKT)
//
  {
    Geometry input = read(inputStr);
    Geometry testPoint = read(testPtWKT);
    Coordinate testPt = testPoint.getCoordinate();
    bool resultOK = indexOfAfterCheck(input, testPt);
    assertTrue(resultOK);
  }

  void runIndexOfAfterTestWithAfterWkt(
      String inputStr, String testPtWKT, String afterPtWKT)
//
  {
    Geometry input = read(inputStr);
    Geometry testPoint = read(testPtWKT);
    Coordinate testPt = testPoint.getCoordinate();
    Geometry afterPoint = read(afterPtWKT);
    Coordinate afterPt = afterPoint.getCoordinate();
    bool resultOK = indexOfAfterCheckWithAfterCoord(input, testPt, afterPt);
    assertTrue(resultOK);
  }

  /**
   * Checks that the point computed by <tt>indexOfAfter</tt>
   * is the same as the input point.
   * (This should be the case for all except pathological cases, 
   * such as the input test point being beyond the end of the line). 
   * 
   * @param input
   * @param testPt
   * @return true if the result of indexOfAfter is the same as the input point
   */
  bool indexOfAfterCheck(Geometry input, Coordinate testPt);

  bool indexOfAfterCheckWithAfterCoord(
      Geometry input, Coordinate testPt, Coordinate afterPt);

  static final double TOLERANCE_DIST = 0.001;

  void runOffsetTest(String inputWKT, String testPtWKT, double offsetDistance,
      String expectedPtWKT)
//
  {
    Geometry input = read(inputWKT);
    Geometry testPoint = read(testPtWKT);
    Geometry expectedPoint = read(expectedPtWKT);
    Coordinate testPt = testPoint.getCoordinate();
    Coordinate expectedPt = expectedPoint.getCoordinate();
    Coordinate offsetPt = extractOffsetAt(input, testPt, offsetDistance);

    bool isOk = offsetPt.distance(expectedPt) < TOLERANCE_DIST;
    if (!isOk) print("Expected = $expectedPoint  Actual = $offsetPt");
    assertTrue(isOk);
  }

  Coordinate extractOffsetAt(
      Geometry input, Coordinate testPt, double offsetDistance);
}
