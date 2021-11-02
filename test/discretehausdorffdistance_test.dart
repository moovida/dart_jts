/*
 * Copyright (c) 2016 Vivid Solutions.
 *
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License 2.0
 * and Eclipse Distribution License v. 1.0 which accompanies this distribution.
 * The Eclipse Public License is available at http://www.eclipse.org/legal/epl-v20.html
 * and the Eclipse Distribution License is available at
 *
 * http://www.eclipse.org/org/documents/edl-v10.php.
 */

import 'package:dart_jts/dart_jts.dart';
import 'package:test/test.dart';

import 'testing_utilities.dart';

//region Utility functions
final double TOLERANCE = 0.00001;

void runTest(String wkt1, String wkt2, double expectedDistance){
  Geometry g1 = read(wkt1)!;
  Geometry g2 = read(wkt2)!;

  double distance = DiscreteHausdorffDistance.distanceStatic(g1, g2);
  expect(NumberUtils.equalsWithTolerance(expectedDistance, distance, TOLERANCE), true);
}

void runTestDF(String wkt1, String wkt2, double densifyFrac, double expectedDistance){
  Geometry g1 = read(wkt1)!;
  Geometry g2 = read(wkt2)!;

  double distance = DiscreteHausdorffDistance.distanceStaticDF(g1, g2, densifyFrac);
  expect(NumberUtils.equalsWithTolerance(expectedDistance, distance, TOLERANCE), true);
}
// endregion

// region tests

void testLineSegments(){
  runTest("LINESTRING (0 0, 2 1)", "LINESTRING (0 0, 2 0)", 1.0);
}
void testLineSegments2(){
  runTest("LINESTRING (0 0, 2 0)", "LINESTRING (0 1, 1 2, 2 1)", 2.0);
}
void testLinePoints(){
  runTest("LINESTRING (0 0, 2 0)", "MULTIPOINT (0 1, 1 0, 2 1)", 1.0);
}

///
/// Shows effects of limiting HD to vertices
/// Answer is not true Hausdorff distance.
///
void testLinesShowingDiscretenessEffect(){
  runTest("LINESTRING (130 0, 0 0, 0 150)", "LINESTRING (10 10, 10 150, 130 10)", 14.142135623730951);
}
void testLinesShowingDiscretenessEffectDF() {
  runTestDF(
      "LINESTRING (130 0, 0 0, 0 150)", "LINESTRING (10 10, 10 150, 130 10)",
      0.5, 70.0);
}

// endregion

void runAllTests(){
  test("testLineSegments", testLineSegments);
  test("testLineSegments2", testLineSegments2);
  test("testLinePoints", testLinePoints);
  test("testLinesShowingDiscretenessEffect", testLinesShowingDiscretenessEffect);
  test("testLinesShowingDiscretenessEffectDF", testLinesShowingDiscretenessEffectDF);
}

void main() {
  group("testDiscreteHausdorffDistance - ", runAllTests);
}