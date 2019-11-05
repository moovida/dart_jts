import 'package:dart_sfs/dart_sfs.dart';
import 'package:test/test.dart';

assertEquals(actual, matcher) {
  if (actual is double && matcher is double) {
    if (actual.isNaN && matcher.isNaN) return;
  }
  if (actual is Coordinate && matcher is Coordinate) {
    if (actual.equals(matcher)) return;
  }
  expect(actual, matcher);
}

assertEqualsD(double n1, double n2, double tolerance) {
  expect(NumberUtils.equalsWithTolerance(n1, n2, tolerance), true);
}

assertTrue(actual) {
  expect(actual, true);
}

assertTrueMsg(String msg, actual) {
  expect(actual, true, reason: msg);
}
