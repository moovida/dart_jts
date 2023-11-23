import 'package:dart_jts/dart_jts.dart';
import 'package:test/test.dart';
import 'spatial_index_tester.dart';

void main() {
  test('testSpatialIndex', () {
    SpatialIndexTester tester = new SpatialIndexTester();
    tester.setSpatialIndex(new Quadtree());
    tester.init();
    tester.run();
    expect(tester.isSuccess, true);
  });
  test('testNullQuery', () {
    Quadtree qt = new Quadtree();
    List result1 = qt.query(null);
    expect(result1.length == 0, true);

    qt.insert(new Envelope(0, 10, 0, 10), "some data");
    List result2 = qt.query(null);
    expect(result2.length == 0, true);
  });
}
