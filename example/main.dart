import 'package:dart_jts/dart_jts.dart';

void main() {
  String polygon1Wkt = "POLYGON (( 0 0, 100 0, 100 100, 0 100, 0 0 ))";
  String polygon2Wkt = "POLYGON (( 50 50, 150 50, 150 150, 50 150, 50 50 ))";

  WKTReader reader = WKTReader();

  var p1 = reader.read(polygon1Wkt);
  var p2 = reader.read(polygon2Wkt);

  if (p1.intersects(p2)) {
    print("Yeah!");
  } else {
    print("Something's very wrong!");
  }

  print("For many many many examples, have a look at the testcases!");
}
