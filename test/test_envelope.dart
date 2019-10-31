import "package:test/test.dart";
import 'package:dart_sfs/dart_sfs.dart';

void main() {
  group("envelope - ", () {
    test("testEverything", () {
      Envelope e1 = Envelope.empty();
      expect(e1.isNull(), true);
      expect(0, e1.getWidth());
      expect(0, e1.getHeight());
      e1.expandToInclude(100, 101);
      e1.expandToInclude(200, 202);
      e1.expandToInclude(150, 151);
      expect(200, e1.getMaxX());
      expect(202, e1.getMaxY());
      expect(100, e1.getMinX());
      expect(101, e1.getMinY());
      expect(e1.contains(120, 120), true);
      expect(e1.contains(120, 101), true);
      expect(!e1.contains(120, 100), true);
      expect(101, e1.getHeight());
      expect(100, e1.getWidth());
      expect(!e1.isNull(), true);

      Envelope e2 =  Envelope(499, 500, 500, 501);
      expect(!e1.containsEnvelope(e2), true);
      expect(!e1.intersectsEnvelope(e2), true);
      e1.expandToIncludeEnvelope(e2);
      expect(e1.containsEnvelope(e2), true);
      expect(e1.intersectsEnvelope(e2), true);
      expect(500, e1.getMaxX());
      expect(501, e1.getMaxY());
      expect(100, e1.getMinX());
      expect(101, e1.getMinY());

      Envelope e3 =  Envelope(300, 700, 300, 700);
      expect(!e1.containsEnvelope(e3), true);
      expect(e1.intersectsEnvelope(e3), true);

      Envelope e4 =  Envelope(300, 301, 300, 301);
      expect(e1.containsEnvelope(e4), true);
      expect(e1.intersectsEnvelope(e4), true);
    });
    test("empty", () {
    });
  });

}
