import 'package:dart_jts/dart_jts.dart';

class SpatialIndexTester {
  static bool VERBOSE = false;
  static final double CELL_EXTENT = 20.31;
  static final int CELLS_PER_GRID_SIDE = 10;
  static final double FEATURE_EXTENT = 10.1;
  static final double OFFSET = 5.03;
  static final double QUERY_ENVELOPE_EXTENT_1 = 1.009;
  static final double QUERY_ENVELOPE_EXTENT_2 = 11.7;

  SpatialIndex? index;
  List<Envelope> sourceData = List.empty(growable: true);
  bool isSuccess = true;

  SpatialIndexTester() {}

  void setSpatialIndex(SpatialIndex index) {
    this.index = index;
  }

  SpatialIndex? getSpatialIndex() {
    return index;
  }

  void init() {
    sourceData = List.empty(growable: true);
    addSourceData(0, sourceData);
    addSourceData(OFFSET, sourceData);
    if (VERBOSE) {
      //System.out.println("===============================");
      //System.out.println("Grid Extent: " + (CELL_EXTENT * CELLS_PER_GRID_SIDE));
      //System.out.println("Cell Extent: " + CELL_EXTENT);
      //System.out.println("Feature Extent: " + FEATURE_EXTENT);
      //System.out.println("Cells Per Grid Side: " + CELLS_PER_GRID_SIDE);
      //System.out.println("Offset For 2nd Set Of Features: " + OFFSET);
      //System.out.println("Feature Count: " + sourceData.size());
    }
    insert(sourceData, index!);
  }

  void run() {
    doTest(index!, QUERY_ENVELOPE_EXTENT_1, sourceData);
    doTest(index!, QUERY_ENVELOPE_EXTENT_2, sourceData);
  }

  void insert(List<Envelope> sourceData, SpatialIndex index) {
    for (var envelope in sourceData) {
      index.insert(envelope, envelope);
    }
  }

  void addSourceData(double offset, List<Envelope> sourceData) {
    for (int i = 0; i < CELLS_PER_GRID_SIDE; i++) {
      double minx = (i * CELL_EXTENT) + offset;
      double maxx = minx + FEATURE_EXTENT;
      for (int j = 0; j < CELLS_PER_GRID_SIDE; j++) {
        double miny = (j * CELL_EXTENT) + offset;
        double maxy = miny + FEATURE_EXTENT;
        Envelope e = Envelope(minx, maxx, miny, maxy);
        sourceData.add(e);
      }
    }
  }

  void doTest(SpatialIndex index, double queryEnvelopeExtent,
      List<Envelope> sourceData) {
    // int extraMatchCount = 0;
    // int expectedMatchCount = 0;
    // int actualMatchCount = 0;
    // int queryCount = 0;
    for (var x = 0.0;
        x < CELL_EXTENT * CELLS_PER_GRID_SIDE;
        x += queryEnvelopeExtent) {
      for (var y = 0.0;
          y < CELL_EXTENT * CELLS_PER_GRID_SIDE;
          y += queryEnvelopeExtent) {
        Envelope queryEnvelope =
            Envelope(x, x + queryEnvelopeExtent, y, y + queryEnvelopeExtent);
        List<Envelope> expectedMatches =
            intersectingEnvelopes(queryEnvelope, sourceData);
        List<Envelope> actualMatches =
            index.query(queryEnvelope).cast<Envelope>();
        // since index returns candidates only, it may return more than the expected value
        if (expectedMatches.length > actualMatches.length) {
          isSuccess = false;
        }
        // extraMatchCount += (actualMatches.length - expectedMatches.length);
        // expectedMatchCount += expectedMatches.length;
        // actualMatchCount += actualMatches.length;
        // compare(expectedMatches, actualMatches);
        // queryCount++;
      }
    }
    if (VERBOSE) {
      //System.out.println("---------------");
      //System.out.println("Envelope Extent: " + queryEnvelopeExtent);
      //System.out.println("Expected Matches: " + expectedMatchCount);
      //System.out.println("Actual Matches: " + actualMatchCount);
      //System.out.println("Extra Matches: " + extraMatchCount);
      //System.out.println("Query Count: " + queryCount);
      //System.out.println("Average Expected Matches: " + (expectedMatchCount/(double)queryCount));
      //System.out.println("Average Actual Matches: " + (actualMatchCount/(double)queryCount));
      //System.out.println("Average Extra Matches: " + (extraMatchCount/(double)queryCount));
    }
  }

  void compare(
      List<Envelope> expectedEnvelopes, List<Envelope> actualEnvelopes) {
    //Don't use #containsAll because we want to check using
    //==, not #equals. [Jon Aquino]
    for (var expected in expectedEnvelopes) {
      bool found = false;
      for (var actual in actualEnvelopes) {
        if (actual == expected) {
          found = true;
          break;
        }
      }
      if (!found) isSuccess = false;
    }
  }

  List<Envelope> intersectingEnvelopes(
      Envelope queryEnvelope, List<Envelope> envelopes) {
    List<Envelope> intersectingEnvelopes = List.empty(growable: true);
    for (var candidate in envelopes) {
      if (candidate.intersectsEnvelope(queryEnvelope)) {
        intersectingEnvelopes.add(candidate);
      }
    }
    return intersectingEnvelopes;
  }
}
