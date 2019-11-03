import 'package:dart_sfs/dart_sfs.dart';
import "package:test/test.dart";
import "dart:math" as math;

void main() {
  List<Coordinate> COORDS_1 = [Coordinate(1, 1), Coordinate(2, 2), Coordinate(3, 3)];
  List<Coordinate> COORDS_EMPTY = [];

  group("CoordinateArraysTest - ", () {
    test("testPtNotInList", () {
      expect(
          CoordinateArrays.ptNotInList([Coordinate(1, 1), Coordinate(2, 2), Coordinate(3, 3)], [Coordinate(1, 1), Coordinate(1, 2), Coordinate(1, 3)])
              .equals2D(Coordinate(2, 2)),
          true);

      expect(
          CoordinateArrays.ptNotInList([Coordinate(1, 1), Coordinate(2, 2), Coordinate(3, 3)], [Coordinate(1, 1), Coordinate(2, 2), Coordinate(3, 3)]) == null,
          true);
    });

    test("testEnvelopes", () {
      expect(CoordinateArrays.envelope(COORDS_1), Envelope(1, 3, 1, 3));
      expect(CoordinateArrays.envelope(COORDS_EMPTY), Envelope.empty());
      expect(CoordinateArrays.equals(CoordinateArrays.intersection(COORDS_1, Envelope(1, 2, 1, 2)), [Coordinate(1, 1), Coordinate(2, 2)]), true);

      expect(CoordinateArrays.equals(CoordinateArrays.intersection(COORDS_1, Envelope(10, 20, 10, 20)), COORDS_EMPTY), true);

      expect(CoordinateArrays.equals(CoordinateArrays.intersection(COORDS_EMPTY, Envelope(1, 2, 1, 2)), COORDS_EMPTY), true);

      expect(CoordinateArrays.equals(CoordinateArrays.intersection(COORDS_1, Envelope.empty()), COORDS_EMPTY), true);
    });
  });

  group("CoordinateListTest - ", () {
    test("testForward and reverse", () {
      checkValue(coordList([0, 0, 1, 1, 2, 2]), [0, 0, 1, 1, 2, 2]);
      checkValue(coordList([0, 0, 1, 1, 2, 2]).reversed.toList(), [2, 2, 1, 1, 0, 0]);
    });
  });

  group("CoordinateSequencesTest - ", () {
    test("testCopyToLargerDim", () {
      PackedCoordinateSequenceFactory csFactory = PackedCoordinateSequenceFactory();
      CoordinateSequence cs2D = createTestSequence(csFactory, 10, 2);
      CoordinateSequence cs3D = csFactory.createSizeDim(10, 3);
      CoordinateSequences.copy(cs2D, 0, cs3D, 0, cs3D.size());
      expect(CoordinateSequences.isEqual(cs2D, cs3D), true);
    });
    test("testCopyToSmallerDim", () {
      PackedCoordinateSequenceFactory csFactory = PackedCoordinateSequenceFactory();
      CoordinateSequence cs3D = createTestSequence(csFactory, 10, 3);
      CoordinateSequence cs2D = csFactory.createSizeDim(10, 2);
      CoordinateSequences.copy(cs3D, 0, cs2D, 0, cs2D.size());
      expect(CoordinateSequences.isEqual(cs2D, cs3D), true);
    });
    test("testScrollRing", () {
      doTestScrollRing(CoordinateArraySequenceFactory(), 2);
      doTestScrollRing(CoordinateArraySequenceFactory(), 3);
      doTestScrollRing(PackedCoordinateSequenceFactory.DOUBLE_FACTORY, 2);
      doTestScrollRing(PackedCoordinateSequenceFactory.DOUBLE_FACTORY, 4);
    });
    test("testScroll", () {
      doTestScroll(CoordinateArraySequenceFactory(), 2);
      doTestScroll(CoordinateArraySequenceFactory(), 3);
      doTestScroll(PackedCoordinateSequenceFactory.DOUBLE_FACTORY, 2);
      doTestScroll(PackedCoordinateSequenceFactory.DOUBLE_FACTORY, 4);
    });
    test("testIndexOf", () {
      doTestIndexOf(CoordinateArraySequenceFactory(), 2);
      doTestIndexOf(PackedCoordinateSequenceFactory.DOUBLE_FACTORY, 5);
    });
    test("testMinCoordinateIndex", () {
      doTestMinCoordinateIndex(CoordinateArraySequenceFactory(), 2);
      doTestMinCoordinateIndex(PackedCoordinateSequenceFactory.DOUBLE_FACTORY, 5);
    });
    test("testIsRing", () {
      doTestIsRing(CoordinateArraySequenceFactory(), 2);
      doTestIsRing(PackedCoordinateSequenceFactory.DOUBLE_FACTORY, 5);
    });
    test("testCopy", () {
      doTestCopy(CoordinateArraySequenceFactory(), 2);
      doTestCopy(PackedCoordinateSequenceFactory.DOUBLE_FACTORY, 5);
    });
    test("testReverse", () {
      doTestReverse(CoordinateArraySequenceFactory(), 2);
      doTestReverse(PackedCoordinateSequenceFactory.DOUBLE_FACTORY, 5);
    });
  });

//  group("CoordinateArraysTest - ", () {
//    test("", () {
//    });
//  });
}

void checkValue(List<Coordinate> coordArray, List<double> ords) {
  expect(coordArray.length * 2, ords.length);

  for (int i = 0; i < coordArray.length; i += 2) {
    Coordinate pt = coordArray[i];
    expect(pt.x, ords[2 * i]);
    expect(pt.y, ords[2 * i + 1]);
  }
}

List<Coordinate> coordList(List<double> ords) {
  List<Coordinate> cl = [];
  for (int i = 0; i < ords.length; i += 2) {
    cl.add(Coordinate(ords[i], ords[i + 1]));
//    CollectionsUtils.addIfNotEqualToLast(cl, Coordinate(ords[i], ords[i + 1]));
  }
  return cl;
}

final List<List<double>> ordinateValues = [
  [75.76, 77.43],
  [41.35, 90.75],
  [73.74, 41.67],
  [20.87, 86.49],
  [17.49, 93.59],
  [67.75, 80.63],
  [63.01, 52.57],
  [32.9, 44.44],
  [79.36, 29.8],
  [38.17, 88.0],
  [19.31, 49.71],
  [57.03, 19.28],
  [63.76, 77.35],
  [45.26, 85.15],
  [51.71, 50.38],
  [92.16, 19.85],
  [64.18, 27.7],
  [64.74, 65.1],
  [80.07, 13.55],
  [55.54, 94.07]
];

CoordinateSequence createSequenceFromOrdinates(CoordinateSequenceFactory csFactory, int dim) {
  CoordinateSequence sequence = csFactory.createSizeDim(ordinateValues.length, dim);
  for (int i = 0; i < ordinateValues.length; i++) {
    sequence.setOrdinate(i, 0, ordinateValues[i][0]);
    sequence.setOrdinate(i, 1, ordinateValues[i][1]);
  }
  return fillNonPlanarDimensions(sequence);
}

CoordinateSequence createTestSequence(CoordinateSequenceFactory csFactory, int size, int dim) {
  CoordinateSequence cs = csFactory.createSizeDim(size, dim);
  // initialize with a data signature where coords look like [1, 10, 100, ...]
  for (int i = 0; i < size; i++) {
    for (int d = 0; d < dim; d++) {
      cs.setOrdinate(i, d, (i * math.pow(10, d)).toDouble());
    }
  }
  return cs;
}

/// @deprecated only use to update in conjunction with {@link this.ttestCreateRandomOrdinates}
CoordinateSequence createRandomTestSequence(CoordinateSequenceFactory csFactory, int size, int dim, math.Random rnd, Envelope range, PrecisionModel pm) {
  CoordinateSequence cs = csFactory.createSizeDim(size, dim);
  for (int i = 0; i < size; i++) {
    cs.setOrdinate(i, 0, pm.makePrecise(range.getWidth() * rnd.nextDouble() + range.getMinX()));
    cs.setOrdinate(i, 1, pm.makePrecise(range.getHeight() * rnd.nextDouble() + range.getMinY()));
  }

  return fillNonPlanarDimensions(cs);
}

void doTestReverse(CoordinateSequenceFactory factory, int dimension) {
  // arrange
  CoordinateSequence sequence = createSequenceFromOrdinates(factory, dimension);
  CoordinateSequence reversed = sequence.copy();

  // act
  CoordinateSequences.reverse(reversed);

  // assert
  for (int i = 0; i < sequence.size(); i++) checkCoordinateAt(sequence, i, reversed, sequence.size() - i - 1, dimension);
}

void doTestCopy(CoordinateSequenceFactory factory, int dimension) {
  // arrange
  CoordinateSequence sequence = createSequenceFromOrdinates(factory, dimension);
  if (sequence.size() <= 7) {
    print("sequence has a size of ${sequence.size()}. Execution of this test needs a sequence " + "with more than 6 coordinates.");
    return;
  }

  CoordinateSequence fullCopy = factory.createSizeDim(sequence.size(), dimension);
  CoordinateSequence partialCopy = factory.createSizeDim(sequence.size() - 5, dimension);

  // act
  CoordinateSequences.copy(sequence, 0, fullCopy, 0, sequence.size());
  CoordinateSequences.copy(sequence, 2, partialCopy, 0, partialCopy.size());

  // assert
  for (int i = 0; i < fullCopy.size(); i++) checkCoordinateAt(sequence, i, fullCopy, i, dimension);
  for (int i = 0; i < partialCopy.size(); i++) checkCoordinateAt(sequence, 2 + i, partialCopy, i, dimension);

  // ToDo test if dimensions don't match
}

void doTestIsRing(CoordinateSequenceFactory factory, int dimension) {
  // arrange
  CoordinateSequence ring = createCircle(factory, dimension, Coordinate.emptyDefault(), 5);
  CoordinateSequence noRing = createCircularString(factory, dimension, Coordinate.emptyDefault(), 5, 0.1, 22);
  CoordinateSequence empty = createAlmostRing(factory, dimension, 0);
  CoordinateSequence incomplete1 = createAlmostRing(factory, dimension, 1);
  CoordinateSequence incomplete2 = createAlmostRing(factory, dimension, 2);
  CoordinateSequence incomplete3 = createAlmostRing(factory, dimension, 3);
  CoordinateSequence incomplete4a = createAlmostRing(factory, dimension, 4);
  CoordinateSequence incomplete4b = CoordinateSequences.ensureValidRing(factory, incomplete4a);

  // act
  bool isRingRing = CoordinateSequences.isRing(ring);
  bool isRingNoRing = CoordinateSequences.isRing(noRing);
  bool isRingEmpty = CoordinateSequences.isRing(empty);
  bool isRingIncomplete1 = CoordinateSequences.isRing(incomplete1);
  bool isRingIncomplete2 = CoordinateSequences.isRing(incomplete2);
  bool isRingIncomplete3 = CoordinateSequences.isRing(incomplete3);
  bool isRingIncomplete4a = CoordinateSequences.isRing(incomplete4a);
  bool isRingIncomplete4b = CoordinateSequences.isRing(incomplete4b);

  // assert
  expect(isRingRing, true);
  expect(!isRingNoRing, true);
  expect(isRingEmpty, true);
  expect(!isRingIncomplete1, true);
  expect(!isRingIncomplete2, true);
  expect(!isRingIncomplete3, true);
  expect(!isRingIncomplete4a, true);
  expect(isRingIncomplete4b, true);
}

void doTestIndexOf(CoordinateSequenceFactory factory, int dimension) {
  // arrange
  CoordinateSequence sequence = createSequenceFromOrdinates(factory, dimension);

  // act & assert
  List<Coordinate> coordinates = sequence.toCoordinateArray();
  for (int i = 0; i < sequence.size(); i++) expect(i, CoordinateSequences.indexOf(coordinates[i], sequence));
}

void doTestMinCoordinateIndex(CoordinateSequenceFactory factory, int dimension) {
  CoordinateSequence sequence = createSequenceFromOrdinates(factory, dimension);
  if (sequence.size() <= 6) {
    print("sequence has a size of ${sequence.size()}. Execution of this test needs a sequence " + "with more than 5 coordinates.");
    return;
  }

  int minIndex = (sequence.size() / 2).round();
  sequence.setOrdinate(minIndex, 0, 5);
  sequence.setOrdinate(minIndex, 1, 5);

  expect(minIndex, CoordinateSequences.minCoordinateIndex(sequence));
  expect(minIndex, CoordinateSequences.minCoordinateIndexWithRange(sequence, 2, sequence.size() - 2));
}

void doTestScroll(CoordinateSequenceFactory factory, int dimension) {
  // arrange
  CoordinateSequence sequence = createCircularString(factory, dimension, Coordinate(20, 20), 7.0, 0.1, 22);
  CoordinateSequence scrolled = sequence.copy();

  // act
  CoordinateSequences.scrollWithIndex(scrolled, 12);

  // assert
  int io = 12;
  for (int iis = 0; iis < scrolled.size() - 1; iis++) {
    checkCoordinateAt(sequence, io, scrolled, iis, dimension);
    io++;
    io %= scrolled.size();
  }
}

void doTestScrollRing(CoordinateSequenceFactory factory, int dimension) {
  // arrange
  //System.out.println("Testing '" + factory.getClass().getSimpleName() + "' with dim=" +dimension );
  CoordinateSequence sequence = createCircle(factory, dimension, Coordinate(10, 10), 9.0);
  CoordinateSequence scrolled = sequence.copy();

  // act
  CoordinateSequences.scrollWithIndex(scrolled, 12);

  // assert
  int io = 12;
  for (int iis = 0; iis < scrolled.size() - 1; iis++) {
    checkCoordinateAt(sequence, io, scrolled, iis, dimension);
    io++;
    io %= scrolled.size() - 1;
  }
  checkCoordinateAt(scrolled, 0, scrolled, scrolled.size() - 1, dimension);
}

void checkCoordinateAt(CoordinateSequence seq1, int pos1, CoordinateSequence seq2, int pos2, int dim) {
  expect(seq1.getOrdinate(pos1, 0), seq2.getOrdinate(pos2, 0));
  expect(seq1.getOrdinate(pos1, 1), seq2.getOrdinate(pos2, 1));

  // check additional ordinates
  for (int j = 2; j < dim; j++) {
    expect(seq1.getOrdinate(pos1, j), seq2.getOrdinate(pos2, j));
  }
}

CoordinateSequence createAlmostRing(CoordinateSequenceFactory factory, int dimension, int num) {
  if (num > 4) num = 4;

  CoordinateSequence sequence = factory.createSizeDim(num, dimension);
  if (num == 0) return fillNonPlanarDimensions(sequence);

  sequence.setOrdinate(0, 0, 10);
  sequence.setOrdinate(0, 0, 10);
  if (num == 1) return fillNonPlanarDimensions(sequence);

  sequence.setOrdinate(0, 0, 20);
  sequence.setOrdinate(0, 0, 10);
  if (num == 2) return fillNonPlanarDimensions(sequence);

  sequence.setOrdinate(0, 0, 20);
  sequence.setOrdinate(0, 0, 20);
  if (num == 3) return fillNonPlanarDimensions(sequence);

  sequence.setOrdinate(0, 0, 10.0000000000001);
  sequence.setOrdinate(0, 0, 9.9999999999999);
  return fillNonPlanarDimensions(sequence);
}

CoordinateSequence fillNonPlanarDimensions(CoordinateSequence seq) {
  if (seq.getDimension() < 3) return seq;

  for (int i = 0; i < seq.size(); i++) {
    for (int j = 2; j < seq.getDimension(); j++) {
      seq.setOrdinate(i, j, (i * math.pow(10, j - 1)).toDouble());
    }
  }

  return seq;
}

CoordinateSequence createCircle(CoordinateSequenceFactory factory, int dimension, Coordinate center, double radius) {
  // Get a complete circular string
  CoordinateSequence res = createCircularString(factory, dimension, center, radius, 0, 49);

  // ensure it is closed
  for (int i = 0; i < dimension; i++) res.setOrdinate(48, i, res.getOrdinate(0, i));

  return res;
}

CoordinateSequence createCircularString(CoordinateSequenceFactory factory, int dimension, Coordinate center, double radius, double startAngle, int numPoints) {
  final int numSegmentsCircle = 48;
  final double angleCircle = 2 * math.pi;
  final double angleStep = angleCircle / numSegmentsCircle;

  CoordinateSequence sequence = factory.createSizeDim(numPoints, dimension);
  PrecisionModel pm = PrecisionModel.fixedPrecision(100);
  double angle = startAngle;
  for (int i = 0; i < numPoints; i++) {
    double dx = math.cos(angle) * radius;
    sequence.setOrdinate(i, 0, pm.makePrecise(center.x + dx));
    double dy = math.sin(angle) * radius;
    sequence.setOrdinate(i, 1, pm.makePrecise(center.y + dy));

    // set other ordinate values to predictable values
    for (int j = 2; j < dimension; j++) {
      sequence.setOrdinate(i, j, (math.pow(10, j - 1) * i).toDouble());
    }

    angle += angleStep;
    angle %= angleCircle;
  }

  return sequence;
}
