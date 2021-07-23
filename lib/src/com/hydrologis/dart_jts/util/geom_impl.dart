part of dart_jts;

/// Utility functions for manipulating {@link CoordinateSequence}s
///
/// @version 1.7
class CoordinateSequences {
  /// Reverses the coordinates in a sequence in-place.
  static void reverse(CoordinateSequence seq) {
    int last = seq.size() - 1;
    int mid = last ~/ 2;
    for (int i = 0; i <= mid; i++) {
      swap(seq, i, last - i);
    }
  }

  /// Swaps two coordinates in a sequence.
  ///
  /// @param seq the sequence to modify
  /// @param i the index of a coordinate to swap
  /// @param j the index of a coordinate to swap
  static void swap(CoordinateSequence seq, int i, int j) {
    if (i == j) return;
    for (int dim = 0; dim < seq.getDimension(); dim++) {
      double tmp = seq.getOrdinate(i, dim);
      seq.setOrdinate(i, dim, seq.getOrdinate(j, dim));
      seq.setOrdinate(j, dim, tmp);
    }
  }

  /// Copies a section of a {@link CoordinateSequence} to another {@link CoordinateSequence}.
  /// The sequences may have different dimensions;
  /// in this case only the common dimensions are copied.
  ///
  /// @param src the sequence to copy from
  /// @param srcPos the position in the source sequence to start copying at
  /// @param dest the sequence to copy to
  /// @param destPos the position in the destination sequence to copy to
  /// @param length the number of coordinates to copy
  static void copy(CoordinateSequence src, int srcPos, CoordinateSequence dest,
      int destPos, int length) {
    for (int i = 0; i < length; i++) {
      copyCoord(src, srcPos + i, dest, destPos + i);
    }
  }

  /// Copies a coordinate of a {@link CoordinateSequence} to another {@link CoordinateSequence}.
  /// The sequences may have different dimensions;
  /// in this case only the common dimensions are copied.
  ///
  /// @param src the sequence to copy from
  /// @param srcPos the source coordinate to copy
  /// @param dest the sequence to copy to
  /// @param destPos the destination coordinate to copy to
  static void copyCoord(CoordinateSequence src, int srcPos,
      CoordinateSequence dest, int destPos) {
    int minDim = math.min(src.getDimension(), dest.getDimension());
    for (int dim = 0; dim < minDim; dim++) {
      dest.setOrdinate(destPos, dim, src.getOrdinate(srcPos, dim));
    }
  }

  /// Tests whether a {@link CoordinateSequence} forms a valid {@link LinearRing},
  /// by checking the sequence length and closure
  /// (whether the first and last points are identical in 2D).
  /// Self-intersection is not checked.
  ///
  /// @param seq the sequence to test
  /// @return true if the sequence is a ring
  /// @see LinearRing
  static bool isRing(CoordinateSequence seq) {
    int n = seq.size();
    if (n == 0) return true;
    // too few points
    if (n <= 3) return false;
    // test if closed
    return seq.getOrdinate(0, CoordinateSequence.X) ==
            seq.getOrdinate(n - 1, CoordinateSequence.X) &&
        seq.getOrdinate(0, CoordinateSequence.Y) ==
            seq.getOrdinate(n - 1, CoordinateSequence.Y);
  }

  /// Ensures that a CoordinateSequence forms a valid ring,
  /// returning a new closed sequence of the correct length if required.
  /// If the input sequence is already a valid ring, it is returned
  /// without modification.
  /// If the input sequence is too short or is not closed,
  /// it is extended with one or more copies of the start point.
  ///
  /// @param fact the CoordinateSequenceFactory to use to create the new sequence
  /// @param seq the sequence to test
  /// @return the original sequence, if it was a valid ring, or a new sequence which is valid.
  static CoordinateSequence ensureValidRing(
      CoordinateSequenceFactory fact, CoordinateSequence seq) {
    int n = seq.size();
    // empty sequence is valid
    if (n == 0) return seq;
    // too short - make a new one
    if (n <= 3) return _createClosedRing(fact, seq, 4);

    bool isClosed = seq.getOrdinate(0, CoordinateSequence.X) ==
            seq.getOrdinate(n - 1, CoordinateSequence.X) &&
        seq.getOrdinate(0, CoordinateSequence.Y) ==
            seq.getOrdinate(n - 1, CoordinateSequence.Y);
    if (isClosed) return seq;
    // make a new closed ring
    return _createClosedRing(fact, seq, n + 1);
  }

  static CoordinateSequence _createClosedRing(
      CoordinateSequenceFactory fact, CoordinateSequence seq, int size) {
    CoordinateSequence newseq = fact.createSizeDim(size, seq.getDimension());
    int n = seq.size();
    copy(seq, 0, newseq, 0, n);
    // fill remaining coordinates with start point
    for (int i = n; i < size; i++) copy(seq, 0, newseq, i, 1);
    return newseq;
  }

  static CoordinateSequence extend(
      CoordinateSequenceFactory fact, CoordinateSequence seq, int size) {
    CoordinateSequence newseq = fact.createSizeDim(size, seq.getDimension());
    int n = seq.size();
    copy(seq, 0, newseq, 0, n);
    // fill remaining coordinates with end point, if it exists
    if (n > 0) {
      for (int i = n; i < size; i++) {
        copy(seq, n - 1, newseq, i, 1);
      }
    }
    return newseq;
  }

  /// Tests whether two {@link CoordinateSequence}s are equal.
  /// To be equal, the sequences must be the same length.
  /// They do not need to be of the same dimension,
  /// but the ordinate values for the smallest dimension of the two
  /// must be equal.
  /// Two <code>NaN</code> ordinates values are considered to be equal.
  ///
  /// @param cs1 a CoordinateSequence
  /// @param cs2 a CoordinateSequence
  /// @return true if the sequences are equal in the common dimensions
  static bool isEqual(CoordinateSequence cs1, CoordinateSequence cs2) {
    int cs1Size = cs1.size();
    int cs2Size = cs2.size();
    if (cs1Size != cs2Size) return false;
    int dim = math.min(cs1.getDimension(), cs2.getDimension());
    for (int i = 0; i < cs1Size; i++) {
      for (int d = 0; d < dim; d++) {
        double v1 = cs1.getOrdinate(i, d);
        double v2 = cs2.getOrdinate(i, d);
        if (cs1.getOrdinate(i, d) == cs2.getOrdinate(i, d)) continue;
        // special check for NaNs
        if (v1.isNaN && v2.isNaN) continue;
        return false;
      }
    }
    return true;
  }

  /// Creates a string representation of a {@link CoordinateSequence}.
  /// The format is:
  /// <pre>
  ///   ( ord0,ord1.. ord0,ord1,...  ... )
  /// </pre>
  ///
  /// @param cs the sequence to output
  /// @return the string representation of the sequence
  static String asString(CoordinateSequence cs) {
    int size = cs.size();
    if (size == 0) return "()";
    int dim = cs.getDimension();
    StringBuffer builder = new StringBuffer();
    builder.write('(');
    for (int i = 0; i < size; i++) {
      if (i > 0) builder.write(" ");
      for (int d = 0; d < dim; d++) {
        if (d > 0) builder.write(",");
        builder.write(cs.getOrdinate(i, d).toStringAsFixed(1));
      }
    }
    builder.write(')');
    return builder.toString();
  }

  ///  Returns the minimum coordinate, using the usual lexicographic comparison.
  ///
  ///@param  seq  the coordinate sequence to search
  ///@return  the minimum coordinate in the sequence, found using <code>compareTo</code>
  ///@see Coordinate#compareTo(Object)
  static Coordinate minCoordinate(CoordinateSequence seq) {
    Coordinate? minCoord;
    for (int i = 0; i < seq.size(); i++) {
      Coordinate testCoord = seq.getCoordinate(i);
      if (minCoord == null || minCoord.compareTo(testCoord) > 0) {
        minCoord = testCoord;
      }
    }
    return minCoord!;
  }

  ///  Returns the index of the minimum coordinate of the whole
  ///  coordinate sequence, using the usual lexicographic comparison.
  ///
  ///@param  seq  the coordinate sequence to search
  ///@return  the index of the minimum coordinate in the sequence, found using <code>compareTo</code>
  ///@see Coordinate#compareTo(Object)
  static int minCoordinateIndex(CoordinateSequence seq) {
    return minCoordinateIndexWithRange(seq, 0, seq.size() - 1);
  }

  ///  Returns the index of the minimum coordinate of a part of
  ///  the coordinate sequence (defined by {@code from} and {@code to},
  ///  using the usual lexicographic comparison.
  ///
  ///@param  seq   the coordinate sequence to search
  ///@param  from  the lower search index
  ///@param  to    the upper search index
  ///@return  the index of the minimum coordinate in the sequence, found using <code>compareTo</code>
  ///@see Coordinate#compareTo(Object)
  static int minCoordinateIndexWithRange(
      CoordinateSequence seq, int from, int to) {
    int minCoordIndex = -1;
    Coordinate? minCoord;
    for (int i = from; i <= to; i++) {
      Coordinate testCoord = seq.getCoordinate(i);
      if (minCoord == null || minCoord.compareTo(testCoord) > 0) {
        minCoord = testCoord;
        minCoordIndex = i;
      }
    }
    return minCoordIndex;
  }

  ///  Shifts the positions of the coordinates until <code>firstCoordinate</code>
  ///  is first.
  ///
  ///@param  seq      the coordinate sequence to rearrange
  ///@param  firstCoordinate  the coordinate to make first
  static void scroll(CoordinateSequence seq, Coordinate firstCoordinate) {
    int i = indexOf(firstCoordinate, seq);
    if (i <= 0) return;
    scrollWithIndex(seq, i);
  }

  ///  Shifts the positions of the coordinates until the coordinate at  <code>firstCoordinateIndex</code>
  ///  is first.
  ///
  ///@param  seq      the coordinate sequence to rearrange
  ///@param  indexOfFirstCoordinate  the index of the coordinate to make first
  static void scrollWithIndex(
      CoordinateSequence seq, int indexOfFirstCoordinate) {
    scrollWithIndexAndRingcheck(
        seq, indexOfFirstCoordinate, CoordinateSequences.isRing(seq));
  }

  ///  Shifts the positions of the coordinates until the coordinate at  <code>firstCoordinateIndex</code>
  ///  is first.
  ///
  ///@param  seq      the coordinate sequence to rearrange
  ///@param  indexOfFirstCoordinate
  ///                 the index of the coordinate to make first
  ///@param  ensureRing
  ///                 makes sure that {@code} will be a closed ring upon exit
  static void scrollWithIndexAndRingcheck(
      CoordinateSequence seq, int indexOfFirstCoordinate, bool ensureRing) {
    int i = indexOfFirstCoordinate;
    if (i <= 0) return;

    // make a copy of the sequence
    CoordinateSequence copy = seq.copy();

    // test if ring, determine last index
    int last = ensureRing ? seq.size() - 1 : seq.size();

    // fill in values
    for (int j = 0; j < last; j++) {
      for (int k = 0; k < seq.getDimension(); k++)
        seq.setOrdinate(
            j, k, copy.getOrdinate((indexOfFirstCoordinate + j) % last, k));
    }

    // Fix the ring (first == last)
    if (ensureRing) {
      for (int k = 0; k < seq.getDimension(); k++)
        seq.setOrdinate(last, k, seq.getOrdinate(0, k));
    }
  }

  ///  Returns the index of <code>coordinate</code> in a {@link CoordinateSequence}
  ///  The first position is 0; the second, 1; etc.
  ///
  ///@param  coordinate   the <code>Coordinate</code> to search for
  ///@param  seq  the coordinate sequence to search
  ///@return              the position of <code>coordinate</code>, or -1 if it is
  ///      not found
  static int indexOf(Coordinate coordinate, CoordinateSequence seq) {
    for (int i = 0; i < seq.size(); i++) {
      if (coordinate.x == seq.getOrdinate(i, CoordinateSequence.X) &&
          coordinate.y == seq.getOrdinate(i, CoordinateSequence.Y)) {
        return i;
      }
    }
    return -1;
  }
}

/// Useful utility functions for handling Coordinate arrays
///
/// @version 1.7
class CoordinateArrays {
  static final List<Coordinate> coordArrayType = [];

  /// A {@link Comparator} for {@link Coordinate} arrays
  /// in the forward direction of their coordinates,
  /// using lexicographic ordering.
  static Comparator<List<Coordinate>> forwardComparator =
      (o1, o2) => CoordinateArrays.compare(o1, o2);

  /// A {@link Comparator} for {@link Coordinate} arrays
  /// modulo their directionality.
  /// E.g. if two coordinate arrays are identical but reversed
  /// they will compare as equal under this ordering.
  /// If the arrays are not equal, the ordering returned
  /// is the ordering in the forward direction.
  ///
  static Comparator<List<Coordinate>> bidirectionalComparator = (pts1, pts2) {
    if (pts1.length < pts2.length) return -1;
    if (pts1.length > pts2.length) return 1;

    if (pts1.isEmpty) return 0;

    int forwardComp = CoordinateArrays.compare(pts1, pts2);
    bool isEqualRev = isEqualReversed(pts1, pts2);
    if (isEqualRev) {
      return 0;
    }
    return forwardComp;
  };

  CoordinateArrays();

  /// Determine dimension based on subclass of {@link Coordinate}.
  ///
  /// @param pts supplied coordinates
  /// @return number of ordinates recorded
  static int dimension(List<Coordinate>? pts) {
    if (pts == null || pts.length == 0) {
      return 3; // unknown, assume default
    }
    int dimension = 0;
    pts.forEach((coordinate) {
      dimension = math.max(dimension, Coordinates.dimension(coordinate));
    });
    return dimension;
  }

  /// Determine number of measures based on subclass of {@link Coordinate}.
  ///
  /// @param pts supplied coordinates
  /// @return number of measures recorded
  static int measures(List<Coordinate>? pts) {
    if (pts == null || pts.isEmpty) {
      return 0; // unknown, assume default
    }
    int measures = 0;
    pts.forEach((coordinate) {
      measures = math.max(measures, Coordinates.measures(coordinate));
    });
    return measures;
  }

  /// Tests whether an array of {@link Coordinate}s forms a ring,
  /// by checking length and closure.
  /// Self-intersection is not checked.
  ///
  /// @param pts an array of Coordinates
  /// @return true if the coordinate form a ring.
  static isRing(List<Coordinate> pts) {
    if (pts.length < 4) return false;
    if (!pts[0].equals2D(pts[pts.length - 1])) return false;
    return true;
  }

  /// Finds a point in a list of points which is not contained in another list of points
  /// @param testPts the {@link Coordinate}s to test
  /// @param pts an array of {@link Coordinate}s to test the input points against
  /// @return a {@link Coordinate} from <code>testPts</code> which is not in <code>pts</code>, '
  /// or <code>null</code>
  static Coordinate? ptNotInList(
      List<Coordinate> testPts, List<Coordinate> pts) {
    for (int i = 0; i < testPts.length; i++) {
      Coordinate testPt = testPts[i];
      if (CoordinateArrays.indexOf(testPt, pts) < 0) return testPt;
    }
    return null;
  }

  /// Compares two {@link Coordinate} arrays
  /// in the forward direction of their coordinates,
  /// using lexicographic ordering.
  ///
  /// @param pts1
  /// @param pts2
  /// @return an integer indicating the order
  static int compare(List<Coordinate> pts1, List<Coordinate> pts2) {
    int i = 0;
    while (i < pts1.length && i < pts2.length) {
      int compare = pts1[i].compareTo(pts2[i]);
      if (compare != 0) return compare;
      i++;
    }
    // handle situation when arrays are of different length
    if (i < pts2.length) return -1;
    if (i < pts1.length) return 1;

    return 0;
  }

  /// Determines which orientation of the {@link Coordinate} array
  /// is (overall) increasing.
  /// In other words, determines which end of the array is "smaller"
  /// (using the standard ordering on {@link Coordinate}).
  /// Returns an integer indicating the increasing direction.
  /// If the sequence is a palindrome, it is defined to be
  /// oriented in a positive direction.
  ///
  /// @param pts the array of Coordinates to test
  /// @return <code>1</code> if the array is smaller at the start
  /// or is a palindrome,
  /// <code>-1</code> if smaller at the end
  static int increasingDirection(List<Coordinate> pts) {
    for (int i = 0; i < pts.length / 2; i++) {
      int j = pts.length - 1 - i;
      // skip equal points on both ends
      int comp = pts[i].compareTo(pts[j]);
      if (comp != 0) return comp;
    }
    // array must be a palindrome - defined to be in positive direction
    return 1;
  }

  /// Determines whether two {@link Coordinate} arrays of equal length
  /// are equal in opposite directions.
  ///
  /// @param pts1
  /// @param pts2
  /// @return <code>true</code> if the two arrays are equal in opposite directions.
  static isEqualReversed(List<Coordinate> pts1, List<Coordinate> pts2) {
    for (int i = 0; i < pts1.length; i++) {
      Coordinate p1 = pts1[i];
      Coordinate p2 = pts2[pts1.length - i - 1];
      if (p1.compareTo(p2) != 0) return false;
    }
    return true;
  }

  /// Creates a deep copy of the argument {@link Coordinate} array.
  ///
  /// @param coordinates an array of Coordinates
  /// @return a deep copy of the input
  static List<Coordinate> copyDeep(List<Coordinate> coordinates) {
    List<Coordinate> copy = []; //..length = (coordinates.length);
    for (int i = 0; i < coordinates.length; i++) {
      copy.add(coordinates[i].copy());
      // copy[i] = coordinates[i].copy();
    }
    return copy;
  }

  /// Creates a deep copy of a given section of a source {@link Coordinate} array
  /// into a destination Coordinate array.
  /// The destination array must be an appropriate size to receive
  /// the copied coordinates.
  ///
  /// @param src an array of Coordinates
  /// @param srcStart the index to start copying from
  /// @param dest the
  /// @param destStart the destination index to start copying to
  /// @param length the number of items to copy
  static void copyDeepWithLength(List<Coordinate> src, int srcStart,
      List<Coordinate> dest, int destStart, int length) {
    for (int i = 0; i < length; i++) {
      dest[destStart + i] = src[srcStart + i].copy();
    }
  }

  /// Returns either the given coordinate array if its length is greater than the
  /// given amount, or an empty coordinate array.
  static List<Coordinate> atLeastNCoordinatesOrNothing(
      int n, List<Coordinate> c) {
    return c.length >= n ? c : <Coordinate>[];
  }

  /// If the coordinate array argument has repeated points,
  /// constructs a new array containing no repeated points.
  /// Otherwise, returns the argument.
  /// @see #hasRepeatedPoints(List<Coordinate>)
  static List<Coordinate> removeRepeatedPoints(List<Coordinate> coord) {
    if (!CollectionsUtils.hasRepeated(coord)) return coord;
    return CollectionsUtils.removeRepeated(coord);
  }

  /// Collapses a coordinate array to remove all null elements.
  ///
  /// @param coord the coordinate array to collapse
  /// @return an array containing only non-null elements
  static List<Coordinate> removeNull(List<Coordinate?> coord) {
    int nonNull = 0;
    for (int i = 0; i < coord.length; i++) {
      if (coord[i] != null) nonNull++;
    }
    List<Coordinate> newCoord = []; //..length = (nonNull);
    // empty case
    if (nonNull == 0) return newCoord;

    // int j = 0;
    for (int i = 0; i < coord.length; i++) {
      if (coord[i] != null) newCoord.add(coord[i]!);
      // if (coord[i] != null) newCoord[j++] = coord[i]!;
    }
    return newCoord;
  }

  /// Reverses the coordinates in an array in-place.
  static void reverse(List<Coordinate> coord) {
    int last = coord.length - 1;
    int mid = (last / 2).floor();
    for (int i = 0; i <= mid; i++) {
      Coordinate tmp = coord[i];
      coord[i] = coord[last - i];
      coord[last - i] = tmp;
    }
  }

  /// Returns true if the two arrays are identical, both null, or pointwise
  /// equal (as compared using Coordinate#equals)
  /// @see Coordinate#equals(Object)
  static bool equals(List<Coordinate>? coord1, List<Coordinate>? coord2) {
    if (coord1 == coord2) return true;
    if (coord1 == null || coord2 == null) return false;
    if (coord1.length != coord2.length) return false;
    for (int i = 0; i < coord1.length; i++) {
      if (coord1[i] != coord2[i]) return false;
    }
    return true;
  }

  /// Returns true if the two arrays are identical, both null, or pointwise
  /// equal, using a user-defined {@link Comparator} for {@link Coordinate} s
  ///
  /// @param coord1 an array of Coordinates
  /// @param coord2 an array of Coordinates
  /// @param coordinateComparator a Comparator for Coordinates
  static bool equalsWithComparator(List<Coordinate>? coord1,
      List<Coordinate>? coord2, Comparator coordinateComparator) {
    if (coord1 == coord2) return true;
    if (coord1 == null || coord2 == null) return false;
    if (coord1.length != coord2.length) return false;
    for (int i = 0; i < coord1.length; i++) {
      if (coordinateComparator(coord1[i], coord2[i]) != 0) return false;
    }
    return true;
  }

  ///  Returns the minimum coordinate, using the usual lexicographic comparison.
  ///
  ///@param  coordinates  the array to search
  ///@return              the minimum coordinate in the array, found using <code>compareTo</code>
  ///@see Coordinate#compareTo(Object)
  static Coordinate minCoordinate(List<Coordinate> coordinates) {
    Coordinate? minCoord;
    for (int i = 0; i < coordinates.length; i++) {
      if (minCoord == null || minCoord.compareTo(coordinates[i]) > 0) {
        minCoord = coordinates[i];
      }
    }
    return minCoord!;
  }

  ///  Shifts the positions of the coordinates until <code>firstCoordinate</code>
  ///  is first.
  ///
  ///@param  coordinates      the array to rearrange
  ///@param  firstCoordinate  the coordinate to make first
  static void scroll(List<Coordinate> coordinates, Coordinate firstCoordinate) {
    int i = coordinates.indexOf(firstCoordinate);
    if (i < 0) return;

    var newCoordinates = CollectionsUtils.shiftToFirst(coordinates, i)!;
    coordinates.clear();
    coordinates.addAll(newCoordinates);
  }

  ///  Returns the index of <code>coordinate</code> in <code>coordinates</code>.
  ///  The first position is 0; the second, 1; etc.
  ///
  ///@param  coordinate   the <code>Coordinate</code> to search for
  ///@param  coordinates  the array to search
  ///@return              the position of <code>coordinate</code>, or -1 if it is
  ///      not found
  static int indexOf(Coordinate coordinate, List<Coordinate> coordinates) {
    for (int i = 0; i < coordinates.length; i++) {
      if (coordinate == coordinates[i]) {
        return i;
      }
    }
    return -1;
  }

  /// Extracts a subsequence of the input {@link Coordinate} array
  /// from indices <code>start</code> to
  /// <code>end</code> (inclusive).
  /// The input indices are clamped to the array size;
  /// If the end index is less than the start index,
  /// the extracted array will be empty.
  ///
  /// @param pts the input array
  /// @param start the index of the start of the subsequence to extract
  /// @param end the index of the end of the subsequence to extract
  /// @return a subsequence of the input array
  static List<Coordinate> extract(List<Coordinate> pts, int start, int end) {
    start = MathUtils.clamp(start, 0, pts.length).toInt();
    end = MathUtils.clamp(end, -1, pts.length).toInt();

    int npts = end - start + 1;
    if (end < 0) npts = 0;
    if (start >= pts.length) npts = 0;
    if (end < start) npts = 0;

    List<Coordinate> extractPts = []; //..length = (npts);
    if (npts == 0) return extractPts;

    // int iPts = 0;
    for (int i = start; i <= end; i++) {
      extractPts.add(pts[i]);
      // extractPts[iPts++] = pts[i];
    }
    return extractPts;
  }

  /// Computes the envelope of the coordinates.
  ///
  /// @param coordinates the coordinates to scan
  /// @return the envelope of the coordinates
  static Envelope envelope(List<Coordinate> coordinates) {
    Envelope env = Envelope.empty();
    for (int i = 0; i < coordinates.length; i++) {
      env.expandToIncludeCoordinate(coordinates[i]);
    }
    return env;
  }

  /// Extracts the coordinates which intersect an {@link Envelope}.
  ///
  /// @param coordinates the coordinates to scan
  /// @param env the envelope to intersect with
  /// @return an array of the coordinates which intersect the envelope
  static List<Coordinate> intersection(
      List<Coordinate> coordinates, Envelope env) {
    List<Coordinate> coordList = [];
    for (int i = 0; i < coordinates.length; i++) {
      if (env.intersectsCoordinate(coordinates[i]))
        coordList.add(coordinates[i]);
    }
    return coordList;
  }
}

/**
 * Useful utility functions for handling Coordinate objects.
 */
class Coordinates {
  /**
   * Factory method providing access to common Coordinate implementations.
   *
   * @param dimension
   * @return created coordinate
   */
  static Coordinate create(int dimension) {
    return createWithMeasure(dimension, 0);
  }

  /**
   * Factory method providing access to common Coordinate implementations.
   *
   * @param dimension
   * @param measures
   * @return created coordinate
   */
  static Coordinate createWithMeasure(int dimension, int measures) {
    if (dimension == 2) {
      return new CoordinateXY();
    } else if (dimension == 3 && measures == 0) {
      return new Coordinate.empty2D();
    } else if (dimension == 3 && measures == 1) {
      return new CoordinateXYM.empty();
    } else if (dimension == 4 && measures == 1) {
      return new CoordinateXYZM.empty();
    }
    return new Coordinate.empty2D();
  }

  /**
   * Determine dimension based on subclass of {@link Coordinate}.
   *
   * @param coordinate supplied coordinate
   * @return number of ordinates recorded
   */
  static int dimension(Coordinate coordinate) {
    if (coordinate is CoordinateXY) {
      return 2;
    } else if (coordinate is CoordinateXYM) {
      return 3;
    } else if (coordinate is CoordinateXYZM) {
      return 4;
    } else if (coordinate is Coordinate) {
      return 3;
    }
    return 3;
  }

  /**
   * Determine number of measures based on subclass of {@link Coordinate}.
   *
   * @param coordinate supplied coordinate
   * @return number of measures recorded
   */
  static int measures(Coordinate coordinate) {
    if (coordinate is CoordinateXY) {
      return 0;
    } else if (coordinate is CoordinateXYM) {
      return 1;
    } else if (coordinate is CoordinateXYZM) {
      return 1;
    } else if (coordinate is Coordinate) {
      return 0;
    }
    return 0;
  }
}

/// Creates {@link CoordinateSequence}s represented as an array of {@link Coordinate}s.
///
/// @version 1.7
class CoordinateArraySequenceFactory implements CoordinateSequenceFactory {
  static final CoordinateArraySequenceFactory _instanceObject =
      CoordinateArraySequenceFactory._internal();

  factory CoordinateArraySequenceFactory() {
    return _instanceObject;
  }

  CoordinateArraySequenceFactory._internal();

  /// Returns a {@link CoordinateArraySequence} based on the given array (the array is
  /// not copied).
  ///
  /// @param coordinates
  ///            the coordinates, which may not be null nor contain null
  ///            elements
  CoordinateSequence create(List<Coordinate> coordinates) {
    return CoordinateArraySequence(coordinates);
  }

  /// @see org.locationtech.jts.geom.CoordinateSequenceFactory#create(org.locationtech.jts.geom.CoordinateSequence)
  CoordinateSequence createFromSequence(CoordinateSequence coordSeq) {
    return CoordinateArraySequence.fromSequence(coordSeq);
  }

  /// The created sequence dimension is clamped to be &lt;= 3.
  ///
  /// @see org.locationtech.jts.geom.CoordinateSequenceFactory#create(int, int)
  ///
  CoordinateSequence createSizeDim(int size, int dimension) {
    if (dimension > 3) dimension = 3;
//throw IllegalArgumentException("dimension must be <= 3");

// handle bogus dimension
    if (dimension < 2) dimension = 2;

    return CoordinateArraySequence.fromSizeDimension(size, dimension);
  }

  CoordinateSequence createSizeDimMeas(int size, int dimension, int measures) {
    int spatial = dimension - measures;

    if (measures > 1) {
      measures = 1; // clip measures
//throw IllegalArgumentException("measures must be <= 1");
    }
    if ((spatial) > 3) {
      spatial = 3; // clip spatial dimension
//throw IllegalArgumentException("spatial dimension must be <= 3");
    }

    if (spatial < 2) spatial = 2; // handle bogus spatial dimension

    return CoordinateArraySequence.fromSizeDimensionMeasures(
        size, spatial + measures, measures);
  }
}

/// A {@link CoordinateSequence} backed by an array of {@link Coordinate}s.
/// This is the implementation that {@link Geometry}s use by default.
/// Coordinates returned by #toArray and #getCoordinate are live --
/// modifications to them are actually changing the
/// CoordinateSequence's underlying data.
/// A dimension may be specified for the coordinates in the sequence,
/// which may be 2 or 3.
/// The actual coordinates will always have 3 ordinates,
/// but the dimension is useful as metadata in some situations.
///
/// @version 1.7
class CoordinateArraySequence extends CoordinateSequence {
  /// The actual dimension of the coordinates in the sequence.
  /// Allowable values are 2, 3 or 4.
  int dimension = 3;

  /// The number of measures of the coordinates in the sequence.
  /// Allowable values are 0 or 1.
  int measures = 0;

  List<Coordinate>? coordinates;

  /// Constructs a sequence based on the given array
  /// of {@link Coordinate}s (the
  /// array is not copied).
  /// The coordinate dimension defaults to 3.
  ///
  /// @param coordinates the coordinate array that will be referenced.
  CoordinateArraySequence(List<Coordinate> coordinates)
      : this.withDimensionMeasures(
            coordinates,
            CoordinateArrays.dimension(coordinates),
            CoordinateArrays.measures(coordinates));

  /// Constructs a sequence based on the given array
  /// of {@link Coordinate}s (the
  /// array is not copied).
  ///
  /// @param coordinates the coordinate array that will be referenced.
  /// @param dimension the dimension of the coordinates
  CoordinateArraySequence.withDimension(
      List<Coordinate> coordinates, int dimension)
      : this.withDimensionMeasures(
            coordinates, dimension, CoordinateArrays.measures(coordinates));

  /// Constructs a sequence based on the given array
  /// of {@link Coordinate}s (the
  /// array is not copied).
  ///
  /// @param coordinates the coordinate array that will be referenced.
  /// @param dimension the dimension of the coordinates
  CoordinateArraySequence.withDimensionMeasures(
      this.coordinates, int dimension, int measures) {
    this.dimension = dimension;
    this.measures = measures;
    if (coordinates == null) {
      this.coordinates = [];
    }
    enforceArrayConsistency(this.coordinates!);
  }

  /// Constructs a sequence of a given size, populated
  /// with {@link Coordinate}s.
  ///
  /// @param size the size of the sequence to create
  CoordinateArraySequence.fromSize(int size) {
    coordinates = []; //..length = (size);
    for (int i = 0; i < size; i++) {
      coordinates!.add(Coordinate.empty2D());
      // coordinates![i] = Coordinate.empty2D();
    }
  }

  /// Constructs a sequence of a given size, populated
  /// with {@link Coordinate}s.
  ///
  /// @param size the size of the sequence to create
  /// @param dimension the dimension of the coordinates
  CoordinateArraySequence.fromSizeDimension(int size, int dimension) {
    coordinates = []; //..length = (size);
    this.dimension = dimension;
    for (int i = 0; i < size; i++) {
      coordinates!.add(Coordinates.create(dimension));
      // coordinates![i] = Coordinates.create(dimension);
    }
  }

  /// Constructs a sequence of a given size, populated
  /// with {@link Coordinate}s.
  ///
  /// @param size the size of the sequence to create
  /// @param dimension the dimension of the coordinates
  CoordinateArraySequence.fromSizeDimensionMeasures(
      int size, int dimension, int measures) {
    coordinates = []; //..length = (size);
    this.dimension = dimension;
    this.measures = measures;
    for (int i = 0; i < size; i++) {
      coordinates!.add(createCoordinate());
      // coordinates![i] = createCoordinate();
    }
  }

  /// Creates a sequence based on a deep copy of the given {@link CoordinateSequence}.
  /// The coordinate dimension is set to equal the dimension of the input.
  ///
  /// @param coordSeq the coordinate sequence that will be copied.
  CoordinateArraySequence.fromSequence(CoordinateSequence? coordSeq) {
    // NOTE: this will make a sequence of the default dimension
    if (coordSeq == null) {
      coordinates = [];
      return;
    }
    dimension = coordSeq.getDimension();
    measures = coordSeq.getMeasures();
    coordinates = []; //..length = (coordSeq.size());

    for (int i = 0; i < coordSeq.size(); i++) {
      coordinates!.add(coordSeq.getCoordinateCopy(i));
      // coordinates![i] = coordSeq.getCoordinateCopy(i);
    }
  }

  /// Ensure array contents of the same type, making use of {@link #createCoordinate()} as needed.
  ///
  /// @param array array is modified in place as needed
  void enforceArrayConsistency(List<Coordinate?> array) {
    Coordinate sample = createCoordinate();

    String type = sample.runtimeType.toString();
    for (int i = 0; i < array.length; i++) {
      Coordinate? coordinate = array[i];
      if (coordinate == null) {
        array[i] = createCoordinate();
      } else if (coordinate.runtimeType.toString() != type) {
        Coordinate duplicate = createCoordinate();
        duplicate.setCoordinate(coordinate);
        array[i] = duplicate;
      }
    }
  }

  /// @see org.locationtech.jts.geom.CoordinateSequence#getDimension()
  int getDimension() {
    return dimension;
  }

  int getMeasures() {
    return measures;
  }

  /// Get the Coordinate with index i.
  ///
  /// @param i
  ///                  the index of the coordinate
  /// @return the requested Coordinate instance
  Coordinate getCoordinate(int i) {
    return coordinates![i];
  }

  /// Get a copy of the Coordinate with index i.
  ///
  /// @param i  the index of the coordinate
  /// @return a copy of the requested Coordinate
  Coordinate getCoordinateCopy(int i) {
    Coordinate copy = createCoordinate();
    copy.setCoordinate(coordinates![i]);
    return copy;
  }

  /// @see org.locationtech.jts.geom.CoordinateSequence#getX(int)
  void getCoordinateInto(int index, Coordinate coord) {
    coord.x = coordinates![index].x;
    coord.y = coordinates![index].y;
    if (hasZ()) {
      coord.z = coordinates![index].z;
    }
//    if (hasM()) {
//      coord.setMm = coordinates[index].m;
//    }
  }

  /// @see org.locationtech.jts.geom.CoordinateSequence#getX(int)
  double getX(int index) {
    return coordinates![index].x;
  }

  /// @see org.locationtech.jts.geom.CoordinateSequence#getY(int)
  double getY(int index) {
    return coordinates![index].y;
  }

  /// @see org.locationtech.jts.geom.CoordinateSequence#getZ(int)
  double getZ(int index) {
    if (hasZ()) {
      return coordinates![index].z;
    } else {
      return double.nan;
    }
  }

  /// @see org.locationtech.jts.geom.CoordinateSequence#getM(int)
  double getM(int index) {
    if (hasM()) {
      return coordinates![index].getM();
    } else {
      return double.nan;
    }
  }

  /// @see org.locationtech.jts.geom.CoordinateSequence#getOrdinate(int, int)
  double getOrdinate(int index, int ordinateIndex) {
    switch (ordinateIndex) {
      case CoordinateSequence.X:
        return coordinates![index].x;
      case CoordinateSequence.Y:
        return coordinates![index].y;
      default:
        return coordinates![index].getOrdinate(ordinateIndex);
    }
  }

  /// Creates a deep copy of the Object
  ///
  /// @return The deep copy
  /// @deprecated
  Object clone() {
    return copy();
  }

  /// Creates a deep copy of the CoordinateArraySequence
  ///
  /// @return The deep copy
  CoordinateArraySequence copy() {
    List<Coordinate> cloneCoordinates = []; //..length = (size());
    for (int i = 0; i < size(); i++) {
      Coordinate duplicate = createCoordinate();
      duplicate.setCoordinate(coordinates![i]);
      cloneCoordinates.add(duplicate);
      // cloneCoordinates[i] = duplicate;
    }
    return CoordinateArraySequence.withDimensionMeasures(
        cloneCoordinates, dimension, measures);
  }

  /// Returns the size of the coordinate sequence
  ///
  /// @return the number of coordinates
  int size() {
    return coordinates!.length;
  }

  /// @see org.locationtech.jts.geom.CoordinateSequence#setOrdinate(int, int, double)
  void setOrdinate(int index, int ordinateIndex, double value) {
    switch (ordinateIndex) {
      case CoordinateSequence.X:
        coordinates![index].x = value;
        break;
      case CoordinateSequence.Y:
        coordinates![index].y = value;
        break;
      default:
        coordinates![index].setOrdinate(ordinateIndex, value);
    }
  }

  /// This method exposes the internal Array of Coordinate Objects
  ///
  /// @return the List<Coordinate> array.
  List<Coordinate> toCoordinateArray() {
    return coordinates!;
  }

  Envelope expandEnvelope(Envelope env) {
    for (int i = 0; i < coordinates!.length; i++) {
      env.expandToIncludeCoordinate(coordinates![i]);
    }
    return env;
  }

  /// Returns the string Representation of the coordinate array
  ///
  /// @return The string
  String toString() {
    if (coordinates!.isNotEmpty) {
      StringBuffer strBuilder = StringBuffer();
      strBuilder.write('(');
      strBuilder.write(coordinates![0]);
      for (int i = 1; i < coordinates!.length; i++) {
        strBuilder.write(", ");
        strBuilder.write(coordinates![i]);
      }
      strBuilder.write(')');
      return strBuilder.toString();
    } else {
      return "()";
    }
  }
}

/// A {@link CoordinateSequence} implementation based on a packed arrays.
/// In this implementation, {@link Coordinate}s returned by #toArray and #get are copies
/// of the internal values.
/// To change the actual values, use the provided setters.
/// <p>
/// For efficiency, created Coordinate arrays
/// are cached using a soft reference.
/// The cache is cleared each time the coordinate sequence contents are
/// modified through a setter method.
///
/// @version 1.7
abstract class PackedCoordinateSequence extends CoordinateSequence {
  /// The dimensions of the coordinates held in the packed array
  int dimension;

  /// The number of measures of the coordinates held in the packed array.
  int measures;

  /// Creates an instance of this class
  /// @param dimension the total number of ordinates that make up a {@link Coordinate} in this sequence.
  /// @param measures the number of measure-ordinates each {@link Coordinate} in this sequence has.
  PackedCoordinateSequence(this.dimension, this.measures) {
    if (dimension - measures < 2) {
      throw ArgumentError("Must have at least 2 spatial dimensions");
    }
  }

  /// A soft reference to the List<Coordinate> representation of this sequence.
  /// Makes repeated coordinate array accesses more efficient.
  List<Coordinate>? coordRef;

  /// @see CoordinateSequence#getDimension()
  int getDimension() {
    return this.dimension;
  }

  /// @see CoordinateSequence#getMeasures()
  int getMeasures() {
    return this.measures;
  }

  /// @see CoordinateSequence#getCoordinate(int)
  Coordinate getCoordinate(int i) {
    List<Coordinate>? coords = getCachedCoords();
    if (coords != null)
      return coords[i];
    else
      return getCoordinateInternal(i);
  }

  /// @see CoordinateSequence#getCoordinate(int)
  Coordinate getCoordinateCopy(int i) {
    return getCoordinateInternal(i);
  }

  /// @see CoordinateSequence#getCoordinate(int)
  void getCoordinateInto(int i, Coordinate coord) {
    coord.x = getOrdinate(i, 0);
    coord.y = getOrdinate(i, 1);
    if (hasZ()) {
      coord.z = getZ(i);
    }
    if (hasM()) {
      coord.setM(getM(i));
    }
  }

  /// @see CoordinateSequence#toCoordinateArray()
  List<Coordinate> toCoordinateArray() {
    List<Coordinate>? coords = getCachedCoords();
// testing - never cache
    if (coords != null) return coords;

    coords = []; //..length = (size());
    for (int i = 0; i < size(); i++) {
      coords.add(getCoordinateInternal(i));
      // coords[i] = getCoordinateInternal(i);
    }
    coordRef = coords;

    return coords;
  }

  List<Coordinate>? getCachedCoords() {
    if (coordRef != null) {
      List<Coordinate>? coords = coordRef;
      if (coords != null) {
        return coords;
      } else {
        // System.out.print("-");
        coordRef = null;
        return null;
      }
    } else {
      // System.out.print("-");
      return null;
    }
  }

  /// @see CoordinateSequence#getX(int)
  double getX(int index) {
    return getOrdinate(index, 0);
  }

  /// @see CoordinateSequence#getY(int)
  double getY(int index) {
    return getOrdinate(index, 1);
  }

  /// @see CoordinateSequence#getOrdinate(int, int)
  double getOrdinate(int index, int ordinateIndex);

  /// Sets the first ordinate of a coordinate in this sequence.
  ///
  /// @param index  the coordinate index
  /// @param value  the new ordinate value
  void setX(int index, double value) {
    coordRef = null;
    setOrdinate(index, 0, value);
  }

  /// Sets the second ordinate of a coordinate in this sequence.
  ///
  /// @param index  the coordinate index
  /// @param value  the new ordinate value
  void setY(int index, double value) {
    coordRef = null;
    setOrdinate(index, 1, value);
  }

  String asString() {
    return CoordinateSequences.asString(this);
  }

  /// Returns a Coordinate representation of the specified coordinate, by always
  /// building a new Coordinate object
  ///
  /// @param index  the coordinate index
  /// @return  the {@link Coordinate} at the given index
  Coordinate getCoordinateInternal(int index);

  /// @see java.lang.Object#clone()
  /// @see CoordinateSequence#clone()
  /// @deprecated
  Object clone();

  PackedCoordinateSequence copy();

  /// Sets the ordinate of a coordinate in this sequence.
  /// <br>
  /// Warning: for performance reasons the ordinate index is not checked
  /// - if it is over dimensions you may not get an exception but a meaningless value.
  ///
  /// @param index
  ///          the coordinate index
  /// @param ordinate
  ///          the ordinate index in the coordinate, 0 based, smaller than the
  ///          number of dimensions
  /// @param value
  ///          the new ordinate value
  void setOrdinate(int index, int ordinate, double value);
}

/// Packed coordinate sequence implementation based on doubles
class Double extends PackedCoordinateSequence {
  /// The packed coordinate array
  late List<double> coords;

  /// Builds a new packed coordinate sequence
  ///
  /// @param coords  an array of <c>double</c> values that contains the ordinate values of the sequence
  /// @param dimension the total number of ordinates that make up a {@link Coordinate} in this sequence.
  /// @param measures the number of measure-ordinates each {@link Coordinate} in this sequence has.
  Double(List<double> coords, int dimension, int measures)
      : super(dimension, measures) {
    if (coords.length % dimension != 0) {
      throw ArgumentError("Packed array does not contain " +
          "an integral number of coordinates");
    }
    this.coords = coords;
  }

  /// Builds a new packed coordinate sequence out of a coordinate array
  ///
  /// @param coordinates an array of {@link Coordinate}s
  /// @param dimension the total number of ordinates that make up a {@link Coordinate} in this sequence.
  Double.fromCoordinatesDim(List<Coordinate> coordinates, int dimension)
      : super(dimension, 0);

  /// Builds a new packed coordinate sequence out of a coordinate array
  ///
  /// @param coordinates an array of {@link Coordinate}s
  /// @param dimension the total number of ordinates that make up a {@link Coordinate} in this sequence.
  /// @param measures the number of measure-ordinates each {@link Coordinate} in this sequence has.
  Double.fromCoordinatesDimMeas(
      List<Coordinate>? coordinates, int dimension, int measures)
      : super(dimension, measures) {
    if (coordinates == null) {
      coordinates = [];
    }

    coords = List.filled(coordinates.length * this.dimension, 0.0);
    for (int i = 0; i < coordinates.length; i++) {
      int offset = i * dimension;
      coords[offset] = coordinates[i].x;
      coords[offset + 1] = coordinates[i].y;
      if (dimension >= 3) {
        coords[offset + 2] = coordinates[i].getOrdinate(2);
      } // Z or M
      if (dimension >= 4) {
        coords[offset + 3] = coordinates[i].getOrdinate(3);
      } // M
    }
  }

  /// Builds a new packed coordinate sequence out of a coordinate array
  ///
  /// @param coordinates an array of {@link Coordinate}s
  Double.fromCoordinates(List<Coordinate> coordinates)
      : this.fromCoordinatesDimMeas(coordinates, 3, 0);

  /// Builds a new empty packed coordinate sequence of a given size and dimension
  ///
  /// @param size the number of coordinates in this sequence
  /// @param dimension the total number of ordinates that make up a {@link Coordinate} in this sequence.
  /// @param measures the number of measure-ordinates each {@link Coordinate} in this sequence has.
  Double.fromSizeDimMeas(int size, int dimension, int measures)
      : super(dimension, measures) {
    coords = List.filled(size * this.dimension, 0.0);
  }

  /// @see PackedCoordinateSequence#getCoordinate(int)
  Coordinate getCoordinateInternal(int i) {
    double x = coords[i * dimension];
    double y = coords[i * dimension + 1];
    if (dimension == 2 && measures == 0) {
      return new CoordinateXY.fromXY(x, y);
    } else if (dimension == 3 && measures == 0) {
      double z = coords[i * dimension + 2];
      return new Coordinate.fromXYZ(x, y, z);
    } else if (dimension == 3 && measures == 1) {
      double m = coords[i * dimension + 2];
      return new CoordinateXYM(x, y, m);
    } else if (dimension == 4 && measures == 1) {
      double z = coords[i * dimension + 2];
      double m = coords[i * dimension + 3];
      return new CoordinateXYZM(x, y, z, m);
    }
    return new CoordinateXY.fromXY(x, y);
  }

  /// Gets the underlying array containing the coordinate values.
  ///
  /// @return the array of coordinate values
  List<double?> getRawCoordinates() {
    return coords;
  }

  /// @see CoordinateSequence#size()
  int size() {
    return coords.length ~/ dimension;
  }

  /// @see java.lang.Object#clone()
  /// @see PackedCoordinateSequence#clone()
  /// @deprecated
  Object clone() {
    return copy();
  }

  /// @see PackedCoordinateSequence#size()
  Double copy() {
    List<double> clone = List.from(coords);
    return Double(clone, dimension, measures);
  }

  /// @see PackedCoordinateSequence#getOrdinate(int, int)
  ///      Beware, for performance reasons the ordinate index is not checked, if
  ///      it's over dimensions you may not get an exception but a meaningless
  ///      value.
  double getOrdinate(int index, int ordinate) {
    return coords[index * dimension + ordinate];
  }

  /// @see PackedCoordinateSequence#setOrdinate(int, int, double)
  void setOrdinate(int index, int ordinate, double value) {
    coordRef = null;
    coords[index * dimension + ordinate] = value;
  }

  /// @see CoordinateSequence#expandEnvelope(Envelope)
  Envelope expandEnvelope(Envelope env) {
    for (int i = 0; i < coords.length; i += dimension) {
      env.expandToInclude(coords[i], coords[i + 1]);
    }
    return env;
  }
}

/// Builds packed array coordinate sequences.
/// The array data type can be either
/// <code>double</code> or <code>float</code>,
/// and defaults to <code>double</code>.
class PackedCoordinateSequenceFactory extends CoordinateSequenceFactory {
  /// Type code for arrays of type <code>double</code>.
  static final int DOUBLE = 0;

  /// A factory using array type {@link #DOUBLE}
  static final PackedCoordinateSequenceFactory DOUBLE_FACTORY =
      PackedCoordinateSequenceFactory.withType(DOUBLE);

  static final int DEFAULT_MEASURES = 0;

  static final int DEFAULT_DIMENSION = 3;

  int type = DOUBLE;

  /// Creates a new PackedCoordinateSequenceFactory
  /// of type DOUBLE.
  PackedCoordinateSequenceFactory() : this.withType(DOUBLE);

  /// Creates a new PackedCoordinateSequenceFactory
  /// of the given type.
  /// Acceptable type values are
  /// {@linkplain PackedCoordinateSequenceFactory#FLOAT}or
  /// {@linkplain PackedCoordinateSequenceFactory#DOUBLE}
  PackedCoordinateSequenceFactory.withType(int type) {
    this.type = type;
  }

  /// Gets the type of packed coordinate sequence this factory builds, either
  /// {@linkplain PackedCoordinateSequenceFactory#FLOAT} or
  /// {@linkplain PackedCoordinateSequenceFactory#DOUBLE}
  ///
  /// @return the type of packed array built
  int getType() {
    return type;
  }

  /// @see CoordinateSequenceFactory#create(List<Coordinate>)
  CoordinateSequence create(List<Coordinate> coordinates) {
    int dimension = DEFAULT_DIMENSION;
    int measures = DEFAULT_MEASURES;
    if (coordinates != null &&
        coordinates.length > 0 &&
        coordinates[0] != null) {
      Coordinate first = coordinates[0];
      dimension = Coordinates.dimension(first);
      measures = Coordinates.measures(first);
    }
    if (type == DOUBLE) {
      return Double.fromCoordinatesDimMeas(coordinates, dimension, measures);
    } else {
      throw ArgumentError("Only PackedCoordinateSequence Double supported");
    }
  }

  /// @see CoordinateSequenceFactory#create(CoordinateSequence)
  CoordinateSequence createFromSequence(CoordinateSequence coordSeq) {
    int dimension = coordSeq.getDimension();
    int measures = coordSeq.getMeasures();
    if (type == DOUBLE) {
      return Double.fromCoordinatesDimMeas(
          coordSeq.toCoordinateArray(), dimension, measures);
    } else {
      throw ArgumentError("Only PackedCoordinateSequence Double supported");
    }
  }

  /// Creates a packed coordinate sequence of type {@link #DOUBLE}
  /// from the provided array
  /// using the given coordinate dimension and a measure count of 0.
  ///
  /// @param packedCoordinates the array containing coordinate values
  /// @param dimension the coordinate dimension
  /// @return a packed coordinate sequence of type {@link #DOUBLE}
  CoordinateSequence createFromDoublesDim(
      List<double> packedCoordinates, int dimension) {
    return createFromDoublesDimMeas(
        packedCoordinates, dimension, DEFAULT_MEASURES);
  }

  /// Creates a packed coordinate sequence of type {@link #DOUBLE}
  /// from the provided array
  /// using the given coordinate dimension and measure count.
  ///
  /// @param packedCoordinates the array containing coordinate values
  /// @param dimension the coordinate dimension
  /// @param measures the coordinate measure count
  /// @return a packed coordinate sequence of type {@link #DOUBLE}
  CoordinateSequence createFromDoublesDimMeas(
      List<double> packedCoordinates, int dimension, int measures) {
    if (type == DOUBLE) {
      return Double(packedCoordinates, dimension, measures);
    } else {
      throw ArgumentError("Only PackedCoordinateSequence Double supported");
    }
  }

  /// @see org.locationtech.jts.geom.CoordinateSequenceFactory#create(int, int)
  CoordinateSequence createSizeDim(int size, int dimension) {
    if (type == DOUBLE) {
      return Double.fromSizeDimMeas(size, dimension, DEFAULT_MEASURES);
    } else {
      throw ArgumentError("Only PackedCoordinateSequence Double supported");
    }
  }

  /// @see org.locationtech.jts.geom.CoordinateSequenceFactory#create(int, int, int)
  CoordinateSequence createFromSizeDimMeas(
      int size, int dimension, int measures) {
    if (type == DOUBLE) {
      return Double.fromSizeDimMeas(size, dimension, measures);
    } else {
      throw ArgumentError("Only PackedCoordinateSequence Double supported");
    }
  }
}
