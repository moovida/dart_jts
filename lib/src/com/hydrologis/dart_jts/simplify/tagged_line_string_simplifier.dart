part of dart_jts;

class TaggedLineStringSimplifier {
  ///
  /// Simplifies a TaggedLineString, preserving topology
  /// (in the sense that no new intersections are introduced).
  /// Uses the recursive Douglas-Peucker algorithm.
  ///
  /// @author Martin Davis
  /// @version 1.7
  ///
  LineIntersector li = new RobustLineIntersector();
  LineSegmentIndex inputIndex = new LineSegmentIndex();
  LineSegmentIndex outputIndex = new LineSegmentIndex();
  late TaggedLineString line;
  List<Coordinate> linePts = List.empty(growable: true);
  double distanceTolerance = 0.0;
  TaggedLineStringSimplifier(this.inputIndex, this.outputIndex);

  ///
  /// Sets the distance tolerance for the simplification.
  /// All vertices in the simplified geometry will be within this
  /// distance of the original geometry.
  ///
  /// @param distanceTolerance the approximation tolerance to use
  ///
  void setDistanceTolerance(double distanceTolerance) {
    this.distanceTolerance = distanceTolerance;
  }

  ///
  /// Simplifies the given {@link TaggedLineString}
  /// using the distance tolerance specified.
  ///
  /// @param line the linestring to simplify
  ///
  void simplify(TaggedLineString line) {
    this.line = line;
    linePts = line.getParentCoordinates();
    simplifySection(0, linePts.length - 1, 0);
  }

  void simplifySection(int i, int j, int depth) {
    depth += 1;
    List<int> sectionIndex = List.filled(2, 0);
    if ((i + 1) == j) {
      LineSegment newSeg = line.getSegment(i);
      line.addToResult(newSeg);
      // leave this segment in the input index, for efficiency
      return;
    }

    bool isValidToSimplify = true;

    ///
    /// Following logic ensures that there is enough points in the output line.
    /// If there is already more points than the minimum, there's nothing to check.
    /// Otherwise, if in the worst case there wouldn't be enough points,
    /// don't flatten this segment (which avoids the worst case scenario)
    ////
    if (line.getResultSize() < line.getMinimumSize()) {
      int worstCaseSize = depth + 1;
      if (worstCaseSize < line.getMinimumSize()) isValidToSimplify = false;
    }

    List<double> distance = List.empty(growable: true);
    int furthestPtIndex = findFurthestPoint(linePts, i, j, distance);
    // flattening must be less than distanceTolerance
    if (distance[0] > distanceTolerance) isValidToSimplify = false;
    // test if flattened section would cause intersection
    LineSegment candidateSeg = LineSegment.empty();
    candidateSeg.p0 = linePts[i];
    candidateSeg.p1 = linePts[j];
    sectionIndex[0]=i;
    sectionIndex[1]=j;
    if (hasBadIntersection(line, sectionIndex, candidateSeg)) {
      isValidToSimplify = false;
    }

    if (isValidToSimplify) {
      LineSegment newSeg = flatten(i, j);
      line.addToResult(newSeg);
      return;
    }
    simplifySection(i, furthestPtIndex, depth);
    simplifySection(furthestPtIndex, j, depth);
  }

  int findFurthestPoint(
      List<Coordinate> pts, int i, int j, List<double> maxDistance) {
    LineSegment seg = LineSegment.empty();
    seg.p0 = pts[i];
    seg.p1 = pts[j];
    double maxDist = -1.0;
    int maxIndex = i;
    for (int k = i + 1; k < j; k++) {
      Coordinate midPt = pts[k];
      double distance = seg.distanceCoord(midPt);
      if (distance > maxDist) {
        maxDist = distance;
        maxIndex = k;
      }
    }
    maxDistance.add(maxDist);
    return maxIndex;
  }

  ///
  /// Flattens a section of the line between
  /// indexes <code>start</code> and <code>end</code>,
  /// replacing them with a line between the endpoints.
  /// The input and output indexes are updated
  /// to reflect this.
  ///
  /// @param start the start index of the flattened section
  /// @param end the end index of the flattened section
  /// @return the new segment created
  ////
  LineSegment flatten(int start, int end) {
    // make a new segment for the simplified geometry
    Coordinate p0 = linePts[start];
    Coordinate p1 = linePts[end];
    LineSegment newSeg = LineSegment.fromCoordinates(p0, p1);
    // update the indexes
    remove(line, start, end);
    outputIndex.addSegment(newSeg);
    return newSeg;
  }

  bool hasBadIntersection(TaggedLineString parentLine, List<int> sectionIndex,
      LineSegment candidateSeg) {
    if (hasBadOutputIntersection(candidateSeg)) return true;
    if (hasBadInputIntersection(parentLine, sectionIndex, candidateSeg))
      return true;
    return false;
  }

  bool hasBadOutputIntersection(LineSegment candidateSeg) {
    List<dynamic> querySegs = outputIndex.query(candidateSeg);
    for (var segment in querySegs) {
      LineSegment querySeg = segment;
      if (hasInvalidIntersection(querySeg, candidateSeg)) {
        return true;
      }
    }
    return false;
  }

  bool hasBadInputIntersection(TaggedLineString parentLine,
      List<int> sectionIndex, LineSegment candidateSeg) {
    List<dynamic> querySegs = inputIndex.query(candidateSeg);
    for (var segment in querySegs) {
      TaggedLineSegment querySeg = segment;
      if (hasInvalidIntersection(querySeg, candidateSeg)) {
        //-- don't fail if the segment is part of parent line
        if (isInLineSection(parentLine, sectionIndex, querySeg)) continue;
        return true;
      }
    }
    return false;
  }

  ///
  /// Tests whether a segment is in a section of a TaggedLineString
  /// @param line
  /// @param sectionIndex
  /// @param seg
  /// @return
  ///
  static bool isInLineSection(
      TaggedLineString line, List<int> sectionIndex, TaggedLineSegment seg) {
    // not in this line
    if (seg.getParent() != line.getParent()) return false;
    int? segIndex = seg.getIndex();
    if (segIndex != null &&
        segIndex >= sectionIndex[0] &&
        segIndex < sectionIndex[1]) {
      return true;
    }
    return false;
  }

  bool hasInvalidIntersection(LineSegment seg0, LineSegment seg1) {
    //-- segments must not be equal
    if (seg0.equalsTopo(seg1)) return true;
    li.computeIntersection(seg0.p0, seg0.p1, seg1.p0, seg1.p1);
    return li.isInteriorIntersection();
  }

  ///
  /// Remove the segs in the section of the line
  /// @param line
  /// @param pts
  /// @param sectionStartIndex
  /// @param sectionEndIndex
  ///
  void remove(TaggedLineString line, int start, int end) {
    for (int i = start; i < end; i++) {
      TaggedLineSegment seg = line.getSegment(i);
      inputIndex.remove(seg);
    }
  }
}
