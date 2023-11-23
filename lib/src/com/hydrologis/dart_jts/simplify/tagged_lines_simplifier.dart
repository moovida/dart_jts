part of dart_jts;
class TaggedLinesSimplifier {
  ///
  /// Simplifies a collection of TaggedLineStrings, preserving topology
  /// (in the sense that no new intersections are introduced).
  /// This class is essentially just a container for the common
  /// indexes used by {@link TaggedLineStringSimplifier}.
  ///
  LineSegmentIndex inputIndex = LineSegmentIndex();
  LineSegmentIndex outputIndex = LineSegmentIndex();
  double distanceTolerance = 0.0;
  late TaggedLineString line;
  TaggedLinesSimplifier();
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
  /// Simplify a collection of TaggedLineStrings
  ///
  /// @param taggedLines the collection of lines to simplify
  ///
  void simplify(List<TaggedLineString> taggedLines) {
    for (TaggedLineString line in taggedLines) {
      inputIndex.add(line);
    }
    for (TaggedLineString line in taggedLines) {
    TaggedLineStringSimplifier tlss
    = new TaggedLineStringSimplifier(inputIndex, outputIndex);
    tlss.setDistanceTolerance(distanceTolerance);
    tlss.simplify(line);
    }
  }
}