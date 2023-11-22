part of dart_jts;

class TaggedLineSegment extends LineSegment {
  ///
  /// A {@link LineSegment} which is tagged with its location in a parent {@link Geometry}.
  /// Used to index the segments in a geometry and recover the segment locations
  /// from the index.
  ////
  Geometry? parent;
  int? index;
  TaggedLineSegment(Coordinate p0, Coordinate p1, Geometry? parent, int index) :
    super.fromCoordinates(p0, p1){
    this.parent = parent;
    this.index = index;
  }
  TaggedLineSegment.fromCoordinates(Coordinate p0,Coordinate p1):super.fromCoordinates(p0,p1){
    this.index = -1;
    this.parent = null;
  }
  Geometry? getParent() { return parent; }
  int? getIndex() { return index; }
}
