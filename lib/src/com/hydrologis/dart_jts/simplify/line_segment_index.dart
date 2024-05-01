part of dart_jts;

class LineSegmentIndex {
  ///
  /// An spatial index on a set of {@link LineSegment}s.
  /// Supports adding and removing items.
  ///
  /// @author Martin Davis
  ///

  Quadtree index = Quadtree();
  LineSegmentIndex();
  void add(TaggedLineString line) {
    List<TaggedLineSegment> segs = line.getSegments();
    for (int i = 0; i < segs.length; i++) {
      TaggedLineSegment seg = segs[i];
      addSegment(seg);
    }
  }

  void addSegment(LineSegment seg) {
    index.insert(Envelope.fromCoordinates(seg.p0, seg.p1), seg);
  }

  void remove(LineSegment seg) {
    index.remove(Envelope.fromCoordinates(seg.p0, seg.p1), seg);
  }

  List query(LineSegment querySeg) {
    Envelope env = new Envelope.fromCoordinates(querySeg.p0, querySeg.p1);

    LineSegmentVisitor visitor = new LineSegmentVisitor(querySeg);
    index.queryWithVisitor(env, visitor);
    List itemsFound = visitor.getItems();
    return itemsFound;
  }
}

///
/// ItemVisitor subclass to reduce volume of query results.
///
class LineSegmentVisitor implements ItemVisitor {
  LineSegment querySeg;
  List<Object> items = List.empty(growable: true);

  LineSegmentVisitor(this.querySeg);

  void visitItem(Object item) {
    LineSegment seg = item as LineSegment;
    if (Envelope.intersectsEnvelopeCoords(
        seg.p0, seg.p1, querySeg.p0, querySeg.p1)) {
      items.add(item);
    }
  }

  List<Object> getItems() {
    return items;
  }
}
