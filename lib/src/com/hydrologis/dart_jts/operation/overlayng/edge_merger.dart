part of dart_jts;

/**
 * Performs merging on the noded edges of the input geometries.
 * Merging takes place on edges which are coincident
 * (i.e. have the same coordinate list, modulo direction).
 * The following situations can occur:
 * <ul>
 * <li>Coincident edges from different input geometries have their labels combined
 * <li>Coincident edges from the same area geometry indicate a topology collapse.
 * In this case the topology locations are "summed" to provide a final
 * assignment of side location
 * <li>Coincident edges from the same linear geometry can simply be merged
 * using the same ON location
 * </ul>
 *
 * The merging attempts to preserve the direction of linear
 * edges if possible (which is the case if there is
 * no other coincident edge, or if all coincident edges have the same direction).
 * This ensures that the overlay output line direction will be as consistent
 * as possible with input lines.
 * <p>
 * The merger also preserves the order of the edges in the input.
 * This means that for polygon-line overlay
 * the result lines will be in the same order as in the input
 * (possibly with multiple result lines for a single input line).
 *
 * @author mdavis
 *
 */
class EdgeMerger {
  static List<EdgeNG> merge(List<EdgeNG> edges) {
    // use a list to collect the final edges, to preserve order
    List<EdgeNG> mergedEdges = List.empty(growable: true);
    Map<EdgeKey, EdgeNG> edgeMap = {};

    for (EdgeNG edge in edges) {
      EdgeKey edgeKey = EdgeKey.create(edge);
      EdgeNG? baseEdge = edgeMap[edgeKey];
      if (baseEdge == null) {
        // this is the first (and maybe only) edge for this line
        edgeMap.putIfAbsent(edgeKey, () => edge);
        mergedEdges.add(edge);
      } else {
        // found an existing edge

        // Assert: edges are identical (up to direction)
        // this is a fast (but incomplete) sanity check
        Assert.isTrue(baseEdge.size() == edge.size(),
            "Merge of edges of different sizes - probable noding error.");

        baseEdge.merge(edge);
        //Debug.println("edge merged: " + existing);
        //Debug.println(edge.toLineString());
      }
    }
    return mergedEdges;
  }
}

/**
 * A key for sorting and comparing edges in a noded arrangement.
 * Relies on the fact that in a correctly noded arrangement
 * edges are identical (up to direction)
 * if they have their first segment in common.
 *
 * @author mdavis
 *
 */
class EdgeKey implements Comparable<EdgeKey> {
  static EdgeKey create(EdgeNG edge) {
    return EdgeKey(edge);
  }

  late double p0x;
  late double p0y;
  late double p1x;
  late double p1y;

  EdgeKey(EdgeNG edge) {
    initPoints(edge);
  }

  void initPoints(EdgeNG edge) {
    bool direction = edge.direction();
    if (direction) {
      _init(edge.getCoordinate(0), edge.getCoordinate(1));
    } else {
      int len = edge.size();
      _init(edge.getCoordinate(len - 1), edge.getCoordinate(len - 2));
    }
  }

  void _init(Coordinate p0, Coordinate p1) {
    p0x = p0.getX();
    p0y = p0.getY();
    p1x = p1.getX();
    p1y = p1.getY();
  }

  int compareToX(EdgeKey ek) {
    if (p0x < ek.p0x) return -1;
    if (p0x > ek.p0x) return 1;
    if (p0y < ek.p0y) return -1;
    if (p0y > ek.p0y) return 1;
    // first points are equal, compare second
    if (p1x < ek.p1x) return -1;
    if (p1x > ek.p1x) return 1;
    if (p1y < ek.p1y) return -1;
    if (p1y > ek.p1y) return 1;
    return 0;
  }

  @override
  bool operator ==(Object other) =>
      identical(this, other) ||
      other is EdgeKey &&
          runtimeType == other.runtimeType &&
          p0x == other.p0x &&
          p0y == other.p0y &&
          p1x == other.p1x &&
          p1y == other.p1y;

  @override
  int get hashCode => p0x.hashCode ^ p0y.hashCode ^ p1x.hashCode ^ p1y.hashCode;

  String toString() {
    return "EdgeKey(" + format(p0x, p0y) + ", " + format(p1x, p1y) + ")";
  }

  String format(double x, double y) {
    return x.toStringAsPrecision(9) + " " + y.toStringAsPrecision(9);
  }

  @override
  int compareTo(EdgeKey ek) {
    if (p0x < ek.p0x) return -1;
    if (p0x > ek.p0x) return 1;
    if (p0y < ek.p0y) return -1;
    if (p0y > ek.p0y) return 1;
    // first points are equal, compare second
    if (p1x < ek.p1x) return -1;
    if (p1x > ek.p1x) return 1;
    if (p1y < ek.p1y) return -1;
    if (p1y > ek.p1y) return 1;
    return 0;
  }
}
