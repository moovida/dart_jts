part of dart_sfs;

/**
 * Represents a point on an
 * edge which intersects with another edge.
 * <p>
 * The intersection may either be a single point, or a line segment
 * (in which case this point is the start of the line segment)
 * The intersection point must be precise.
 *
 * @version 1.7
 */
class EdgeIntersection implements Comparable {
  Coordinate coord; // the point of intersection
  int segmentIndex; // the index of the containing line segment in the parent edge
  double dist; // the edge distance of this point along the containing line segment

  EdgeIntersection(Coordinate coord, int segmentIndex, double dist) {
    this.coord = new Coordinate.fromCoordinate(coord);
    this.segmentIndex = segmentIndex;
    this.dist = dist;
  }

  Coordinate getCoordinate() {
    return coord;
  }

  int getSegmentIndex() {
    return segmentIndex;
  }

  double getDistance() {
    return dist;
  }

  int compareTo(Object obj) {
    EdgeIntersection other = obj;
    return compare(other.segmentIndex, other.dist);
  }

  /**
   * @return -1 this EdgeIntersection is located before the argument location
   * @return 0 this EdgeIntersection is at the argument location
   * @return 1 this EdgeIntersection is located after the argument location
   */
  int compare(int segmentIndex, double dist) {
    if (this.segmentIndex < segmentIndex) return -1;
    if (this.segmentIndex > segmentIndex) return 1;
    if (this.dist < dist) return -1;
    if (this.dist > dist) return 1;
    return 0;
  }

  bool isEndPoint(int maxSegmentIndex) {
    if (segmentIndex == 0 && dist == 0.0) return true;
    if (segmentIndex == maxSegmentIndex) return true;
    return false;
  }

//   void print(PrintStream out)
//  {
//    out.print(coord);
//    out.print(" seg # = " + segmentIndex);
//    out.println(" dist = " + dist);
//  }
  String toString() {
    return "$coord seg # = $segmentIndex dist = $dist";
  }
}

/**
 * A list of edge intersections along an {@link Edge}.
 * Implements splitting an edge with intersections
 * into multiple resultant edges.
 *
 * @version 1.7
 */
class EdgeIntersectionList {
  // a Map <EdgeIntersection, EdgeIntersection>
  Map nodeMap = new SplayTreeMap();
  Edge edge; // the parent edge

  EdgeIntersectionList(Edge edge) {
    this.edge = edge;
  }

  /**
   * Adds an intersection into the list, if it isn't already there.
   * The input segmentIndex and dist are expected to be normalized.
   * @return the EdgeIntersection found or added
   */
  EdgeIntersection add(Coordinate intPt, int segmentIndex, double dist) {
    EdgeIntersection eiNew = new EdgeIntersection(intPt, segmentIndex, dist);
    EdgeIntersection ei = nodeMap[eiNew];
    if (ei != null) {
      return ei;
    }
    nodeMap[eiNew] = eiNew;
    return eiNew;
  }

  /**
   * Returns an iterator of {@link EdgeIntersection}s
   *
   * @return an Iterator of EdgeIntersections
   */
  Iterator iterator() {
    return nodeMap.values.iterator;
  }

  /**
   * Tests if the given point is an edge intersection
   *
   * @param pt the point to test
   * @return true if the point is an intersection
   */
  bool isIntersection(Coordinate pt) {
    for (Iterator it = iterator(); it.moveNext();) {
      EdgeIntersection ei = it.current;
      if (ei.coord.equals(pt)) return true;
    }
    return false;
  }

  /**
   * Adds entries for the first and last points of the edge to the list
   */
  void addEndpoints() {
    int maxSegIndex = edge.pts.length - 1;
    add(edge.pts[0], 0, 0.0);
    add(edge.pts[maxSegIndex], maxSegIndex, 0.0);
  }

  /**
   * Creates new edges for all the edges that the intersections in this
   * list split the parent edge into.
   * Adds the edges to the input list (this is so a single list
   * can be used to accumulate all split edges for a Geometry).
   *
   * @param edgeList a list of EdgeIntersections
   */
  void addSplitEdges(List edgeList) {
    // ensure that the list has entries for the first and last point of the edge
    addEndpoints();

    Iterator it = iterator();
    // there should always be at least two entries in the list
    EdgeIntersection eiPrev = it.current;
    while (it.moveNext()) {
      EdgeIntersection ei = it.current;
      Edge newEdge = createSplitEdge(eiPrev, ei);
      edgeList.add(newEdge);

      eiPrev = ei;
    }
  }

  /**
   * Create a new "split edge" with the section of points between
   * (and including) the two intersections.
   * The label for the new edge is the same as the label for the parent edge.
   */
  Edge createSplitEdge(EdgeIntersection ei0, EdgeIntersection ei1) {
//Debug.print("\ncreateSplitEdge"); Debug.print(ei0); Debug.print(ei1);
    int npts = ei1.segmentIndex - ei0.segmentIndex + 2;

    Coordinate lastSegStartPt = edge.pts[ei1.segmentIndex];
    // if the last intersection point is not equal to the its segment start pt,
    // add it to the points list as well.
    // (This check is needed because the distance metric is not totally reliable!)
    // The check for point equality is 2D only - Z values are ignored
    bool useIntPt1 = ei1.dist > 0.0 || !ei1.coord.equals2D(lastSegStartPt);
    if (!useIntPt1) {
      npts--;
    }

    List<Coordinate> pts = List(npts);
    int ipt = 0;
    pts[ipt++] = new Coordinate.fromCoordinate(ei0.coord);
    for (int i = ei0.segmentIndex + 1; i <= ei1.segmentIndex; i++) {
      pts[ipt++] = edge.pts[i];
    }
    if (useIntPt1) pts[ipt] = ei1.coord;
    return new Edge(pts, new Label.fromLabel(edge.label));
  }

//   void print(PrintStream out)
//  {
//    out.println("Intersections:");
//    for (Iterator it = iterator(); it.hasNext(); ) {
//      EdgeIntersection ei = (EdgeIntersection) it.next();
//      ei.print(out);
//    }
//  }
}

/**
 * A Depth object records the topological depth of the sides
 * of an Edge for up to two Geometries.
 * @version 1.7
 */
class Depth {
  static final int NULL_VALUE = -1;

  static int depthAtLocation(int location) {
    if (location == Location.EXTERIOR) return 0;
    if (location == Location.INTERIOR) return 1;
    return NULL_VALUE;
  }

  List<List<int>> depth = MatrixUtils.createMatrix(2, 3, 0);

  Depth() {
    // initialize depth array to a sentinel value
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 3; j++) {
        depth[i][j] = NULL_VALUE;
      }
    }
  }

  int getDepth(int geomIndex, int posIndex) {
    return depth[geomIndex][posIndex];
  }

  void setDepth(int geomIndex, int posIndex, int depthValue) {
    depth[geomIndex][posIndex] = depthValue;
  }

  int getLocation(int geomIndex, int posIndex) {
    if (depth[geomIndex][posIndex] <= 0) return Location.EXTERIOR;
    return Location.INTERIOR;
  }

  void add3(int geomIndex, int posIndex, int location) {
    if (location == Location.INTERIOR) depth[geomIndex][posIndex]++;
  }

  /**
   * A Depth object is null (has never been initialized) if all depths are null.
   */
  bool isNull() {
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 3; j++) {
        if (depth[i][j] != NULL_VALUE) return false;
      }
    }
    return true;
  }

  bool isNull1(int geomIndex) {
    return depth[geomIndex][1] == NULL_VALUE;
  }

  bool isNull2(int geomIndex, int posIndex) {
    return depth[geomIndex][posIndex] == NULL_VALUE;
  }

  void add(Label lbl) {
    for (int i = 0; i < 2; i++) {
      for (int j = 1; j < 3; j++) {
        int loc = lbl.getLocationWithPosIndex(i, j);
        if (loc == Location.EXTERIOR || loc == Location.INTERIOR) {
          // initialize depth if it is null, otherwise add this location value
          if (isNull2(i, j)) {
            depth[i][j] = depthAtLocation(loc);
          } else
            depth[i][j] += depthAtLocation(loc);
        }
      }
    }
  }

  int getDelta(int geomIndex) {
    return depth[geomIndex][Position.RIGHT] - depth[geomIndex][Position.LEFT];
  }

  /**
   * Normalize the depths for each geometry, if they are non-null.
   * A normalized depth
   * has depth values in the set { 0, 1 }.
   * Normalizing the depths
   * involves reducing the depths by the same amount so that at least
   * one of them is 0.  If the remaining value is &gt; 0, it is set to 1.
   */
  void normalize() {
    for (int i = 0; i < 2; i++) {
      if (!isNull1(i)) {
        int minDepth = depth[i][1];
        if (depth[i][2] < minDepth) minDepth = depth[i][2];

        if (minDepth < 0) minDepth = 0;
        for (int j = 1; j < 3; j++) {
          int newValue = 0;
          if (depth[i][j] > minDepth) newValue = 1;
          depth[i][j] = newValue;
        }
      }
    }
  }

  String toString() {
    return "A: ${depth[0][1]} ,${depth[0][2]} B: ${depth[1][1]},${depth[1][2]}";
  }
}


/**
 * @version 1.7
 */
 class Node
    extends GraphComponent
{
   Coordinate coord; // only non-null if this node is precise
   EdgeEndStar edges;

   Node(Coordinate coord, EdgeEndStar edges)
  {
    this.coord = coord;
    this.edges = edges;
    label = new Label.args2(0, Location.NONE);
  }

   Coordinate getCoordinate() { return coord; }
   EdgeEndStar getEdges() { return edges; }

  /**
   * Tests whether any incident edge is flagged as
   * being in the result.
   * This test can be used to determine if the node is in the result,
   * since if any incident edge is in the result, the node must be in the result as well.
   *
   * @return <code>true</code> if any incident edge in the in the result
   */
   bool isIncidentEdgeInResult()
  {
    for (Iterator it = getEdges().getEdges().iterator(); it.hasNext(); ) {
      DirectedEdge de = (DirectedEdge) it.next();
      if (de.getEdge().isInResult())
        return true;
    }
    return false;
  }

   bool isIsolated()
  {
    return (label.getGeometryCount() == 1);
  }
  /**
   * Basic nodes do not compute IMs
   */
   void computeIM(IntersectionMatrix im) {}
  /**
   * Add the edge to the list of edges at this node
   */
   void add(EdgeEnd e)
  {
    // Assert: start pt of e is equal to node point
    edges.insert(e);
    e.setNode(this);
  }

   void mergeLabelFromNode(Node n)
  {
    mergeLabel(n.label);
  }

  /**
   * To merge labels for two nodes,
   * the merged location for each LabelElement is computed.
   * The location for the corresponding node LabelElement is set to the result,
   * as long as the location is non-null.
   */

   void mergeLabel(Label label2)
  {
    for (int i = 0; i < 2; i++) {
      int loc = computeMergedLocation(label2, i);
      int thisLoc = label.getLocation(i);
      if (thisLoc == Location.NONE) label.setLocation(i, loc);
    }
  }

   void setLabelWithIndex(int argIndex, int onLocation)
  {
    if (label == null) {
      label = new Label.args2(argIndex, onLocation);
    }
    else
      label.setLocationWithIndex(argIndex, onLocation);
  }

  /**
   * Updates the label of a node to BOUNDARY,
   * obeying the mod-2 boundaryDetermination rule.
   */
   void setLabelBoundary(int argIndex)
  {
    if (label == null) return;

    // determine the current location for the point (if any)
    int loc = Location.NONE;
    if (label != null)
      loc = label.getLocation(argIndex);
    // flip the loc
    int newLoc;
    switch (loc) {
      case Location.BOUNDARY: newLoc = Location.INTERIOR; break;
      case Location.INTERIOR: newLoc = Location.BOUNDARY; break;
      default: newLoc = Location.BOUNDARY;  break;
    }
    label.setLocationWithIndex(argIndex, newLoc);
  }

  /**
   * The location for a given eltIndex for a node will be one
   * of { null, INTERIOR, BOUNDARY }.
   * A node may be on both the boundary and the interior of a geometry;
   * in this case, the rule is that the node is considered to be in the boundary.
   * The merged location is the maximum of the two input values.
   */
  int computeMergedLocation(Label label2, int eltIndex)
  {
    int loc = Location.NONE;
    loc = label.getLocation(eltIndex);
    if (! label2.isNull(eltIndex)) {
      int nLoc = label2.getLocation(eltIndex);
      if (loc != Location.BOUNDARY) loc = nLoc;
    }
    return loc;
  }

//   void print(PrintStream out)
//  {
//    out.println("node " + coord + " lbl: " + label);
//  }
}

/**
 * @version 1.7
 */
class Edge extends GraphComponent {
  /**
   * Updates an IM from the label for an edge.
   * Handles edges from both L and A geometries.
   */
  static void updateIMStatic(Label label, IntersectionMatrix im) {
    im.setAtLeastIfValid(label.getLocationWithPosIndex(0, Position.ON), label.getLocationWithPosIndex(1, Position.ON), 1);
    if (label.isArea()) {
      im.setAtLeastIfValid(label.getLocationWithPosIndex(0, Position.LEFT), label.getLocationWithPosIndex(1, Position.LEFT), 2);
      im.setAtLeastIfValid(label.getLocationWithPosIndex(0, Position.RIGHT), label.getLocationWithPosIndex(1, Position.RIGHT), 2);
    }
  }

  List<Coordinate> pts;
  Envelope env;
  EdgeIntersectionList eiList;
  String name;
  MonotoneChainEdge mce;
  bool _isIsolated = true;
  Depth depth = new Depth();
  int depthDelta = 0; // the change in area depth from the R to L side of this edge

  Edge(List<Coordinate> pts, Label label) {
    eiList = new EdgeIntersectionList(this);
    this.pts = pts;
    this.label = label;
  }

  Edge.fromList(List<Coordinate> pts) : this(pts, null);

  int getNumPoints() {
    return pts.length;
  }

  void setName(String name) {
    this.name = name;
  }

  List<Coordinate> getCoordinates() {
    return pts;
  }

  Coordinate getCoordinateWithIndex(int i) {
    return pts[i];
  }

  Coordinate getCoordinate() {
    if (pts.length > 0) return pts[0];
    return null;
  }

  Envelope getEnvelope() {
    // compute envelope lazily
    if (env == null) {
      env = new Envelope.empty();
      for (int i = 0; i < pts.length; i++) {
        env.expandToIncludeCoordinate(pts[i]);
      }
    }
    return env;
  }

  Depth getDepth() {
    return depth;
  }

  /**
   * The depthDelta is the change in depth as an edge is crossed from R to L
   * @return the change in depth as the edge is crossed from R to L
   */
  int getDepthDelta() {
    return depthDelta;
  }

  void setDepthDelta(int depthDelta) {
    this.depthDelta = depthDelta;
  }

  int getMaximumSegmentIndex() {
    return pts.length - 1;
  }

  EdgeIntersectionList getEdgeIntersectionList() {
    return eiList;
  }

  MonotoneChainEdge getMonotoneChainEdge() {
    if (mce == null) mce = new MonotoneChainEdge(this);
    return mce;
  }

  bool isClosed() {
    return pts[0].equals(pts[pts.length - 1]);
  }

  /**
   * An Edge is collapsed if it is an Area edge and it consists of
   * two segments which are equal and opposite (eg a zero-width V).
   */
  bool isCollapsed() {
    if (!label.isArea()) return false;
    if (pts.length != 3) return false;
    if (pts[0].equals(pts[2])) return true;
    return false;
  }

  Edge getCollapsedEdge() {
    List<Coordinate> newPts = [];
    newPts[0] = pts[0];
    newPts[1] = pts[1];
    Edge newe = new Edge(newPts, Label.toLineLabel(label));
    return newe;
  }

  void setIsolated(bool isIsolated) {
    this._isIsolated = isIsolated;
  }

  bool isIsolated() {
    return _isIsolated;
  }

  /**
   * Adds EdgeIntersections for one or both
   * intersections found for a segment of an edge to the edge intersection list.
   */
  void addIntersections(LineIntersector li, int segmentIndex, int geomIndex) {
    for (int i = 0; i < li.getIntersectionNum(); i++) {
      addIntersection(li, segmentIndex, geomIndex, i);
    }
  }

  /**
   * Add an EdgeIntersection for intersection intIndex.
   * An intersection that falls exactly on a vertex of the edge is normalized
   * to use the higher of the two possible segmentIndexes
   */
  void addIntersection(LineIntersector li, int segmentIndex, int geomIndex, int intIndex) {
    Coordinate intPt = new Coordinate.fromCoordinate(li.getIntersection(intIndex));
    int normalizedSegmentIndex = segmentIndex;
    double dist = li.getEdgeDistance(geomIndex, intIndex);
//Debug.println("edge intpt: " + intPt + " dist: " + dist);
    // normalize the intersection point location
    int nextSegIndex = normalizedSegmentIndex + 1;
    if (nextSegIndex < pts.length) {
      Coordinate nextPt = pts[nextSegIndex];
//Debug.println("next pt: " + nextPt);

      // Normalize segment index if intPt falls on vertex
      // The check for point equality is 2D only - Z values are ignored
      if (intPt.equals2D(nextPt)) {
//Debug.println("normalized distance");
        normalizedSegmentIndex = nextSegIndex;
        dist = 0.0;
      }
    }
    /**
     * Add the intersection point to edge intersection list.
     */
    EdgeIntersection ei = eiList.add(intPt, normalizedSegmentIndex, dist);
//ei.print(System.out);
  }

  /**
   * Update the IM with the contribution for this component.
   * A component only contributes if it has a labelling for both parent geometries
   */
  void computeIM(IntersectionMatrix im) {
    updateIMStatic(label, im);
  }

  /**
   * equals is defined to be:
   * <p>
   * e1 equals e2
   * <b>iff</b>
   * the coordinates of e1 are the same or the reverse of the coordinates in e2
   */
  bool equals(Object o) {
    if (!(o is Edge)) return false;
    Edge e = o;

    if (pts.length != e.pts.length) return false;

    bool isEqualForward = true;
    bool isEqualReverse = true;
    int iRev = pts.length;
    for (int i = 0; i < pts.length; i++) {
      if (!pts[i].equals2D(e.pts[i])) {
        isEqualForward = false;
      }
      if (!pts[i].equals2D(e.pts[--iRev])) {
        isEqualReverse = false;
      }
      if (!isEqualForward && !isEqualReverse) return false;
    }
    return true;
  }

  /**
   * @return true if the coordinate sequences of the Edges are identical
   */
  bool isPointwiseEqual(Edge e) {
    if (pts.length != e.pts.length) return false;

    for (int i = 0; i < pts.length; i++) {
      if (!pts[i].equals2D(e.pts[i])) {
        return false;
      }
    }
    return true;
  }

  String toString() {
    StringBuffer builder = new StringBuffer();
    builder.write("edge " + name + ": ");
    builder.write("LINESTRING (");
    for (int i = 0; i < pts.length; i++) {
      if (i > 0) builder.write(",");
      builder.write(pts[i].x);
      builder.write(" ");
      builder.write(pts[i].y);
    }
    builder.write(")  $label $depthDelta");
    return builder.toString();
  }

//   void print(PrintStream out)
//  {
//    out.print("edge " + name + ": ");
//    out.print("LINESTRING (");
//    for (int i = 0; i < pts.length; i++) {
//      if (i > 0) out.print(",");
//      out.print(pts[i].x + " " + pts[i].y);
//    }
//    out.print(")  " + label + " " + depthDelta);
//  }
//   void printReverse(PrintStream out)
//  {
//    out.print("edge " + name + ": ");
//    for (int i = pts.length - 1; i >= 0; i--) {
//      out.print(pts[i] + " ");
//    }
//    out.println("");
//  }

}

/**
 * A TopologyLocation is the labelling of a
 * GraphComponent's topological relationship to a single Geometry.
 * <p>
 * If the parent component is an area edge, each side and the edge itself
 * have a topological location.  These locations are named
 * <ul>
 * <li> ON: on the edge
 * <li> LEFT: left-hand side of the edge
 * <li> RIGHT: right-hand side
 * </ul>
 * If the parent component is a line edge or node, there is a single
 * topological relationship attribute, ON.
 * <p>
 * The possible values of a topological location are
 * {Location.NONE, Location.EXTERIOR, Location.BOUNDARY, Location.INTERIOR}
 * <p>
 * The labelling is stored in an array location[j] where
 * where j has the values ON, LEFT, RIGHT
 * @version 1.7
 */
class TopologyLocation {
  List<int> location;

  TopologyLocation.fromList(List<int> location) {
    init(location.length);
  }

  /**
   * Constructs a TopologyLocation specifying how points on, to the left of, and to the
   * right of some GraphComponent relate to some Geometry. Possible values for the
   * parameters are Location.NULL, Location.EXTERIOR, Location.BOUNDARY,
   * and Location.INTERIOR.
   * @see Location
   */
  TopologyLocation(int on, int left, int right) {
    init(3);
    location[Position.ON] = on;
    location[Position.LEFT] = left;
    location[Position.RIGHT] = right;
  }

  TopologyLocation.fromOn(int on) {
    init(1);
    location[Position.ON] = on;
  }

  TopologyLocation.fromTL(TopologyLocation gl) {
    init(gl.location.length);
    if (gl != null) {
      for (int i = 0; i < location.length; i++) {
        location[i] = gl.location[i];
      }
    }
  }

  void init(int size) {
    location = List(size);
    setAllLocations(Location.NONE);
  }

  int get(int posIndex) {
    if (posIndex < location.length) return location[posIndex];
    return Location.NONE;
  }

  /**
   * @return true if all locations are NULL
   */
  bool isNull() {
    for (int i = 0; i < location.length; i++) {
      if (location[i] != Location.NONE) return false;
    }
    return true;
  }

  /**
   * @return true if any locations are NULL
   */
  bool isAnyNull() {
    for (int i = 0; i < location.length; i++) {
      if (location[i] == Location.NONE) return true;
    }
    return false;
  }

  bool isEqualOnSide(TopologyLocation le, int locIndex) {
    return location[locIndex] == le.location[locIndex];
  }

  bool isArea() {
    return location.length > 1;
  }

  bool isLine() {
    return location.length == 1;
  }

  void flip() {
    if (location.length <= 1) return;
    int temp = location[Position.LEFT];
    location[Position.LEFT] = location[Position.RIGHT];
    location[Position.RIGHT] = temp;
  }

  void setAllLocations(int locValue) {
    for (int i = 0; i < location.length; i++) {
      location[i] = locValue;
    }
  }

  void setAllLocationsIfNull(int locValue) {
    for (int i = 0; i < location.length; i++) {
      if (location[i] == Location.NONE) location[i] = locValue;
    }
  }

  void setLocationWithIndex(int locIndex, int locValue) {
    location[locIndex] = locValue;
  }

  void setLocation(int locValue) {
    setLocationWithIndex(Position.ON, locValue);
  }

  List<int> getLocations() {
    return location;
  }

  void setLocations(int on, int left, int right) {
    location[Position.ON] = on;
    location[Position.LEFT] = left;
    location[Position.RIGHT] = right;
  }

  bool allPositionsEqual(int loc) {
    for (int i = 0; i < location.length; i++) {
      if (location[i] != loc) return false;
    }
    return true;
  }

  /**
   * merge updates only the NULL attributes of this object
   * with the attributes of another.
   */
  void merge(TopologyLocation gl) {
    // if the src is an Area label & and the dest is not, increase the dest to be an Area
    if (gl.location.length > location.length) {
      List<int> newLoc = List(3);
      newLoc[Position.ON] = location[Position.ON];
      newLoc[Position.LEFT] = Location.NONE;
      newLoc[Position.RIGHT] = Location.NONE;
      location = newLoc;
    }
    for (int i = 0; i < location.length; i++) {
      if (location[i] == Location.NONE && i < gl.location.length) location[i] = gl.location[i];
    }
  }

  String toString() {
    StringBuffer buf = new StringBuffer();
    if (location.length > 1) buf.write(Location.toLocationSymbol(location[Position.LEFT]));
    buf.write(Location.toLocationSymbol(location[Position.ON]));
    if (location.length > 1) buf.write(Location.toLocationSymbol(location[Position.RIGHT]));
    return buf.toString();
  }
}

/**
 * A Position indicates the position of a Location relative to a graph component
 * (Node, Edge, or Area).
 * @version 1.7
 */
class Position {
  /** An indicator that a Location is <i>on</i> a GraphComponent */
  static final int ON = 0;

  /** An indicator that a Location is to the <i>left</i> of a GraphComponent */
  static final int LEFT = 1;

  /** An indicator that a Location is to the <i>right</i> of a GraphComponent */
  static final int RIGHT = 2;

  /**
   * Returns LEFT if the position is RIGHT, RIGHT if the position is LEFT, or the position
   * otherwise.
   */
  static int opposite(int position) {
    if (position == LEFT) return RIGHT;
    if (position == RIGHT) return LEFT;
    return position;
  }
}

/**
 * A <code>Label</code> indicates the topological relationship of a component
 * of a topology graph to a given <code>Geometry</code>.
 * This class supports labels for relationships to two <code>Geometry</code>s,
 * which is sufficient for algorithms for binary operations.
 * <P>
 * Topology graphs support the concept of labeling nodes and edges in the graph.
 * The label of a node or edge specifies its topological relationship to one or
 * more geometries.  (In fact, since JTS operations have only two arguments labels
 * are required for only two geometries).  A label for a node or edge has one or
 * two elements, depending on whether the node or edge occurs in one or both of the
 * input <code>Geometry</code>s.  Elements contain attributes which categorize the
 * topological location of the node or edge relative to the parent
 * <code>Geometry</code>; that is, whether the node or edge is in the interior,
 * boundary or exterior of the <code>Geometry</code>.  Attributes have a value
 * from the set <code>{Interior, Boundary, Exterior}</code>.  In a node each
 * element has  a single attribute <code>&lt;On&gt;</code>.  For an edge each element has a
 * triplet of attributes <code>&lt;Left, On, Right&gt;</code>.
 * <P>
 * It is up to the client code to associate the 0 and 1 <code>TopologyLocation</code>s
 * with specific geometries.
 * @version 1.7
 *
 */
class Label {
  // converts a Label to a Line label (that is, one with no side Locations)
  static Label toLineLabel(Label label) {
    Label lineLabel = new Label(Location.NONE);
    for (int i = 0; i < 2; i++) {
      lineLabel.setLocationWithIndex(i, label.getLocation(i));
    }
    return lineLabel;
  }

  List<TopologyLocation> elt = List(2);

  /**
   * Construct a Label with a single location for both Geometries.
   * Initialize the locations to Null
   */
  Label(int onLoc) {
    elt[0] = new TopologyLocation.fromOn(onLoc);
    elt[1] = new TopologyLocation.fromOn(onLoc);
  }

  /**
   * Construct a Label with a single location for both Geometries.
   * Initialize the location for the Geometry index.
   */
  Label.args2(int geomIndex, int onLoc) {
    elt[0] = new TopologyLocation.fromOn(Location.NONE);
    elt[1] = new TopologyLocation.fromOn(Location.NONE);
    elt[geomIndex].setLocation(onLoc);
  }

  /**
   * Construct a Label with On, Left and Right locations for both Geometries.
   * Initialize the locations for both Geometries to the given values.
   */
  Label.args3(int onLoc, int leftLoc, int rightLoc) {
    elt[0] = new TopologyLocation(onLoc, leftLoc, rightLoc);
    elt[1] = new TopologyLocation(onLoc, leftLoc, rightLoc);
  }

  /**
   * Construct a Label with On, Left and Right locations for both Geometries.
   * Initialize the locations for the given Geometry index.
   */
  Label.args4(int geomIndex, int onLoc, int leftLoc, int rightLoc) {
    elt[0] = new TopologyLocation(Location.NONE, Location.NONE, Location.NONE);
    elt[1] = new TopologyLocation(Location.NONE, Location.NONE, Location.NONE);
    elt[geomIndex].setLocations(onLoc, leftLoc, rightLoc);
  }

  /**
   * Construct a Label with the same values as the argument Label.
   */
  Label.fromLabel(Label lbl) {
    elt[0] = new TopologyLocation.fromTL(lbl.elt[0]);
    elt[1] = new TopologyLocation.fromTL(lbl.elt[1]);
  }

  void flip() {
    elt[0].flip();
    elt[1].flip();
  }

  int getLocationWithPosIndex(int geomIndex, int posIndex) {
    return elt[geomIndex].get(posIndex);
  }

  int getLocation(int geomIndex) {
    return elt[geomIndex].get(Position.ON);
  }

  void setLocation(int geomIndex, int posIndex, int location) {
    elt[geomIndex].setLocationWithIndex(posIndex, location);
  }

  void setLocationWithIndex(int geomIndex, int location) {
    elt[geomIndex].setLocationWithIndex(Position.ON, location);
  }

  void setAllLocations(int geomIndex, int location) {
    elt[geomIndex].setAllLocations(location);
  }

  void setAllLocationsIfNullWithIndex(int geomIndex, int location) {
    elt[geomIndex].setAllLocationsIfNull(location);
  }

  void setAllLocationsIfNull(int location) {
    setAllLocationsIfNullWithIndex(0, location);
    setAllLocationsIfNullWithIndex(1, location);
  }

  /**
   * Merge this label with another one.
   * Merging updates any null attributes of this label with the attributes from lbl
   */
  void merge(Label lbl) {
    for (int i = 0; i < 2; i++) {
      if (elt[i] == null && lbl.elt[i] != null) {
        elt[i] = new TopologyLocation.fromTL(lbl.elt[i]);
      } else {
        elt[i].merge(lbl.elt[i]);
      }
    }
  }

  int getGeometryCount() {
    int count = 0;
    if (!elt[0].isNull()) count++;
    if (!elt[1].isNull()) count++;
    return count;
  }

  bool isNull(int geomIndex) {
    return elt[geomIndex].isNull();
  }

  bool isAnyNull(int geomIndex) {
    return elt[geomIndex].isAnyNull();
  }

  bool isArea() {
    return elt[0].isArea() || elt[1].isArea();
  }

  bool isAreaWithIndex(int geomIndex) {
    /*  Testing
  	if (elt[0].getLocations().length != elt[1].getLocations().length) {
  		System.out.println(this);
  	}
  		*/
    return elt[geomIndex].isArea();
  }

  bool isLine(int geomIndex) {
    return elt[geomIndex].isLine();
  }

  bool isEqualOnSide(Label lbl, int side) {
    return this.elt[0].isEqualOnSide(lbl.elt[0], side) && this.elt[1].isEqualOnSide(lbl.elt[1], side);
  }

  bool allPositionsEqual(int geomIndex, int loc) {
    return elt[geomIndex].allPositionsEqual(loc);
  }

  /**
   * Converts one GeometryLocation to a Line location
   */
  void toLine(int geomIndex) {
    if (elt[geomIndex].isArea()) elt[geomIndex] = new TopologyLocation.fromOn(elt[geomIndex].location[0]);
  }

  String toString() {
    StringBuffer buf = new StringBuffer();
    if (elt[0] != null) {
      buf.write("A:");
      buf.write(elt[0].toString());
    }
    if (elt[1] != null) {
      buf.write(" B:");
      buf.write(elt[1].toString());
    }
    return buf.toString();
  }
}

/**
 * A GraphComponent is the parent class for the objects'
 * that form a graph.  Each GraphComponent can carry a
 * Label.
 * @version 1.7
 */
abstract class GraphComponent {
  Label label;

  /**
   * isInResult indicates if this component has already been included in the result
   */
  bool _isInResult = false;
  bool _isCovered = false;
  bool _isCoveredSet = false;
  bool _isVisited = false;

  GraphComponent() {}

  GraphComponent.fromLabel(Label label) {
    this.label = label;
  }

  Label getLabel() {
    return label;
  }

  void setLabel(Label label) {
    this.label = label;
  }

  void setInResult(bool isInResult) {
    this._isInResult = isInResult;
  }

  bool isInResult() {
    return _isInResult;
  }

  void setCovered(bool isCovered) {
    this._isCovered = isCovered;
    this._isCoveredSet = true;
  }

  bool isCovered() {
    return _isCovered;
  }

  bool isCoveredSet() {
    return _isCoveredSet;
  }

  bool isVisited() {
    return _isVisited;
  }

  void setVisited(bool isVisited) {
    this._isVisited = isVisited;
  }

  /**
   * @return a coordinate in this component (or null, if there are none)
   */
  Coordinate getCoordinate();

  /**
   * compute the contribution to an IM for this component
   */
  void computeIM(IntersectionMatrix im);

  /**
   * An isolated component is one that does not intersect or touch any other
   * component.  This is the case if the label has valid locations for
   * only a single Geometry.
   *
   * @return true if this component is isolated
   */
  bool isIsolated();

  /**
   * Update the IM with the contribution for this component.
   * A component only contributes if it has a labelling for both parent geometries
   */
  void updateIM(IntersectionMatrix im) {
//  TODO  Assert.isTrue(label.getGeometryCount() >= 2, "found partial label");
    computeIM(im);
  }
}
