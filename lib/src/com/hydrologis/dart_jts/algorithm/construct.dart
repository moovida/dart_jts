part of dart_jts;
class IndexedDistanceToPoint {

  Geometry targetGeometry;
  IndexedFacetDistance? facetDistance;
  IndexedPointInPolygonsLocator? ptLocater;

  IndexedDistanceToPoint(this.targetGeometry);


  void init() {
    if (facetDistance != null)
      return;
    facetDistance = new IndexedFacetDistance(targetGeometry);
    ptLocater = new IndexedPointInPolygonsLocator(targetGeometry);
  }

  /**
   * Computes the distance from a point to the geometry.
   *
   * @param pt the input point
   * @return the distance to the geometry
   */
  double distance(Point pt) {
    init();
    //-- distance is 0 if point is inside a target polygon
    if (isInArea(pt)) {
      return 0;
    }
    return facetDistance!.facetDistance(pt);
  }

  bool isInArea(Point pt) {
    return Location.EXTERIOR != ptLocater!.locate(pt.getCoordinate()!);
  }

  /**
   * Gets the nearest locations between the geometry and a point.
   * The first location lies on the geometry,
   * and the second location is the provided point.
   *
   * @param pt the point to compute the nearest location for
   * @return a pair of locations
   */
  List<Coordinate?>? nearestPoints(Point pt) {
    init();
    if (isInArea(pt)) {
      Coordinate? p = pt.getCoordinate();
      return [ p!.copy(), p.copy() ];
    }
    return facetDistance!.nearestPointsToGeometry(pt);
  }
}
/**
 * Determines the location of a point in the polygonal elements of a geometry.
 * Uses spatial indexing to provide efficient performance.
 *
 * @author mdavis
 *
 */
class IndexedPointInPolygonsLocator implements PointOnGeometryLocator {

  Geometry geometry;
  STRtree? index;

  IndexedPointInPolygonsLocator(this.geometry);

  void init() {
    if (index != null)
      return;
    List<Geometry> polys = PolygonExtracter.getPolygons(geometry).cast<Geometry>();
    index = new STRtree();
    for (int i = 0; i < polys.length; i++) {
      Geometry poly = polys[i];
      index!.insert(poly.getEnvelopeInternal(), new IndexedPointInAreaLocator(poly));
    }
  }

  @override
  int locate(Coordinate p) {
    init();

    List<IndexedPointInAreaLocator> results = index!.query(Envelope.fromCoordinate(p)).cast<IndexedPointInAreaLocator>();
    for (IndexedPointInAreaLocator ptLocater in results) {
      int loc = ptLocater.locate(p);
      if (loc != Location.EXTERIOR)
        return loc;
    }
    return Location.EXTERIOR;
  }
}
/**
 * Constructs the Largest Empty Circle for a set
 * of obstacle geometries, up to a given accuracy distance tolerance.
 * The obstacles may be any combination of point, linear and polygonal geometries.
 * <p>
 * The Largest Empty Circle (LEC) is the largest circle
 * whose interior does not intersect with any obstacle
 * and whose center lies within a polygonal boundary.
 * The circle center is the point in the interior of the boundary
 * which has the farthest distance from the obstacles
 * (up to the accuracy of the distance tolerance).
 * The circle itself is determined by the center point
 * and a point lying on an obstacle determining the circle radius.
 * <p>
 * The polygonal boundary may be supplied explicitly.
 * If it is not specified the convex hull of the obstacles is used as the boundary.
 * <p>
 * To compute an LEC which lies <i>wholly</i> within
 * a polygonal boundary, include the boundary of the polygon(s) as a linear obstacle.
 * <p>
 * The implementation uses a successive-approximation technique
 * over a grid of square cells covering the obstacles and boundary.
 * The grid is refined using a branch-and-bound algorithm.
 * Point containment and distance are computed in a performant
 * way by using spatial indexes.
 *
 * @author Martin Davis
 *
 * @see MaximumInscribedCircle
 * @see InteriorPoint
 * @see Centroid
 */
class LargestEmptyCircle {

  /**
   * Computes the center point of the Largest Empty Circle
   * interior-disjoint to a set of obstacles,
   * with accuracy to a given tolerance distance.
   * The obstacles may be any collection of points, lines and polygons.
   * The center of the LEC lies within the convex hull of the obstacles.
   *
   * @param obstacles a geometry representing the obstacles
   * @param tolerance the distance tolerance for computing the center point
   * @return the center point of the Largest Empty Circle
   */
  static Point getCenter(Geometry obstacles, double tolerance) {
    return getCenterBoundary(obstacles, null, tolerance);
  }

  /**
   * Computes the center point of the Largest Empty Circle
   * interior-disjoint to a set of obstacles and within a polygonal boundary,
   * with accuracy to a given tolerance distance.
   * The obstacles may be any collection of points, lines and polygons.
   * The center of the LEC lies within the given boundary.
   *
   * @param obstacles a geometry representing the obstacles
   * @param boundary a polygonal geometry to contain the LEC center
   * @param tolerance the distance tolerance for computing the center point
   * @return the center point of the Largest Empty Circle
   */
  static Point getCenterBoundary(Geometry obstacles, Geometry? boundary, double tolerance) {
    LargestEmptyCircle lec = new LargestEmptyCircle(obstacles, boundary!, tolerance);
    return lec.getCenterPoint();
  }

  /**
   * Computes a radius line of the Largest Empty Circle
   * interior-disjoint to a set of obstacles,
   * with accuracy to a given tolerance distance.
   * The obstacles may be any collection of points, lines and polygons.
   * The center of the LEC lies within the convex hull of the obstacles.
   *
   * @param obstacles a geometry representing the obstacles
   * @param tolerance the distance tolerance for computing the center point
   * @return a line from the center of the circle to a point on the edge
   */
  static LineString getRadiusLine(Geometry obstacles, double tolerance) {
    return getRadiusLineBoundary(obstacles, null, tolerance);
  }

  /**
   * Computes a radius line of the Largest Empty Circle
   * interior-disjoint to a set of obstacles and within a polygonal boundary,
   * with accuracy to a given tolerance distance.
   * The obstacles may be any collection of points, lines and polygons.
   * The center of the LEC lies within the given boundary.
   *
   * @param obstacles a geometry representing the obstacles
   * @param boundary a polygonal geometry to contain the LEC center
   * @param tolerance the distance tolerance for computing the center point
   * @return a line from the center of the circle to a point on the edge
   */
  static LineString getRadiusLineBoundary(Geometry obstacles, Geometry? boundary, double tolerance) {
    LargestEmptyCircle lec = LargestEmptyCircle(obstacles, boundary!, tolerance);
    return lec.getRadiusLineString();
  }

  Geometry obstacles;
  Geometry boundary;
  double tolerance;

  late GeometryFactory factory;
  late IndexedDistanceToPoint obstacleDistance;
  late IndexedPointInAreaLocator boundaryPtLocater;
  late IndexedFacetDistance boundaryDistance;
  late Envelope gridEnv;
  late CellLEC farthestCell;

  late CellLEC centerCell;
  late Coordinate centerPt;
  late Point centerPoint;
  late Coordinate radiusPt;
  late Point radiusPoint;
  late Geometry bounds;

  /**
   * Creates a new instance of a Largest Empty Circle construction,
   * interior-disjoint to a set of obstacle geometries
   * and having its center within a polygonal boundary.
   * The obstacles may be any collection of points, lines and polygons.
   * If the boundary is null or empty the convex hull
   * of the obstacles is used as the boundary.
   *
   * @param obstacles a non-empty geometry representing the obstacles
   * @param boundary a polygonal geometry (may be null or empty)
   * @param tolerance a distance tolerance for computing the circle center point (a positive value)
   */
  LargestEmptyCircle(this.obstacles, this.boundary, this.tolerance) {
    if (obstacles.isEmpty()) {
      throw new ArgumentError("Obstacles geometry is empty");
    }
    if (! (boundary is Polygonal)) {
      throw new ArgumentError("Boundary must be polygonal");
    }
    if (tolerance <= 0) {
      throw new ArgumentError('Accuracy tolerance is non-positive: $tolerance');
    }
    this.factory = obstacles.getFactory();
    obstacleDistance = new IndexedDistanceToPoint( obstacles );
  }

  /**
   * Gets the center point of the Largest Empty Circle
   * (up to the tolerance distance).
   *
   * @return the center point of the Largest Empty Circle
   */
  Point getCenterPoint() {
    compute();
    return centerPoint;
  }

  /**
   * Gets a point defining the radius of the Largest Empty Circle.
   * This is a point on the obstacles which is
   * nearest to the computed center of the Largest Empty Circle.
   * The line segment from the center to this point
   * is a radius of the constructed circle, and this point
   * lies on the boundary of the circle.
   *
   * @return a point defining the radius of the Largest Empty Circle
   */
  Point getRadiusPoint() {
    compute();
    return radiusPoint;
  }

  /**
   * Gets a line representing a radius of the Largest Empty Circle.
   *
   * @return a line from the center of the circle to a point on the edge
   */
  LineString getRadiusLineString() {
    compute();
    LineString radiusLine = factory.createLineString(
        [ centerPt.copy(), radiusPt.copy() ]);
    return radiusLine;
  }

  /**
   * Computes the signed distance from a point to the constraints
   * (obstacles and boundary).
   * Points outside the boundary polygon are assigned a negative distance.
   * Their containing cells will be last in the priority queue
   * (but will still end up being tested since they may be refined).
   *
   * @param p the point to compute the distance for
   * @return the signed distance to the constraints (negative indicates outside the boundary)
   */
  double distanceToConstraintsPoint(Point p) {
    bool isOutide = Location.EXTERIOR == boundaryPtLocater.locate(p.getCoordinate()!);
    if (isOutide) {
      double boundaryDist = boundaryDistance.facetDistance(p);
      return -boundaryDist;
    }
    double dist = obstacleDistance.distance(p);
    return dist;
  }

  double distanceToConstraintsXY(double x, double y) {
    Coordinate coord = new Coordinate(x, y);
    Point pt = factory.createPoint(coord);
    return distanceToConstraintsPoint(pt);
  }

  void initBoundary() {
    bounds = this.boundary;
    if (bounds.isEmpty()) {
      bounds = obstacles.convexHull();
    }
    //-- the centre point must be in the extent of the boundary
    gridEnv = bounds.getEnvelopeInternal();
    // if bounds does not enclose an area cannot create a ptLocater
    if (bounds.getDimension() >= 2) {
      boundaryPtLocater = new IndexedPointInAreaLocator( bounds );
      boundaryDistance = new IndexedFacetDistance( bounds );
    }
  }

  void compute() {
    initBoundary();

    // if boundaryPtLocater is not present then result is degenerate (represented as zero-radius circle)
    // if (boundaryPtLocater == null) {
    //   Coordinate? pt = obstacles.getCoordinate();
    //   centerPt = pt!.copy();
    //   centerPoint = factory.createPoint(pt);
    //   radiusPt = pt.copy();
    //   radiusPoint = factory.createPoint(pt);
    //   return;
    // }

    // Priority queue of cells, ordered by decreasing distance from constraints
    PriorityQueue cellQueue = PriorityQueue();

    //-- grid covers extent of obstacles and boundary (if any)
    createInitialGrid(gridEnv, cellQueue);

    // use the area centroid as the initial candidate center point
    farthestCell = createCentroidCell(obstacles);
    //int totalCells = cellQueue.size();

    /**
     * Carry out the branch-and-bound search
     * of the cell space
     */
    int maxIter = MaximumInscribedCircle.computeMaximumIterations(bounds, tolerance);
    int iter = 0;
    while (! cellQueue.isEmpty() && iter < maxIter) {
      iter++;
      // pick the cell with greatest distance from the queue
      CellLEC cell = cellQueue.poll() as CellLEC;
      //System.out.println(iter + "] Dist: " + cell.getDistance() + " Max D: " + cell.getMaxDistance() + " size: " + cell.getHSide());

      // update the center cell if the candidate is further from the constraints
      if (cell.getDistance() > farthestCell.getDistance()) {
        farthestCell = cell;
      }

      /**
       * If this cell may contain a better approximation to the center
       * of the empty circle, then refine it (partition into subcells
       * which are added into the queue for further processing).
       * Otherwise the cell is pruned (not investigated further),
       * since no point in it can be further than the current farthest distance.
       */
      if (mayContainCircleCenter(cell)) {
        // split the cell into four sub-cells
        double h2 = cell.getHSide() / 2;
        cellQueue.add( createCell( cell.getX() - h2, cell.getY() - h2, h2));
        cellQueue.add( createCell( cell.getX() + h2, cell.getY() - h2, h2));
        cellQueue.add( createCell( cell.getX() - h2, cell.getY() + h2, h2));
        cellQueue.add( createCell( cell.getX() + h2, cell.getY() + h2, h2));
        //totalCells += 4;
      }
    }
    // the farthest cell is the best approximation to the LEC center
    centerCell = farthestCell;
    // compute center point
    centerPt = new Coordinate(centerCell.getX(), centerCell.getY());
    centerPoint = factory.createPoint(centerPt);
    // compute radius point
    List<Coordinate?>? nearestPts = obstacleDistance.nearestPoints(centerPoint);
    if(nearestPts != null) {
      radiusPt = nearestPts[0]!.copy();
      radiusPoint = factory.createPoint(radiusPt);
    }
  }

  /**
   * Tests whether a cell may contain the circle center,
   * and thus should be refined (split into subcells
   * to be investigated further.)
   *
   * @param cell the cell to test
   * @return true if the cell might contain the circle center
   */
  bool mayContainCircleCenter(CellLEC cell) {
    /**
     * Every point in the cell lies outside the boundary,
     * so they cannot be the center point
     */
    if (cell.isFullyOutside())
      return false;

    /**
     * The cell is outside, but overlaps the boundary
     * so it may contain a point which should be checked.
     * This is only the case if the potential overlap distance
     * is larger than the tolerance.
     */
    if (cell.isOutside()) {
      bool isOverlapSignificant = cell.getMaxDistance() > tolerance;
      return isOverlapSignificant;
    }

    /**
     * Cell is inside the boundary. It may contain the center
     * if the maximum possible distance is greater than the current distance
     * (up to tolerance).
     */
    double potentialIncrease = cell.getMaxDistance() - farthestCell.getDistance();
    return potentialIncrease > tolerance;
  }

  /**
   * Initializes the queue with a cell covering
   * the extent of the area.
   *
   * @param env the area extent to cover
   * @param cellQueue the queue to initialize
   */
  void createInitialGrid(Envelope env, PriorityQueue cellQueue) {
    double cellSize = math.max(env.getWidth(), env.getHeight());
    double hSide = cellSize / 2.0;

    // Check for flat collapsed input and if so short-circuit
    // Result will just be centroid
    if (cellSize == 0) return;

    Coordinate? centre = env.centre();
    if(centre != null) {
      cellQueue.add(createCell(centre.x, centre.y, hSide));
    }
  }

  CellLEC createCell(double x, double y, double h) {
    return new CellLEC(x, y, h, distanceToConstraintsXY(x, y));
  }

  // create a cell centered on area centroid
  CellLEC createCentroidCell(Geometry geom) {
    Point p = geom.getCentroid();
    return CellLEC(p.getX(), p.getY(), 0, distanceToConstraintsPoint(p));
  }


}
class MaximumInscribedCircle {

  /**
   * Computes the center point of the Maximum Inscribed Circle
   * of a polygonal geometry, up to a given tolerance distance.
   *
   * @param polygonal a polygonal geometry
   * @param tolerance the distance tolerance for computing the center point
   * @return the center point of the maximum inscribed circle
   */
  static Point getCenter(Geometry polygonal, double tolerance) {
    MaximumInscribedCircle mic = new MaximumInscribedCircle(polygonal, tolerance);
    return mic.getCenterPoint();
  }

  /**
   * Computes a radius line of the Maximum Inscribed Circle
   * of a polygonal geometry, up to a given tolerance distance.
   *
   * @param polygonal a polygonal geometry
   * @param tolerance the distance tolerance for computing the center point
   * @return a line from the center to a point on the circle
   */
  static LineString getRadiusLine(Geometry polygonal, double tolerance) {
    MaximumInscribedCircle mic = new MaximumInscribedCircle(polygonal, tolerance);
    return mic.getRadiusLineString();
  }

  /**
   * Computes the maximum number of iterations allowed.
   * Uses a heuristic based on the size of the input geometry
   * and the tolerance distance.
   * A smaller tolerance distance allows more iterations.
   * This is a rough heuristic, intended
   * to prevent huge iterations for very thin geometries.
   *
   * @param geom the input geometry
   * @param toleranceDist the tolerance distance
   * @return the maximum number of iterations allowed
   */
  static int computeMaximumIterations(Geometry geom, double toleranceDist) {
    double diam = geom.getEnvelopeInternal().getDiameter();
    double ncells = diam / toleranceDist;
    //-- Using log of ncells allows control over number of iterations
    int factor = math.log(ncells).toInt();
    if (factor < 1) factor = 1;
    return 2000 + 2000 * factor;
  }

  late Geometry inputGeom;
  late double tolerance;

  late GeometryFactory factory;
  late IndexedPointInAreaLocator ptLocater;
  late IndexedFacetDistance indexedDistance;
  late CellMIC centerCell;
  late Coordinate centerPt;
  late Coordinate radiusPt;
  late Point centerPoint;
  late Point radiusPoint;

  /**
   * Creates a new instance of a Maximum Inscribed Circle computation.
   *
   * @param polygonal an areal geometry
   * @param tolerance the distance tolerance for computing the centre point (must be positive)
   * @throws IllegalArgumentException if the tolerance is non-positive, or the input geometry is non-polygonal or empty.
   */
  MaximumInscribedCircle(Geometry polygonal, double tolerance) {
    if (tolerance <= 0) {
      throw new ArgumentError("Tolerance must be positive");
    }
    if (! (polygonal is Polygon || polygonal is MultiPolygon)) {
      throw new ArgumentError("Input geometry must be a Polygon or MultiPolygon");
    }
    if (polygonal.isEmpty()) {
      throw new ArgumentError("Empty input geometry is not supported");
    }

    this.inputGeom = polygonal;
    this.factory = polygonal.getFactory();
    this.tolerance = tolerance;
    ptLocater = new IndexedPointInAreaLocator(polygonal);
    indexedDistance = new IndexedFacetDistance( polygonal.getBoundary() );
  }

  /**
   * Gets the center point of the maximum inscribed circle
   * (up to the tolerance distance).
   *
   * @return the center point of the maximum inscribed circle
   */
  Point getCenterPoint() {
    compute();
    return centerPoint;
  }

  /**
   * Gets a point defining the radius of the Maximum Inscribed Circle.
   * This is a point on the boundary which is
   * nearest to the computed center of the Maximum Inscribed Circle.
   * The line segment from the center to this point
   * is a radius of the constructed circle, and this point
   * lies on the boundary of the circle.
   *
   * @return a point defining the radius of the Maximum Inscribed Circle
   */
  Point getRadiusPoint() {
    compute();
    return radiusPoint;
  }

  /**
   * Gets a line representing a radius of the Largest Empty Circle.
   *
   * @return a line from the center of the circle to a point on the edge
   */
  LineString getRadiusLineString() {
    compute();
    LineString radiusLine = factory.createLineString(
        [ centerPt.copy(), radiusPt.copy() ]);
    return radiusLine;
  }

  /**
   * Computes the signed distance from a point to the area boundary.
   * Points outside the polygon are assigned a negative distance.
   * Their containing cells will be last in the priority queue
   * (but may still end up being tested since they may need to be refined).
   *
   * @param p the point to compute the distance for
   * @return the signed distance to the area boundary (negative indicates outside the area)
   */
  double distanceToBoundaryPoint(Point p) {
    double dist = indexedDistance.facetDistance(p);
    bool isOutide = Location.EXTERIOR == ptLocater.locate(p.getCoordinate()!);
    if (isOutide) return -dist;
    return dist;
  }

  double distanceToBoundaryXY(double x, double y) {
    Coordinate coord = new Coordinate(x, y);
    Point pt = factory.createPoint(coord);
    return distanceToBoundaryPoint(pt);
  }

  void compute() {
    // Priority queue of cells, ordered by maximum distance from boundary
    PriorityQueue cellQueue = PriorityQueue();

    createInitialGrid(inputGeom.getEnvelopeInternal(), cellQueue);

    // initial candidate center point
    CellMIC farthestCell = createInterorPointCell(inputGeom);
    //int totalCells = cellQueue.size();

    /**
     * Carry out the branch-and-bound search
     * of the cell space
     */
    int maxIter = computeMaximumIterations(inputGeom, tolerance);
    int iter = 0;
    while (! cellQueue.isEmpty() && iter < maxIter) {
      iter++;
      // pick the most promising cell from the queue
      CellMIC cell = cellQueue.poll() as CellMIC;

      //System.out.println(factory.toGeometry(cell.getEnvelope()));
      //System.out.println(iter + "] Dist: " + cell.getDistance() + " Max D: " + cell.getMaxDistance() + " size: " + cell.getHSide());
      //TestBuilderProxy.showIndicator(inputGeom.getFactory().toGeometry(cell.getEnvelope()));

      //-- if cell must be closer than furthest, terminate since all remaining cells in queue are even closer.
      if (cell.getMaxDistance() < farthestCell.getDistance())
        break;

      // update the circle center cell if the candidate is further from the boundary
      if (cell.getDistance() > farthestCell.getDistance()) {
        farthestCell = cell;
      }
      /**
       * Refine this cell if the potential distance improvement
       * is greater than the required tolerance.
       * Otherwise the cell is pruned (not investigated further),
       * since no point in it is further than
       * the current farthest distance.
       */
      double potentialIncrease = cell.getMaxDistance() - farthestCell.getDistance();
      if (potentialIncrease > tolerance) {
        // split the cell into four sub-cells
        double h2 = cell.getHSide() / 2;
        cellQueue.add( createCell( cell.getX() - h2, cell.getY() - h2, h2));
        cellQueue.add( createCell( cell.getX() + h2, cell.getY() - h2, h2));
        cellQueue.add( createCell( cell.getX() - h2, cell.getY() + h2, h2));
        cellQueue.add( createCell( cell.getX() + h2, cell.getY() + h2, h2));
        //totalCells += 4;
      }
    }
    // the farthest cell is the best approximation to the MIC center
    centerCell = farthestCell;
    centerPt = new Coordinate(centerCell.getX(), centerCell.getY());
    centerPoint = factory.createPoint(centerPt);
    // compute radius point
    List<Coordinate?>? nearestPts = indexedDistance.nearestPointsToGeometry(centerPoint);
    if(nearestPts != null) {
      radiusPt = nearestPts[0]!.copy();
      radiusPoint = factory.createPoint(radiusPt);
    }
  }

  /**
   * Initializes the queue with a cell covering
   * the extent of the area.
   *
   * @param env the area extent to cover
   * @param cellQueue the queue to initialize
   */
  void createInitialGrid(Envelope env, PriorityQueue cellQueue) {
    double cellSize = math.max(env.getWidth(), env.getHeight());
    double hSide = cellSize / 2.0;

    // Check for flat collapsed input and if so short-circuit
    // Result will just be centroid
    if (cellSize == 0) return;

    Coordinate? centre = env.centre();
    if(centre != null) {
      cellQueue.add(createCell(centre.x, centre.y, hSide));
    }
  }

  CellMIC createCell(double x, double y, double hSide) {
    return new CellMIC(x, y, hSide, distanceToBoundaryXY(x, y));
  }

  // create a cell at an interior point
  CellMIC createInterorPointCell(Geometry geom) {
    Point p = geom.getInteriorPoint();
    return new CellMIC(p.getX(), p.getY(), 0, distanceToBoundaryPoint(p));
  }


}
/**
 * A square grid cell centered on a given point
 * with a given side half-length,
 * and having a given distance from the center point to the constraints.
 * The maximum possible distance from any point in the cell to the
 * constraints can be computed.
 * This is used as the ordering and upper-bound function in
 * the branch-and-bound algorithm.
 */
class CellLEC implements Comparable<CellLEC> {

  static final double SQRT2 = 1.4142135623730951;

  double x;
  double y;
  double hSide;
  double distance;
  late double maxDist;

  CellLEC(this.x, this.y, this.hSide, this.distance) {
    this.x = x; // cell center x
    this.y = y; // cell center y
    this.hSide = hSide; // half the cell size

    /**
     * The maximum possible distance to the constraints for points in this cell
     * is the center distance plus the radius (half the diagonal length).
     */
    this.maxDist = distance + hSide * SQRT2;
  }

  bool isFullyOutside() {
    return getMaxDistance() < 0;
  }

  bool isOutside() {
    return distance < 0;
  }

  double getMaxDistance() {
    return maxDist;
  }

  double getDistance() {
    return distance;
  }

  double getHSide() {
    return hSide;
  }

  double getX() {
    return x;
  }

  double getY() {
    return y;
  }

    int compareTo(CellLEC o) {
    return -maxDist.compareTo(o.maxDist);
  }
}
/**
 * A square grid cell centered on a given point,
 * with a given half-side size, and having a given distance
 * to the area boundary.
 * The maximum possible distance from any point in the cell to the
 * boundary can be computed, and is used
 * as the ordering and upper-bound function in
 * the branch-and-bound algorithm.
 *
 */
class CellMIC implements Comparable<CellMIC> {

  static final double SQRT2 = 1.4142135623730951;

  double x;
  double y;
  double hSide;
  double distance;
  late double maxDist;

  CellMIC(this.x, this.y, this.hSide, this.distance) {
     // the maximum possible distance to area boundary for points in this cell
    this.maxDist = distance + hSide * SQRT2;
  }

  Envelope getEnvelope() {
    return Envelope(x - hSide, x + hSide, y - hSide, y + hSide);
  }

  double getMaxDistance() {
    return maxDist;
  }

  double getDistance() {
    return distance;
  }

  double getHSide() {
    return hSide;
  }

  double getX() {
    return x;
  }

  double getY() {
    return y;
  }

  /**
   * For maximum efficieny sort the PriorityQueue with largest maxDistance at front.
   * Since Java PQ sorts least-first, need to invert the comparison
   */
  @override
  int compareTo(CellMIC o) {
    return -maxDist.compareTo(o.maxDist);
  }

}
