part of dart_jts;
/**
 * Validates that a collection of {@link SegmentString}s is correctly noded.
 * Indexing is used to improve performance.
 * By default validation stops after a single
 * non-noded intersection is detected.
 * Alternatively, it can be requested to detect all intersections
 * by using {@link #setFindAllIntersections(boolean)}.
 * <p>
 * The validator does not check for topology collapse situations
 * (e.g. where two segment strings are fully co-incident).
 * <p>
 * The validator checks for the following situations which indicated incorrect noding:
 * <ul>
 * <li>Proper intersections between segments (i.e. the intersection is interior to both segments)
 * <li>Intersections at an interior vertex (i.e. with an endpoint or another interior vertex)
 * </ul>
 * <p>
 * The client may either test the {@link #isValid()} condition,
 * or request that a suitable {@link TopologyException} be thrown.
 *
 * @version 1.7
 *
 * @see NodingIntersectionFinder
 */
class FastNodingValidator
{
  /**
   * Gets a list of all intersections found.
   * Intersections are represented as {@link Coordinate}s.
   * List is empty if none were found.
   *
   * @param segStrings a collection of SegmentStrings
   * @return a list of Coordinate
   */
  static List computeIntersections(List<SegmentString> segStrings)
  {
    FastNodingValidator nv = FastNodingValidator(segStrings);
    nv.setFindAllIntersections(true);
    nv.isIntersectionValid();
    return nv.getIntersections();
  }

  LineIntersector li = new RobustLineIntersector();

  List<SegmentString> segStrings;
  bool findAllIntersections = false;
  late NodingIntersectionFinder segInt;
  bool isValid = true;

  /**
   * Creates a new noding validator for a given set of linework.
   *
   * @param segStrings a collection of {@link SegmentString}s
   */
  FastNodingValidator(this.segStrings);

  void setFindAllIntersections(bool findAllIntersections)
  {
    this.findAllIntersections = findAllIntersections;
  }

  /**
   * Gets a list of all intersections found.
   * Intersections are represented as {@link Coordinate}s.
   * List is empty if none were found.
   *
   * @return a list of Coordinate
   */
  List getIntersections()
  {
    return segInt.getIntersections();
  }

  /**
   * Checks for an intersection and
   * reports if one is found.
   *
   * @return true if the arrangement contains an interior intersection
   */
  bool isIntersectionValid()
  {
    execute();
    return isValid;
  }

  /**
   * Returns an error message indicating the segments containing
   * the intersection.
   *
   * @return an error message documenting the intersection location
   */
  String getErrorMessage()
  {
    if (isValid) return "no intersections found";

    List<Coordinate> intSegs = segInt.getIntersectionSegments();
    return "found non-noded intersection between "
        + WKTWriter.toLineString([intSegs[0], intSegs[1]])
        + " and "
        + WKTWriter.toLineString([intSegs[2], intSegs[3]]);
  }

  /**
   * Checks for an intersection and throws
   * a TopologyException if one is found.
   *
   * @throws TopologyException if an intersection is found
   */
  void checkValid()
  {
    execute();
    if (! isValid)
      throw TopologyException(getErrorMessage());
  }

  void execute()
  {
    if (segInt.count() > 0)
      return;
    checkInteriorIntersections();
  }

  void checkInteriorIntersections()
  {
    /**
     * MD - It may even be reliable to simply check whether
     * end segments (of SegmentStrings) have an interior intersection,
     * since noding should have split any true interior intersections already.
     */
    isValid = true;
    segInt = NodingIntersectionFinder(li);
    segInt.setFindAllIntersections(findAllIntersections);
    MCIndexNoder noder = new MCIndexNoder(segInt);
    noder.setSegmentIntersector(segInt);
    noder.computeNodes(segStrings);
    if (segInt.hasIntersection()) {
      isValid = false;
      return;
    }
  }

}
