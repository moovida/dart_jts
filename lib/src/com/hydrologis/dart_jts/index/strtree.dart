part of dart_jts;

/**
 * A function method which computes the distance
 * between two {@link ItemBoundable}s in an {@link STRtree}.
 * Used for Nearest Neighbour searches.
 * <p>
 * To make a distance function suitable for
 * querying a single index tree
 * via {@link STRtree#nearestNeighbour(ItemDistance)} ,
 * the function should have a non-zero <i>reflexive distance</i>.
 * That is, if the two arguments are the same object,
 * the distance returned should be non-zero.
 * If it is required that only pairs of <b>distinct</b> items be returned,
 * the distance function must be <i>anti-reflexive</i>,
 * and must return {@link Double#MAX_VALUE} for identical arguments.
 *
 * @author Martin Davis
 *
 */
abstract class ItemDistance {
  /**
   * Computes the distance between two items.
   *
   * @param item1
   * @param item2
   * @return the distance between the items
   *
   * @throws IllegalArgumentException if the metric is not applicable to the arguments
   */
  double distance(ItemBoundable item1, ItemBoundable item2);
}

/**
 * Functions for computing distances between {@link Envelope}s.
 *
 * @author mdavis
 *
 */
class EnvelopeDistance {
  /**
   * Computes the maximum distance between the points defining two envelopes.
   * It is equal to the length of the diagonal of
   * the envelope containing both input envelopes.
   * This is a coarse upper bound on the distance between
   * geometries bounded by the envelopes.
   *
   * @param env1 an envelope
   * @param env2 an envelope
   * @return the maximum distance between the points defining the envelopes
   */
  static double maximumDistance(Envelope env1, Envelope env2) {
    double minx = math.min(env1.getMinX(), env2.getMinX());
    double miny = math.min(env1.getMinY(), env2.getMinY());
    double maxx = math.max(env1.getMaxX(), env2.getMaxX());
    double maxy = math.max(env1.getMaxY(), env2.getMaxY());
    return distance(minx, miny, maxx, maxy);
  }

  static double distance(double x1, double y1, double x2, double y2) {
    double dx = x2 - x1;
    double dy = y2 - y1;
    return math.sqrt(dx * dx + dy * dy);
  }

  /**
   * Computes the Min-Max Distance between two {@link Envelope}s.
   * It is equal to the minimum of the maximum distances between all pairs of
   * edge segments from the two envelopes.
   * This is the tight upper bound on the distance between
   * geometric items bounded by the envelopes.
   * <p>
   * Theoretically this bound can be used in the R-tree nearest-neighbour branch-and-bound search
   * instead of {@link #maximumDistance(Envelope, Envelope)}.
   * However, little performance improvement is observed in practice.
   *
   * @param a an envelope
   * @param b an envelope
   * @return the min-max-distance between the envelopes
   */
  static double minMaxDistance(Envelope a, Envelope b) {
    double aminx = a.getMinX();
    double aminy = a.getMinY();
    double amaxx = a.getMaxX();
    double amaxy = a.getMaxY();
    double bminx = b.getMinX();
    double bminy = b.getMinY();
    double bmaxx = b.getMaxX();
    double bmaxy = b.getMaxY();

    double dist =
        maxDistance(aminx, aminy, aminx, amaxy, bminx, bminy, bminx, bmaxy);
    dist = math.min(dist,
        maxDistance(aminx, aminy, aminx, amaxy, bminx, bminy, bmaxx, bminy));
    dist = math.min(dist,
        maxDistance(aminx, aminy, aminx, amaxy, bmaxx, bmaxy, bminx, bmaxy));
    dist = math.min(dist,
        maxDistance(aminx, aminy, aminx, amaxy, bmaxx, bmaxy, bmaxx, bminy));

    dist = math.min(dist,
        maxDistance(aminx, aminy, amaxx, aminy, bminx, bminy, bminx, bmaxy));
    dist = math.min(dist,
        maxDistance(aminx, aminy, amaxx, aminy, bminx, bminy, bmaxx, bminy));
    dist = math.min(dist,
        maxDistance(aminx, aminy, amaxx, aminy, bmaxx, bmaxy, bminx, bmaxy));
    dist = math.min(dist,
        maxDistance(aminx, aminy, amaxx, aminy, bmaxx, bmaxy, bmaxx, bminy));

    dist = math.min(dist,
        maxDistance(amaxx, amaxy, aminx, amaxy, bminx, bminy, bminx, bmaxy));
    dist = math.min(dist,
        maxDistance(amaxx, amaxy, aminx, amaxy, bminx, bminy, bmaxx, bminy));
    dist = math.min(dist,
        maxDistance(amaxx, amaxy, aminx, amaxy, bmaxx, bmaxy, bminx, bmaxy));
    dist = math.min(dist,
        maxDistance(amaxx, amaxy, aminx, amaxy, bmaxx, bmaxy, bmaxx, bminy));

    dist = math.min(dist,
        maxDistance(amaxx, amaxy, amaxx, aminy, bminx, bminy, bminx, bmaxy));
    dist = math.min(dist,
        maxDistance(amaxx, amaxy, amaxx, aminy, bminx, bminy, bmaxx, bminy));
    dist = math.min(dist,
        maxDistance(amaxx, amaxy, amaxx, aminy, bmaxx, bmaxy, bminx, bmaxy));
    dist = math.min(dist,
        maxDistance(amaxx, amaxy, amaxx, aminy, bmaxx, bmaxy, bmaxx, bminy));

    return dist;
  }

  /**
   * Computes the maximum distance between two line segments.
   *
   * @param ax1 x ordinate of first endpoint of segment 1
   * @param ay1 y ordinate of first endpoint of segment 1
   * @param ax2 x ordinate of second endpoint of segment 1
   * @param ay2 y ordinate of second endpoint of segment 1
   * @param bx1 x ordinate of first endpoint of segment 2
   * @param by1 y ordinate of first endpoint of segment 2
   * @param bx2 x ordinate of second endpoint of segment 2
   * @param by2 y ordinate of second endpoint of segment 2
   * @return maximum distance between the segments
   */
  static double maxDistance(double ax1, double ay1, double ax2, double ay2,
      double bx1, double by1, double bx2, double by2) {
    double dist = distance(ax1, ay1, bx1, by1);
    dist = math.max(dist, distance(ax1, ay1, bx2, by2));
    dist = math.max(dist, distance(ax2, ay2, bx1, by1));
    dist = math.max(dist, distance(ax2, ay2, bx2, by2));
    return dist;
  }
}

/**
 * A pair of {@link Boundable}s, whose leaf items
 * support a distance metric between them.
 * Used to compute the distance between the members,
 * and to expand a member relative to the other
 * in order to produce new branches of the
 * Branch-and-Bound evaluation tree.
 * Provides an ordering based on the distance between the members,
 * which allows building a priority queue by minimum distance.
 *
 * @author Martin Davis
 *
 */
class BoundablePair implements Comparable {
  Boundable boundable1;
  Boundable boundable2;
  late double _distance;
  ItemDistance itemDistance;

  // double maxDistance = -1.0;

  BoundablePair(this.boundable1, this.boundable2, this.itemDistance) {
    this.boundable1 = boundable1;
    this.boundable2 = boundable2;
    this.itemDistance = itemDistance;
    _distance = distance();
  }

  /**
   * Gets one of the member {@link Boundable}s in the pair
   * (indexed by [0, 1]).
   *
   * @param i the index of the member to return (0 or 1)
   * @return the chosen member
   */
  Boundable getBoundable(int i) {
    if (i == 0) return boundable1;
    return boundable2;
  }

  /**
   * Computes the maximum distance between any
   * two items in the pair of nodes.
   *
   * @return the maximum distance between items in the pair
   */
  double maximumDistance() {
    return EnvelopeDistance.maximumDistance(
        boundable1.getBounds() as Envelope, boundable2.getBounds() as Envelope);
  }

  /**
   * Computes the distance between the {@link Boundable}s in this pair.
   * The boundables are either composites or leaves.
   * If either is composite, the distance is computed as the minimum distance
   * between the bounds.
   * If both are leaves, the distance is computed by {@link #itemDistance(ItemBoundable, ItemBoundable)}.
   *
   * @return
   */
  double distance() {
    // if items, compute exact distance
    if (isLeaves()) {
      return itemDistance.distance(
          boundable1 as ItemBoundable, boundable2 as ItemBoundable);
    }
    // otherwise compute distance between bounds of boundables
    return (boundable1.getBounds() as Envelope)
        .distance((boundable2.getBounds() as Envelope));
  }

  /**
   * Gets the minimum possible distance between the Boundables in
   * this pair.
   * If the members are both items, this will be the
   * exact distance between them.
   * Otherwise, this distance will be a lower bound on
   * the distances between the items in the members.
   *
   * @return the exact or lower bound distance for this pair
   */
  double getDistance() {
    return _distance;
  }

  /**
   * Compares two pairs based on their minimum distances
   */
  int compareTo(dynamic o) {
    BoundablePair nd = o as BoundablePair;
    if (_distance < nd._distance) return -1;
    if (_distance > nd._distance) return 1;
    return 0;
  }

  /**
   * Tests if both elements of the pair are leaf nodes
   *
   * @return true if both pair elements are leaf nodes
   */
  bool isLeaves() {
    return !(isComposite(boundable1) || isComposite(boundable2));
  }

  static bool isComposite(Object item) {
    return (item is AbstractNode);
  }

  static double area(Boundable b) {
    return (b.getBounds() as Envelope).getArea();
  }

  /**
   * For a pair which is not a leaf
   * (i.e. has at least one composite boundable)
   * computes a list of new pairs
   * from the expansion of the larger boundable
   * with distance less than minDistance
   * and adds them to a priority queue.
   * <p>
   * Note that expanded pairs may contain
   * the same item/node on both sides.
   * This must be allowed to support distance
   * functions which have non-zero distances
   * between the item and itself (non-zero reflexive distance).
   *
   * @param priQ the priority queue to add the new pairs to
   * @param minDistance the limit on the distance between added pairs
   *
   */
  void expandToQueue(PriorityQueue priQ, double minDistance) {
    bool isComp1 = isComposite(boundable1);
    bool isComp2 = isComposite(boundable2);

    /**
     * HEURISTIC: If both boundable are composite,
     * choose the one with largest area to expand.
     * Otherwise, simply expand whichever is composite.
     */
    if (isComp1 && isComp2) {
      if (area(boundable1) > area(boundable2)) {
        expand(boundable1, boundable2, false, priQ, minDistance);
        return;
      } else {
        expand(boundable2, boundable1, true, priQ, minDistance);
        return;
      }
    } else if (isComp1) {
      expand(boundable1, boundable2, false, priQ, minDistance);
      return;
    } else if (isComp2) {
      expand(boundable2, boundable1, true, priQ, minDistance);
      return;
    }

    throw new ArgumentError("neither boundable is composite");
  }

  void expand(Boundable bndComposite, Boundable bndOther, bool isFlipped,
      PriorityQueue priQ, double minDistance) {
    List children = (bndComposite as AbstractNode).getChildBoundables();
    for (Iterator i = children.iterator; i.moveNext();) {
      Boundable child = i.current as Boundable;
      BoundablePair bp;
      if (isFlipped) {
        bp = new BoundablePair(bndOther, child, itemDistance);
      } else {
        bp = new BoundablePair(child, bndOther, itemDistance);
      }
      // only add to queue if this pair might contain the closest points
      // MD - it's actually faster to construct the object rather than called distance(child, bndOther)!
      if (bp.getDistance() < minDistance) {
        priQ.add(bp);
      }
    }
  }
}

/**
 * A spatial object in an AbstractSTRtree.
 *
 * @version 1.7
 */
abstract class Boundable {
  /**
   * Returns a representation of space that encloses this Boundable, preferably
   * not much bigger than this Boundable's boundary yet fast to test for intersection
   * with the bounds of other Boundables. The class of object returned depends
   * on the subclass of AbstractSTRtree.
   * @return an Envelope (for STRtrees), an Interval (for SIRtrees), or other object
   * (for other subclasses of AbstractSTRtree)
   */
  Object getBounds();
}

/**
 * Boundable wrapper for a non-Boundable spatial object. Used internally by
 * AbstractSTRtree.
 *
 * @version 1.7
 */
class ItemBoundable implements Boundable {
  Object bounds;
  Object item;

  ItemBoundable(this.bounds, this.item);

  Object getBounds() {
    return bounds;
  }

  Object getItem() {
    return item;
  }
}

/**
 * A node of an {@link AbstractSTRtree}. A node is one of:
 * <ul>
 * <li>empty
 * <li>an <i>interior node</i> containing child {@link AbstractNode}s
 * <li>a <i>leaf node</i> containing data items ({@link ItemBoundable}s).
 * </ul>
 * A node stores the bounds of its children, and its level within the index tree.
 *
 * @version 1.7
 */
abstract class AbstractNode implements Boundable {
  /**
   *
   */
  List childBoundables = [];
  Object? bounds = null;
  int level = 0;

  /**
   * Default constructor required for serialization.
   */
  AbstractNode() {}

  /**
   * Constructs an AbstractNode at the given level in the tree
   * @param level 0 if this node is a leaf, 1 if a parent of a leaf, and so on; the
   * root node will have the highest level
   */
  AbstractNode.withLevel(int level) {
    this.level = level;
  }

  /**
   * Returns either child {@link AbstractNode}s, or if this is a leaf node, real data (wrapped
   * in {@link ItemBoundable}s).
   *
   * @return a list of the children
   */
  List getChildBoundables() {
    return childBoundables;
  }

  /**
   * Returns a representation of space that encloses this Boundable,
   * preferably not much bigger than this Boundable's boundary yet fast to
   * test for intersection with the bounds of other Boundables. The class of
   * object returned depends on the subclass of AbstractSTRtree.
   *
   * @return an Envelope (for STRtrees), an Interval (for SIRtrees), or other
   *         object (for other subclasses of AbstractSTRtree)
   * @see AbstractSTRtree.IntersectsOp
   */
  Object computeBounds();

  /**
   * Gets the bounds of this node
   *
   * @return the object representing bounds in this index
   */
  Object getBounds() {
    if (bounds == null) {
      bounds = computeBounds();
    }
    return bounds!;
  }

  /**
   * Returns 0 if this node is a leaf, 1 if a parent of a leaf, and so on; the
   * root node will have the highest level
   *
   * @return the node level
   */
  int getLevel() {
    return level;
  }

  /**
   * Gets the count of the {@link Boundable}s at this node.
   *
   * @return the count of boundables at this node
   */
  int size() {
    return childBoundables.length;
  }

  /**
   * Tests whether there are any {@link Boundable}s at this node.
   *
   * @return true if there are boundables at this node
   */
  bool isEmpty() {
    return childBoundables.isEmpty;
  }

  /**
   * Adds either an AbstractNode, or if this is a leaf node, a data object
   * (wrapped in an ItemBoundable)
   *
   * @param childBoundable the child to add
   */
  void addChildBoundable(Boundable childBoundable) {
    Assert.isTrue(bounds == null);
    childBoundables.add(childBoundable);
  }
}

/**
 * A test for intersection between two bounds, necessary because subclasses
 * of AbstractSTRtree have different implementations of bounds.
 */
abstract class IntersectsOp {
  /**
   * For STRtrees, the bounds will be Envelopes; for SIRtrees, Intervals;
   * for other subclasses of AbstractSTRtree, some other class.
   * @param aBounds the bounds of one spatial object
   * @param bBounds the bounds of another spatial object
   * @return whether the two bounds intersect
   */
  bool intersects(Object aBounds, Object bBounds);
}

class EnvelopeIntersectsOp implements IntersectsOp {
  bool intersects(Object aBounds, Object bBounds) {
    return (aBounds as Envelope).intersectsEnvelope(bBounds as Envelope);
  }
}

/**
 * Base class for STRtree and SIRtree. STR-packed R-trees are described in:
 * P. Rigaux, Michel Scholl and Agnes Voisard. <i>Spatial Databases With
 * Application To GIS.</i> Morgan Kaufmann, San Francisco, 2002.
 * <p>
 * This implementation is based on {@link Boundable}s rather than {@link AbstractNode}s,
 * because the STR algorithm operates on both nodes and
 * data, both of which are treated as Boundables.
 * <p>
 * This class is thread-safe.  Building the tree is synchronized,
 * and querying is stateless.
 *
 * @see STRtree
 * @see SIRtree
 *
 * @version 1.7
 */
abstract class AbstractSTRtree {
  late AbstractNode root;

  bool built = false;

  /**
   * Set to <tt>null</tt> when index is built, to avoid retaining memory.
   */
  List? itemBoundables = [];

  int nodeCapacity = 0;

  static final int DEFAULT_NODE_CAPACITY = 10;

  /**
   * Constructs an AbstractSTRtree with the
   * default node capacity.
   */
  AbstractSTRtree() : this.withCapacity(DEFAULT_NODE_CAPACITY);

  /**
   * Constructs an AbstractSTRtree with the specified maximum number of child
   * nodes that a node may have
   *
   * @param nodeCapacity the maximum number of child nodes in a node
   */
  AbstractSTRtree.withCapacity(int nodeCapacity) {
    Assert.isTrue(nodeCapacity > 1, "Node capacity must be greater than 1");
    this.nodeCapacity = nodeCapacity;
  }

  /**
   * Creates parent nodes, grandparent nodes, and so forth up to the root
   * node, for the data that has been inserted into the tree. Can only be
   * called once, and thus can be called only after all of the data has been
   * inserted into the tree.
   */
  // TODO check how to make this method synchronized
  void build() {
    if (built) return;
    root = itemBoundables!.isEmpty
        ? createNode(0)
        : createHigherLevels(itemBoundables!, -1);
    // the item list is no longer needed
    itemBoundables = null;
    built = true;
  }

  AbstractNode createNode(int level);

  /**
   * Sorts the childBoundables then divides them into groups of size M, where
   * M is the node capacity.
   */
  List createParentBoundables(List childBoundables, int newLevel) {
    Assert.isTrue(!childBoundables.isEmpty);
    List parentBoundables = [];
    parentBoundables.add(createNode(newLevel));
    List sortedChildBoundables = List.from(childBoundables);
    sortedChildBoundables.sort(getComparator());
    for (Iterator i = sortedChildBoundables.iterator; i.moveNext();) {
      Boundable childBoundable = i.current as Boundable;
      if (lastNode(parentBoundables).getChildBoundables().length ==
          getNodeCapacity()) {
        parentBoundables.add(createNode(newLevel));
      }
      lastNode(parentBoundables).addChildBoundable(childBoundable);
    }
    return parentBoundables;
  }

  AbstractNode lastNode(List nodes) {
    return nodes[nodes.length - 1] as AbstractNode;
  }

  static int compareDoubles(double a, double b) {
    return a > b
        ? 1
        : a < b
            ? -1
            : 0;
  }

  /**
   * Creates the levels higher than the given level
   *
   * @param boundablesOfALevel
   *            the level to build on
   * @param level
   *            the level of the Boundables, or -1 if the boundables are item
   *            boundables (that is, below level 0)
   * @return the root, which may be a ParentNode or a LeafNode
   */
  AbstractNode createHigherLevels(List boundablesOfALevel, int level) {
    Assert.isTrue(!boundablesOfALevel.isEmpty);
    List parentBoundables =
        createParentBoundables(boundablesOfALevel, level + 1);
    if (parentBoundables.length == 1) {
      return parentBoundables[0] as AbstractNode;
    }
    return createHigherLevels(parentBoundables, level + 1);
  }

  /**
   * Gets the root node of the tree.
   *
   * @return the root node
   */
  AbstractNode getRoot() {
    build();
    return root;
  }

  /**
   * Returns the maximum number of child nodes that a node may have.
   *
   * @return the node capacity
   */
  int getNodeCapacity() {
    return nodeCapacity;
  }

  /**
   * Tests whether the index contains any items.
   * This method does not build the index,
   * so items can still be inserted after it has been called.
   *
   * @return true if the index does not contain any items
   */
  bool isEmpty() {
    if (!built) return itemBoundables!.isEmpty;
    return root.isEmpty();
  }

  int size() {
    if (isEmpty()) {
      return 0;
    }
    build();
    return sizeWithNode(root);
  }

  int sizeWithNode(AbstractNode node) {
    int size = 0;
    for (Iterator i = node.getChildBoundables().iterator; i.moveNext();) {
      Boundable childBoundable = i.current as Boundable;
      if (childBoundable is AbstractNode) {
        size += sizeWithNode(childBoundable as AbstractNode);
      } else if (childBoundable is ItemBoundable) {
        size += 1;
      }
    }
    return size;
  }

  int depth() {
    if (isEmpty()) {
      return 0;
    }
    build();
    return depthWithNode(root);
  }

  int depthWithNode(AbstractNode node) {
    int maxChildDepth = 0;
    for (Iterator i = node.getChildBoundables().iterator; i.moveNext();) {
      Boundable childBoundable = i.current as Boundable;
      if (childBoundable is AbstractNode) {
        int childDepth = depthWithNode(childBoundable as AbstractNode);
        if (childDepth > maxChildDepth) maxChildDepth = childDepth;
      }
    }
    return maxChildDepth + 1;
  }

  void insertObj(Object bounds, Object item) {
    Assert.isTrue(!built,
        "Cannot insert items into an STR packed R-tree after it has been built.");
    itemBoundables!.add(new ItemBoundable(bounds, item));
  }

  /**
   *  Also builds the tree, if necessary.
   */
  List queryObj(Object searchBounds) {
    build();
    List matches = [];
    if (isEmpty()) {
      //Assert.isTrue(root.getBounds() == null);
      return matches;
    }
    if (getIntersectsOp().intersects(root.getBounds(), searchBounds)) {
      queryInternal(searchBounds, root, matches);
    }
    return matches;
  }

  /**
   *  Also builds the tree, if necessary.
   */
  void queryObjWithVisitor(Object searchBounds, ItemVisitor visitor) {
    build();
    if (isEmpty()) {
      // nothing in tree, so return
      //Assert.isTrue(root.getBounds() == null);
      return;
    }
    if (getIntersectsOp().intersects(root.getBounds(), searchBounds)) {
      queryInternalWithVisitor(searchBounds, root, visitor);
    }
  }

  /**
   * @return a test for intersection between two bounds, necessary because subclasses
   * of AbstractSTRtree have different implementations of bounds.
   * @see IntersectsOp
   */
  IntersectsOp getIntersectsOp();

  void queryInternal(Object searchBounds, AbstractNode node, List matches) {
    List childBoundables = node.getChildBoundables();
    for (int i = 0; i < childBoundables.length; i++) {
      Boundable childBoundable = childBoundables[i] as Boundable;
      if (!getIntersectsOp()
          .intersects(childBoundable.getBounds(), searchBounds)) {
        continue;
      }
      if (childBoundable is AbstractNode) {
        queryInternal(searchBounds, childBoundable as AbstractNode, matches);
      } else if (childBoundable is ItemBoundable) {
        matches.add((childBoundable as ItemBoundable).getItem());
      } else {
        Assert.shouldNeverReachHere();
      }
    }
  }

  void queryInternalWithVisitor(
      Object searchBounds, AbstractNode node, ItemVisitor visitor) {
    List childBoundables = node.getChildBoundables();
    for (int i = 0; i < childBoundables.length; i++) {
      Boundable childBoundable = childBoundables[i] as Boundable;
      if (!getIntersectsOp()
          .intersects(childBoundable.getBounds(), searchBounds)) {
        continue;
      }
      if (childBoundable is AbstractNode) {
        queryInternalWithVisitor(
            searchBounds, childBoundable as AbstractNode, visitor);
      } else if (childBoundable is ItemBoundable) {
        visitor.visitItem((childBoundable as ItemBoundable).getItem());
      } else {
        Assert.shouldNeverReachHere();
      }
    }
  }

  /**
   * Gets a tree structure (as a nested list)
   * corresponding to the structure of the items and nodes in this tree.
   * <p>
   * The returned {@link List}s contain either {@link Object} items,
   * or Lists which correspond to subtrees of the tree
   * Subtrees which do not contain any items are not included.
   * <p>
   * Builds the tree if necessary.
   *
   * @return a List of items and/or Lists
   */
  List itemsTree() {
    build();

    List? valuesTree = itemsTreeWithNode(root);
    if (valuesTree == null) return [];
    return valuesTree;
  }

  List? itemsTreeWithNode(AbstractNode node) {
    List valuesTreeForNode = [];
    for (Iterator i = node.getChildBoundables().iterator; i.moveNext();) {
      Boundable childBoundable = i.current as Boundable;
      if (childBoundable is AbstractNode) {
        List? valuesTreeForChild =
            itemsTreeWithNode(childBoundable as AbstractNode);
        // only add if not null (which indicates an item somewhere in this tree
        if (valuesTreeForChild != null)
          valuesTreeForNode.add(valuesTreeForChild);
      } else if (childBoundable is ItemBoundable) {
        valuesTreeForNode.add((childBoundable as ItemBoundable).getItem());
      } else {
        Assert.shouldNeverReachHere();
      }
    }
    if (valuesTreeForNode.length <= 0) return null;
    return valuesTreeForNode;
  }

  /**
   * Removes an item from the tree.
   * (Builds the tree, if necessary.)
   */
  bool removeObj(Object searchBounds, Object item) {
    build();
    if (getIntersectsOp().intersects(root.getBounds(), searchBounds)) {
      return removeWithNode(searchBounds, root, item);
    }
    return false;
  }

  bool removeItem(AbstractNode node, Object item) {
    Boundable? childToRemove = null;
    for (Iterator i = node.getChildBoundables().iterator; i.moveNext();) {
      Boundable childBoundable = i.current as Boundable;
      if (childBoundable is ItemBoundable) {
        if ((childBoundable as ItemBoundable).getItem() == item)
          childToRemove = childBoundable;
      }
    }
    if (childToRemove != null) {
      node.getChildBoundables().remove(childToRemove);
      return true;
    }
    return false;
  }

  bool removeWithNode(Object searchBounds, AbstractNode node, Object item) {
    // first try removing item from this node
    bool found = removeItem(node, item);
    if (found) return true;

    AbstractNode? childToPrune = null;
    // next try removing item from lower nodes
    for (Iterator i = node.getChildBoundables().iterator; i.moveNext();) {
      Boundable childBoundable = i.current as Boundable;
      if (!getIntersectsOp()
          .intersects(childBoundable.getBounds(), searchBounds)) {
        continue;
      }
      if (childBoundable is AbstractNode) {
        found =
            removeWithNode(searchBounds, childBoundable as AbstractNode, item);
        // if found, record child for pruning and exit
        if (found) {
          childToPrune = childBoundable as AbstractNode;
          break;
        }
      }
    }
    // prune child if possible
    if (childToPrune != null) {
      if (childToPrune.getChildBoundables().isEmpty) {
        node.getChildBoundables().remove(childToPrune);
      }
    }
    return found;
  }

  List boundablesAtLevel(int level) {
    List boundables = [];
    boundablesAtLevelWithNode(level, root, boundables);
    return boundables;
  }

  /**
   * @param level -1 to get items
   */
  void boundablesAtLevelWithNode(int level, AbstractNode top, List boundables) {
    Assert.isTrue(level > -2);
    if (top.getLevel() == level) {
      boundables.add(top);
      return;
    }
    for (Iterator i = top.getChildBoundables().iterator; i.moveNext();) {
      Boundable boundable = i.current as Boundable;
      if (boundable is AbstractNode) {
        boundablesAtLevelWithNode(level, boundable as AbstractNode, boundables);
      } else {
        Assert.isTrue(boundable is ItemBoundable);
        if (level == -1) {
          boundables.add(boundable);
        }
      }
    }
    return;
  }

  Comparator getComparator();
}

class STRtreeNode extends AbstractNode {
  STRtreeNode(int level) : super.withLevel(level);

  Object computeBounds() {
    Envelope? bounds = null;
    for (Iterator i = getChildBoundables().iterator; i.moveNext();) {
      Boundable childBoundable = i.current as Boundable;
      if (bounds == null) {
        bounds =
            new Envelope.fromEnvelope(childBoundable.getBounds() as Envelope);
      } else {
        bounds.expandToIncludeEnvelope(childBoundable.getBounds() as Envelope);
      }
    }
    return bounds!;
  }
}

/**
 *  A query-only R-tree created using the Sort-Tile-Recursive (STR) algorithm.
 *  For two-dimensional spatial data.
 * <P>
 *  The STR packed R-tree is simple to implement and maximizes space
 *  utilization; that is, as many leaves as possible are filled to capacity.
 *  Overlap between nodes is far less than in a basic R-tree. However, once the
 *  tree has been built (explicitly or on the first call to #query), items may
 *  not be added or removed.
 * <P>
 * Described in: P. Rigaux, Michel Scholl and Agnes Voisard.
 * <i>Spatial Databases With Application To GIS</i>.
 * Morgan Kaufmann, San Francisco, 2002.
 * <p>
 * <b>Note that inserting items into a tree is not thread-safe.</b>
 * Inserting performed on more than one thread must be synchronized externally.
 * <p>
 * Querying a tree is thread-safe.
 * The building phase is done synchronously,
 * and querying is stateless.
 *
 * @version 1.7
 */
class STRtree extends AbstractSTRtree implements SpatialIndex {
  /**
   *
   */

  static Comparator<dynamic> xComparator = (o1, o2) {
    return AbstractSTRtree.compareDoubles(
        centreX((o1 as Boundable).getBounds() as Envelope),
        centreX((o2 as Boundable).getBounds() as Envelope));
  };

//  new Comparator() {
//   int compare(Object o1, Object o2) {
//  return compareDoubles(
//  centreX((Envelope)((Boundable)o1).getBounds()),
//  centreX((Envelope)((Boundable)o2).getBounds()));
//  }
//  };
  static Comparator<dynamic> yComparator = (o1, o2) {
    return AbstractSTRtree.compareDoubles(
        centreY((o1 as Boundable).getBounds() as Envelope),
        centreY((o2 as Boundable).getBounds() as Envelope));
  };

//  new Comparator() {
//   int compare(Object o1, Object o2) {
//  return compareDoubles(
//  centreY((Envelope)((Boundable)o1).getBounds()),
//  centreY((Envelope)((Boundable)o2).getBounds()));
//  }
//  };

  static double centreX(Envelope e) {
    return avg(e.getMinX(), e.getMaxX());
  }

  static double centreY(Envelope e) {
    return avg(e.getMinY(), e.getMaxY());
  }

  static double avg(double a, double b) {
    return (a + b) / 2.0;
  }

  static IntersectsOp intersectsOp = EnvelopeIntersectsOp();

  /**
   * Creates the parent level for the given child level. First, orders the items
   * by the x-values of the midpoints, and groups them into vertical slices.
   * For each slice, orders the items by the y-values of the midpoints, and
   * group them into runs of size M (the node capacity). For each run, creates
   * a new (parent) node.
   */
  List createParentBoundables(List childBoundables, int newLevel) {
    Assert.isTrue(!childBoundables.isEmpty);
    int minLeafCount =
        ((childBoundables.length / getNodeCapacity().toDouble())).ceil();
    List sortedChildBoundables = List.from(childBoundables);
    sortedChildBoundables.sort(xComparator);
    List<List> verticalSlicesList =
        verticalSlices(sortedChildBoundables, (math.sqrt(minLeafCount)).ceil());
    return createParentBoundablesFromVerticalSlices(
        verticalSlicesList, newLevel);
  }

  List createParentBoundablesFromVerticalSlices(
      List<List> verticalSlices, int newLevel) {
    Assert.isTrue(verticalSlices.length > 0);
    List parentBoundables = [];
    for (int i = 0; i < verticalSlices.length; i++) {
      parentBoundables.addAll(
          createParentBoundablesFromVerticalSlice(verticalSlices[i], newLevel));
    }
    return parentBoundables;
  }

  List createParentBoundablesFromVerticalSlice(
      List childBoundables, int newLevel) {
    return super.createParentBoundables(childBoundables, newLevel);
  }

  /**
   * @param childBoundables Must be sorted by the x-value of the envelope midpoints
   */
  List<List> verticalSlices(List childBoundables, int sliceCount) {
    int sliceCapacity = (childBoundables.length / sliceCount.toDouble()).ceil();
    List<List> slices = []; //..length = sliceCount;
    Iterator i = childBoundables.iterator;
    for (int j = 0; j < sliceCount; j++) {
      slices.add([]);
      // slices[j] = [];
      int boundablesAddedToSlice = 0;
      while (boundablesAddedToSlice < sliceCapacity && i.moveNext()) {
        Boundable childBoundable = i.current as Boundable;
        slices[j].add(childBoundable);
        boundablesAddedToSlice++;
      }
    }
    return slices;
  }

  static final int DEFAULT_NODE_CAPACITY = 10;

  /**
   * Constructs an STRtree with the default node capacity.
   */
  STRtree() : this.withCapacity(DEFAULT_NODE_CAPACITY);

  /**
   * Constructs an STRtree with the given maximum number of child nodes that
   * a node may have.
   * <p>
   * The minimum recommended capacity setting is 4.
   *
   */
  STRtree.withCapacity(int nodeCapacity) : super.withCapacity(nodeCapacity);

  AbstractNode createNode(int level) {
    return new STRtreeNode(level);
  }

  IntersectsOp getIntersectsOp() {
    return intersectsOp;
  }

  /**
   * Inserts an item having the given bounds into the tree.
   */
  void insert(Envelope itemEnv, Object item) {
    if (itemEnv.isNull()) {
      return;
    }
    super.insertObj(itemEnv, item);
  }

  /**
   * Returns items whose bounds intersect the given envelope.
   */
  List query(Envelope searchEnv) {
    //Yes this method does something. It specifies that the bounds is an
    //Envelope. super.query takes an Object, not an Envelope. [Jon Aquino 10/24/2003]
    return super.queryObj(searchEnv);
  }

  /**
   * Returns items whose bounds intersect the given envelope.
   */
  void queryWithVisitor(Envelope searchEnv, ItemVisitor visitor) {
    //Yes this method does something. It specifies that the bounds is an
    //Envelope. super.query takes an Object, not an Envelope. [Jon Aquino 10/24/2003]
    super.queryObjWithVisitor(searchEnv, visitor);
  }

  /**
   * Removes a single item from the tree.
   *
   * @param itemEnv the Envelope of the item to remove
   * @param item the item to remove
   * @return <code>true</code> if the item was found
   */
  bool remove(Envelope itemEnv, Object item) {
    return super.removeObj(itemEnv, item);
  }

  /**
   * Returns the number of items in the tree.
   *
   * @return the number of items in the tree
   */
  int size() {
    return super.size();
  }

  /**
   * Returns the number of items in the tree.
   *
   * @return the number of items in the tree
   */
  int depth() {
    return super.depth();
  }

  Comparator getComparator() {
    return yComparator;
  }

  /**
   * Finds the two nearest items in the tree,
   * using {@link ItemDistance} as the distance metric.
   * A Branch-and-Bound tree traversal algorithm is used
   * to provide an efficient search.
   * <p>
   * If the tree is empty, the return value is <code>null</code.
   * If the tree contains only one item,
   * the return value is a pair containing that item.
   * <b>
   * If it is required to find only pairs of distinct items,
   * the {@link ItemDistance} function must be <b>anti-reflexive</b>.
   *
   * @param itemDist a distance metric applicable to the items in this tree
   * @return the pair of the nearest items
   *    or <code>null</code> if the tree is empty
   */
  List<Object>? nearestNeighbour(ItemDistance itemDist) {
    if (isEmpty()) return null;

    // if tree has only one item this will return null
    BoundablePair bp =
        new BoundablePair(this.getRoot(), this.getRoot(), itemDist);
    return nearestNeighbourWithPair(bp);
  }

  /**
   * Finds the item in this tree which is nearest to the given {@link Object},
   * using {@link ItemDistance} as the distance metric.
   * A Branch-and-Bound tree traversal algorithm is used
   * to provide an efficient search.
   * <p>
   * The query <tt>object</tt> does <b>not</b> have to be
   * contained in the tree, but it does
   * have to be compatible with the <tt>itemDist</tt>
   * distance metric.
   *
   * @param env the envelope of the query item
   * @param item the item to find the nearest neighbour of
   * @param itemDist a distance metric applicable to the items in this tree and the query item
   * @return the nearest item in this tree
   *    or <code>null</code> if the tree is empty
   */
  Object nearestNeighbourWithEnvelope(
      Envelope env, Object item, ItemDistance itemDist) {
    Boundable bnd = new ItemBoundable(env, item);
    BoundablePair bp = new BoundablePair(this.getRoot(), bnd, itemDist);
    return nearestNeighbourWithPair(bp)![0];
  }

  /**
   * Finds the two nearest items from this tree
   * and another tree,
   * using {@link ItemDistance} as the distance metric.
   * A Branch-and-Bound tree traversal algorithm is used
   * to provide an efficient search.
   * The result value is a pair of items,
   * the first from this tree and the second
   * from the argument tree.
   *
   * @param tree another tree
   * @param itemDist a distance metric applicable to the items in the trees
   * @return the pair of the nearest items, one from each tree
   *    or <code>null</code> if no pair of distinct items can be found
   */
  List<Object>? nearestNeighbourWithTree(STRtree tree, ItemDistance itemDist) {
    if (isEmpty() || tree.isEmpty()) return null;
    BoundablePair bp =
        new BoundablePair(this.getRoot(), tree.getRoot(), itemDist);
    return nearestNeighbourWithPair(bp);
  }

  List<Object>? nearestNeighbourWithPair(BoundablePair initBndPair) {
    double distanceLowerBound = double.infinity;
    BoundablePair? minPair = null;

    // initialize search queue
    PriorityQueue priQ = new PriorityQueue();
    priQ.add(initBndPair);

    while (!priQ.isEmpty() && distanceLowerBound > 0.0) {
      // pop head of queue and expand one side of pair
      BoundablePair bndPair = priQ.poll() as BoundablePair;
      double pairDistance = bndPair.getDistance();

      /**
       * If the distance for the first pair in the queue
       * is >= current minimum distance, other nodes
       * in the queue must also have a greater distance.
       * So the current minDistance must be the true minimum,
       * and we are done.
       */
      if (pairDistance >= distanceLowerBound) break;

      /**
       * If the pair members are leaves
       * then their distance is the exact lower bound.
       * Update the distanceLowerBound to reflect this
       * (which must be smaller, due to the test
       * immediately prior to this).
       */
      if (bndPair.isLeaves()) {
        // assert: currentDistance < minimumDistanceFound
        distanceLowerBound = pairDistance;
        minPair = bndPair;
      } else {
        /**
         * Otherwise, expand one side of the pair,
         * and insert the expanded pairs into the queue.
         * The choice of which side to expand is determined heuristically.
         */
        bndPair.expandToQueue(priQ, distanceLowerBound);
      }
    }
    if (minPair == null) return null;
    // done - return items with min distance
    return [
      (minPair.getBoundable(0) as ItemBoundable).getItem(),
      (minPair.getBoundable(1) as ItemBoundable).getItem()
    ];
  }

  /**
   * Tests whether some two items from this tree and another tree
   * lie within a given distance.
   * {@link ItemDistance} is used as the distance metric.
   * A Branch-and-Bound tree traversal algorithm is used
   * to provide an efficient search.
   *
   * @param tree another tree
   * @param itemDist a distance metric applicable to the items in the trees
   * @param maxDistance the distance limit for the search
   * @return true if there are items within the distance
   */
  bool isWithinDistanceWithTree(
      STRtree tree, ItemDistance itemDist, double maxDistance) {
    BoundablePair bp =
        new BoundablePair(this.getRoot(), tree.getRoot(), itemDist);
    return isWithinDistance(bp, maxDistance);
  }

  /**
   * Performs a withinDistance search on the tree node pairs.
   * This is a different search algorithm to nearest neighbour.
   * It can utilize the {@link BoundablePair#maximumDistance()} between
   * tree nodes to confirm if two internal nodes must
   * have items closer than the maxDistance,
   * and short-circuit the search.
   *
   * @param initBndPair the initial pair containing the tree root nodes
   * @param maxDistance the maximum distance to search for
   * @return true if two items lie within the given distance
   */
  bool isWithinDistance(BoundablePair initBndPair, double maxDistance) {
    double distanceUpperBound = double.infinity;

    // initialize search queue
    PriorityQueue priQ = new PriorityQueue();
    priQ.add(initBndPair);

    while (!priQ.isEmpty()) {
      // pop head of queue and expand one side of pair
      BoundablePair bndPair = priQ.poll() as BoundablePair;
      double pairDistance = bndPair.getDistance();

      /**
       * If the distance for the first pair in the queue
       * is > maxDistance, all other pairs
       * in the queue must have a greater distance as well.
       * So can conclude no items are within the distance
       * and terminate with result = false
       */
      if (pairDistance > maxDistance) return false;

      /**
       * If the maximum distance between the nodes
       * is less than the maxDistance,
       * than all items in the nodes must be
       * closer than the max distance.
       * Then can terminate with result = true.
       *
       * NOTE: using Envelope MinMaxDistance
       * would provide a tighter bound,
       * but not much performance improvement has been observed
       */
      if (bndPair.maximumDistance() <= maxDistance) return true;
      /**
       * If the pair items are leaves
       * then their actual distance is an upper bound.
       * Update the distanceUpperBound to reflect this
       */
      if (bndPair.isLeaves()) {
        // assert: currentDistance < minimumDistanceFound
        distanceUpperBound = pairDistance;

        /**
         * If the items are closer than maxDistance
         * can terminate with result = true.
         */
        if (distanceUpperBound <= maxDistance) return true;
      } else {
        /**
         * Otherwise, expand one side of the pair,
         * and insert the expanded pairs into the queue.
         * The choice of which side to expand is determined heuristically.
         */
        bndPair.expandToQueue(priQ, distanceUpperBound);
      }
    }
    return false;
  }

  /**
   * Finds k items in this tree which are the top k nearest neighbors to the given {@code item},
   * using {@code itemDist} as the distance metric.
   * A Branch-and-Bound tree traversal algorithm is used
   * to provide an efficient search.
   * This method implements the KNN algorithm described in the following paper:
   * <p>
   * Roussopoulos, Nick, Stephen Kelley, and Frédéric Vincent. "Nearest neighbor queries."
   * ACM sigmod record. Vol. 24. No. 2. ACM, 1995.
   * <p>
   * The query {@code item} does <b>not</b> have to be
   * contained in the tree, but it does
   * have to be compatible with the {@code itemDist}
   * distance metric.
   *
   * @param env the envelope of the query item
   * @param item the item to find the nearest neighbour of
   * @param itemDist a distance metric applicable to the items in this tree and the query item
   * @param k the K nearest items in kNearestNeighbour
   * @return the K nearest items in this tree
   */
  List<Object> nearestNeighbourWithEnv(
      Envelope env, Object item, ItemDistance itemDist, int k) {
    Boundable bnd = new ItemBoundable(env, item);
    BoundablePair bp = new BoundablePair(this.getRoot(), bnd, itemDist);
    return nearestNeighbourK(bp, k);
  }

  List<Object> nearestNeighbourK(BoundablePair initBndPair, int k) {
    return nearestNeighbourKWithMaxD(initBndPair, double.infinity, k);
  }

  List<Object> nearestNeighbourKWithMaxD(
      BoundablePair initBndPair, double maxDistance, int k) {
    double distanceLowerBound = maxDistance;

    // initialize internal structures
    PriorityQueue priQ = new PriorityQueue();

    // initialize queue
    priQ.add(initBndPair);

    PriorityQueue kNearestNeighbors = new PriorityQueue();

    while (!priQ.isEmpty() && distanceLowerBound >= 0.0) {
      // pop head of queue and expand one side of pair
      BoundablePair bndPair = priQ.poll() as BoundablePair;
      double pairDistance = bndPair.getDistance();

      /**
       * If the distance for the first node in the queue
       * is >= the current maximum distance in the k queue , all other nodes
       * in the queue must also have a greater distance.
       * So the current minDistance must be the true minimum,
       * and we are done.
       */
      if (pairDistance >= distanceLowerBound) {
        break;
      }
      /**
       * If the pair members are leaves
       * then their distance is the exact lower bound.
       * Update the distanceLowerBound to reflect this
       * (which must be smaller, due to the test
       * immediately prior to this).
       */
      if (bndPair.isLeaves()) {
        // assert: currentDistance < minimumDistanceFound

        if (kNearestNeighbors.size() < k) {
          kNearestNeighbors.add(bndPair);
        } else {
          BoundablePair bp1 = kNearestNeighbors.peek() as BoundablePair;
          if (bp1.getDistance() > pairDistance) {
            kNearestNeighbors.poll();
            kNearestNeighbors.add(bndPair);
          }
          /*
    		   * minDistance should be the farthest point in the K nearest neighbor queue.
    		   */
          BoundablePair bp2 = kNearestNeighbors.peek() as BoundablePair;
          distanceLowerBound = bp2.getDistance();
        }
      } else {
        /**
         * Otherwise, expand one side of the pair,
         * (the choice of which side to expand is heuristically determined)
         * and insert the new expanded pairs into the queue
         */
        bndPair.expandToQueue(priQ, distanceLowerBound);
      }
    }
    // done - return items with min distance

    return getItems(kNearestNeighbors);
  }

  static List<Object> getItems(PriorityQueue kNearestNeighbors) {
    /**
     * Iterate the K Nearest Neighbour Queue and retrieve the item from each BoundablePair
     * in this queue
     */
    List<Object> items = []; //..length = kNearestNeighbors.size();
    // int count = 0;
    while (!kNearestNeighbors.isEmpty()) {
      BoundablePair bp = kNearestNeighbors.poll() as BoundablePair;
      items.add((bp.getBoundable(0) as ItemBoundable).getItem());
      // items[count] = (bp.getBoundable(0) as ItemBoundable).getItem();
      // count++;
    }
    return items;
  }
}
