part of dart_jts;

/**
 * An implementation of a
 * <a href='https://en.wikipedia.org/wiki/K-d_tree'>KD-Tree</a>
 * over two dimensions (X and Y).
 * KD-trees provide fast range searching and fast lookup for point data.
 * The tree is built dynamically by inserting points.
 * The tree supports queries by range and for point equality.
 * For querying an internal stack is used instead of recursion to avoid overflow.
 * <p>
 * This implementation supports detecting and snapping points which are closer
 * than a given distance tolerance.
 * If the same point (up to tolerance) is inserted
 * more than once, it is snapped to the existing node.
 * In other words, if a point is inserted which lies
 * within the tolerance of a node already in the index,
 * it is snapped to that node.
 * When an inserted point is snapped to a node then a new node is not created
 * but the count of the existing node is incremented.
 * If more than one node in the tree is within tolerance of an inserted point,
 * the closest and then lowest node is snapped to.
 * <p>
 * The structure of a KD-Tree depends on the order of insertion of the points.
 * A tree may become unbalanced if the inserted points are coherent
 * (e.g. monotonic in one or both dimensions).
 * A perfectly balanced tree has depth of only log2(N),
 * but an unbalanced tree may be much deeper.
 * This has a serious impact on query efficiency.
 * One solution to this is to randomize the order of points before insertion
 * (e.g. by using <a href="https://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle">Fisher-Yates shuffling</a>).
 *
 * @author David Skea
 * @author Martin Davis
 */
class KdTree {
  /**
   * Converts a collection of {@link KdNode}s to an array of {@link Coordinate}s.
   *
   * @param kdnodes
   *          a collection of nodes
   * @return an array of the coordinates represented by the nodes
   */
  static List<Coordinate> toCoordinates(List<KdNode> kdnodes) {
    return toCoordinatesRepeated(kdnodes, false);
  }

  /**
   * Converts a collection of {@link KdNode}s
   * to an array of {@link Coordinate}s,
   * specifying whether repeated nodes should be represented
   * by multiple coordinates.
   *
   * @param kdnodes a collection of nodes
   * @param includeRepeated true if repeated nodes should
   *   be included multiple times
   * @return an array of the coordinates represented by the nodes
   */
  static List<Coordinate> toCoordinatesRepeated(
      List<KdNode> kdnodes, bool includeRepeated) {
    CoordinateList coord = new CoordinateList();
    for (var node in kdnodes) {
      int count = includeRepeated ? node.getCount() : 1;
      for (int i = 0; i < count; i++) {
        coord.addCoord(node.getCoordinate()!, true);
      }
    }
    return coord.toCoordinateArray(true);
  }

  KdNode? root = null;
  int numberOfNodes = 0;
  double tolerance = 0.0;

  /**
   * Creates a new instance of a KdTree with a snapping tolerance of 0.0. (I.e.
   * distinct points will <i>not</i> be snapped)
   */
  factory KdTree.withZeroTolerance() {
    return KdTree(0.0);
  }

  /**
   * Creates a new instance of a KdTree, specifying a snapping distance
   * tolerance. Points which lie closer than the tolerance to a point already in
   * the tree will be treated as identical to the existing point.
   *
   * @param tolerance
   *          the tolerance distance for considering two points equal
   */
  KdTree(this.tolerance);

  /**
   * Gets the root node of this tree.
   *
   * @return the root node of the tree
   */
  KdNode? getRoot() {
    return root;
  }

  /**
   * Tests whether the index contains any items.
   *
   * @return true if the index does not contain any items
   */
  bool isEmpty() {
    if (root == null) return true;
    return false;
  }

  /**
   * Inserts a new point in the kd-tree, with no data.
   *
   * @param p
   *          the point to insert
   * @return the kdnode containing the point
   */
  KdNode insertCoordinate(Coordinate p) {
    return insert(p, null);
  }

  /**
   * Inserts a new point into the kd-tree.
   *
   * @param p
   *          the point to insert
   * @param data
   *          a data item for the point
   * @return returns a new KdNode if a new point is inserted, else an existing
   *         node is returned with its counter incremented. This can be checked
   *         by testing returnedNode.getCount() &gt; 1.
   */
  KdNode insert(Coordinate p, Object? data) {
    if (root == null) {
      root = new KdNode.fromCoordinate(p, data);
      return root!;
    }

    /**
     * Check if the point is already in the tree, up to tolerance.
     * If tolerance is zero, this phase of the insertion can be skipped.
     */
    if (tolerance > 0) {
      KdNode? matchNode = _findBestMatchNode(p);
      if (matchNode != null) {
        // point already in index - increment counter
        matchNode.increment();
        return matchNode;
      }
    }

    return insertExact(p, data);
  }

  /**
   * Finds the node in the tree which is the best match for a point
   * being inserted.
   * The match is made deterministic by returning the lowest of any nodes which
   * lie the same distance from the point.
   * There may be no match if the point is not within the distance tolerance of any
   * existing node.
   *
   * @param p the point being inserted
   * @return the best matching node
   * @return null if no match was found
   */
  KdNode? _findBestMatchNode(Coordinate p) {
    BestMatchVisitor visitor = new BestMatchVisitor(p, tolerance);
    query(visitor.queryEnvelope(), visitor);
    return visitor.getNode();
  }

  /**
   * Inserts a point known to be beyond the distance tolerance of any existing node.
   * The point is inserted at the bottom of the exact splitting path,
   * so that tree shape is deterministic.
   *
   * @param p the point to insert
   * @param data the data for the point
   * @return the created node
   */
  KdNode insertExact(Coordinate p, Object? data) {
    KdNode? currentNode = root;
    KdNode? leafNode = root;
    bool isXLevel = true;
    bool isLessThan = true;

    /**
   * Traverse the tree, first cutting the plane left-right (by X ordinate)
   * then top-bottom (by Y ordinate)
   */
    while (currentNode != null) {
      bool isInTolerance =
          p.distance(currentNode.getCoordinate()!) <= tolerance;

      // check if point is already in tree (up to tolerance) and if so simply
      // return existing node
      if (isInTolerance) {
        currentNode.increment();
        return currentNode;
      }

      double splitValue = currentNode.splitValue(isXLevel);
      if (isXLevel) {
        isLessThan = p.x < splitValue;
      } else {
        isLessThan = p.y < splitValue;
      }
      leafNode = currentNode;
      if (isLessThan) {
        //System.out.print("L");
        currentNode = currentNode.getLeft();
      } else {
        //System.out.print("R");
        currentNode = currentNode.getRight();
      }

      isXLevel = !isXLevel;
    }
    //System.out.println("<<");
    // no node found, add new leaf node to tree
    numberOfNodes = numberOfNodes + 1;
    KdNode node = new KdNode.fromCoordinate(p, data);
    if (isLessThan) {
      leafNode?.setLeft(node);
    } else {
      leafNode?.setRight(node);
    }
    return node;
  }

  /**
   * Performs a range search of the points in the index and visits all nodes found.
   *
   * @param queryEnv the range rectangle to query
   * @param visitor a visitor to visit all nodes found by the search
   */
  void query(Envelope queryEnv, KdNodeVisitor visitor) {
    //-- Deque is faster than Stack
    Queue<QueryStackFrame> queryStack = Queue();
    KdNode? currentNode = root;
    bool isXLevel = true;

    // search is computed via in-order traversal
    while (true) {
      if (currentNode != null) {
        queryStack.addLast(QueryStackFrame(currentNode, isXLevel));

        bool searchLeft = currentNode.isRangeOverLeft(isXLevel, queryEnv);
        if (searchLeft) {
          currentNode = currentNode.getLeft();
          if (currentNode != null) {
            isXLevel = !isXLevel;
          }
        } else {
          currentNode = null;
        }
      } else if (!queryStack.isEmpty) {
        // currentNode is empty, so pop stack
        QueryStackFrame frame = queryStack.removeFirst();
        currentNode = frame.node;
        isXLevel = frame.isXLevel;

        //-- check if search matches current node
        if (queryEnv.containsCoordinate(currentNode.getCoordinate()!)) {
          visitor.visit(currentNode);
        }

        bool searchRight = currentNode.isRangeOverRight(isXLevel, queryEnv);
        if (searchRight) {
          currentNode = currentNode.getRight();
          if (currentNode != null) {
            isXLevel = !isXLevel;
          }
        } else {
          currentNode = null;
        }
      } else {
        //-- stack is empty and no current node
        return;
      }
    }
  }

  /**
   * Performs a range search of the points in the index.
   *
   * @param queryEnv the range rectangle to query
   * @return a list of the KdNodes found
   */
  List<KdNode> queryEnvelope(Envelope queryEnv) {
    final List<KdNode> result = List.empty(growable: true);
    query(queryEnv, KdNodeVisitorEnvelope(result));
    return result;
  }

  /**
   * Performs a range search of the points in the index.
   *
   * @param queryEnv
   *          the range rectangle to query
   * @param result
   *          a list to accumulate the result nodes into
   */
  void queryList(Envelope queryEnv, final List<KdNode> result) {
    query(queryEnv, KdNodeVisitorEnvelope(result));
  }

  /**
   * Searches for a given point in the index and returns its node if found.
   *
   * @param queryPt the point to query
   * @return the point node, if it is found in the index, or null if not
   */
  KdNode? queryCoordinate(Coordinate queryPt) {
    KdNode? currentNode = root;
    bool isXLevel = true;

    while (currentNode != null) {
      if (currentNode.getCoordinate()!.equals2D(queryPt)) return currentNode;

      bool searchLeft = currentNode.isPointOnLeft(isXLevel, queryPt);
      if (searchLeft) {
        currentNode = currentNode.getLeft();
      } else {
        currentNode = currentNode.getRight();
      }
      isXLevel = !isXLevel;
    }
    //-- point not found
    return null;
  }

  /**
   * Computes the depth of the tree.
   *
   * @return the depth of the tree
   */
  int depth() {
    return _depthNode(root!);
  }

  int _depthNode(KdNode? currentNode) {
    if (currentNode == null) return 0;

    int dL = _depthNode(currentNode.getLeft());
    int dR = _depthNode(currentNode.getRight());
    return 1 + (dL > dR ? dL : dR);
  }

  /**
   * Computes the size (number of items) in the tree.
   *
   * @return the size of the tree
   */
  int size() {
    if (root != null) {
      return _sizeNode(root!);
    } else {
      return 0;
    }
  }

  int _sizeNode(KdNode? currentNode) {
    if (currentNode == null) return 0;

    int sizeL = _sizeNode(currentNode.getLeft()!);
    int sizeR = _sizeNode(currentNode.getRight()!);
    return 1 + sizeL + sizeR;
  }
}

/**
 * A node of a {@link KdTree}, which represents one or more points in the same location.
 *
 * @author dskea
 */
class KdNode {
  Coordinate? p = null;
  Object? data;
  KdNode? left;
  KdNode? right;
  late int count;

  /**
   * Creates a new KdNode.
   *
   * @param _x coordinate of point
   * @param _y coordinate of point
   * @param data a data objects to associate with this node
   */
  KdNode(double _x, double _y, this.data) {
    p = new Coordinate(_x, _y);
    left = null;
    right = null;
    count = 1;
  }

  /**
   * Creates a new KdNode.
   *
   * @param p point location of new node
   * @param data a data objects to associate with this node
   */
  KdNode.fromCoordinate(Coordinate p, this.data) {
    this.p = Coordinate(p.x, p.y);
    left = null;
    right = null;
    count = 1;
    this.data = data;
  }

  /**
   * Returns the X coordinate of the node
   *
   * @return X coordinate of the node
   */
  double getX() {
    return p!.x;
  }

  /**
   * Returns the Y coordinate of the node
   *
   * @return Y coordinate of the node
   */
  double getY() {
    return p!.y;
  }

  /**
   * Gets the split value at a node, depending on
   * whether the node splits on X or Y.
   * The X (or Y) ordinates of all points in the left subtree
   * are less than the split value, and those
   * in the right subtree are greater than or equal to the split value.
   *
   * @param isSplitOnX whether the node splits on X or Y
   * @return the splitting value
   */
  double splitValue(bool isSplitOnX) {
    if (isSplitOnX) {
      return p!.getX();
    }
    return p!.getY();
  }

  /**
   * Returns the location of this node
   *
   * @return p location of this node
   */
  Coordinate? getCoordinate() {
    return p;
  }

  /**
   * Gets the user data object associated with this node.
   * @return user data
   */
  Object? getData() {
    return data;
  }

  /**
   * Returns the left node of the tree
   *
   * @return left node
   */
  KdNode? getLeft() {
    return left;
  }

  /**
   * Returns the right node of the tree
   *
   * @return right node
   */
  KdNode? getRight() {
    return right;
  }

  // Increments counts of points at this location
  void increment() {
    count = count + 1;
  }

  /**
   * Returns the number of inserted points that are coincident at this location.
   *
   * @return number of inserted points that this node represents
   */
  int getCount() {
    return count;
  }

  /**
   * Tests whether more than one point with this value have been inserted (up to the tolerance)
   *
   * @return true if more than one point have been inserted with this value
   */
  bool isRepeated() {
    return count > 1;
  }

  // Sets left node value
  void setLeft(KdNode _left) {
    left = _left;
  }

  // Sets right node value
  void setRight(KdNode _right) {
    right = _right;
  }

  /**
   * Tests whether the node's left subtree may contain values
   * in a given range envelope.
   *
   * @param isSplitOnX whether the node splits on  X or Y
   * @param env the range envelope
   * @return true if the left subtree is in range
   */
  bool isRangeOverLeft(bool isSplitOnX, Envelope env) {
    double envMin;
    if (isSplitOnX) {
      envMin = env.getMinX();
    } else {
      envMin = env.getMinY();
    }
    double splitVal = splitValue(isSplitOnX);
    bool isInRange = envMin < splitVal;
    return isInRange;
  }

  /**
   * Tests whether the node's right subtree may contain values
   * in a given range envelope.
   *
   * @param isSplitOnX whether the node splits on  X or Y
   * @param env the range envelope
   * @return true if the right subtree is in range
   */
  bool isRangeOverRight(bool isSplitOnX, Envelope env) {
    double envMax;
    if (isSplitOnX) {
      envMax = env.getMaxX();
    } else {
      envMax = env.getMaxY();
    }
    double splitVal = splitValue(isSplitOnX);
    bool isInRange = splitVal <= envMax;
    return isInRange;
  }

  /**
   * Tests whether a point is strictly to the left
   * of the splitting plane for this node.
   * If so it may be in the left subtree of this node,
   * Otherwise, the point may be in the right subtree.
   * The point is to the left if its X (or Y) ordinate
   * is less than the split value.
   *
   * @param isSplitOnX whether the node splits on  X or Y
   * @param pt the query point
   * @return true if the point is strictly to the left.
   *
   * @see #splitValue(bool)
   */
  bool isPointOnLeft(bool isSplitOnX, Coordinate pt) {
    double ptOrdinate;
    if (isSplitOnX) {
      ptOrdinate = pt.x;
    } else {
      ptOrdinate = pt.y;
    }
    double splitVal = splitValue(isSplitOnX);
    bool isInRange = (ptOrdinate < splitVal);
    return isInRange;
  }
}

class BestMatchVisitor implements KdNodeVisitor {
  double tolerance;
  KdNode? matchNode = null;
  double matchDist = 0.0;
  Coordinate p;

  BestMatchVisitor(this.p, this.tolerance);

  Envelope queryEnvelope() {
    Envelope queryEnv = Envelope.fromCoordinate(p);
    queryEnv.expandByDistance(tolerance);
    return queryEnv;
  }

  KdNode? getNode() {
    return matchNode;
  }

  void visit(KdNode node) {
    double dist = p.distance(node.getCoordinate()!);
    bool isInTolerance = dist <= tolerance;
    if (!isInTolerance) return;
    bool update = false;
    if (matchNode == null ||
        dist < matchDist
        // if distances are the same, record the lesser coordinate
        ||
        (matchNode != null &&
            dist == matchDist &&
            node.getCoordinate()!.compareTo((matchNode!.getCoordinate())!) <
                1)) {
      update = true;
    }

    if (update) {
      matchNode = node;
      matchDist = dist;
    }
  }
}
/**
 * A visitor for {@link KdNode}s in a {@link KdTree} index.
 *
 * @version 1.7
 */

abstract class KdNodeVisitor {
/**
 * Visits a node.
 *
 * @param node the node to visit
 */
  void visit(KdNode node);
}

class KdNodeVisitorEnvelope extends KdNodeVisitor {
  List<KdNode> result;
  KdNodeVisitorEnvelope(this.result);
  @override
  void visit(KdNode node) {
    result.add(node);
  }
}

class QueryStackFrame {
  KdNode node;
  bool isXLevel = false;

  QueryStackFrame(this.node, this.isXLevel);
}
