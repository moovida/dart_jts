part of dart_jts;

class Quadtree implements SpatialIndex {
  late Root root;
  double minExtent = 1.0;
  Quadtree() {
    root = Root();
  }
  static Envelope ensureExtent(Envelope itemEnv, double minExtent) {
    //The names "ensureExtent" and "minExtent" are misleading -- sounds like
    //this method ensures that the extents are greater than minExtent.
    //Perhaps we should rename them to "ensurePositiveExtent" and "defaultExtent".
    //[Jon Aquino]
    double minx = itemEnv.getMinX();
    double maxx = itemEnv.getMaxX();
    double miny = itemEnv.getMinY();
    double maxy = itemEnv.getMaxY();
    // has a non-zero extent
    if (minx != maxx && miny != maxy) return itemEnv;

    // pad one or both extents
    if (minx == maxx) {
      minx = minx - minExtent / 2.0;
      maxx = maxx + minExtent / 2.0;
    }
    if (miny == maxy) {
      miny = miny - minExtent / 2.0;
      maxy = maxy + minExtent / 2.0;
    }
    return new Envelope(minx, maxx, miny, maxy);
  }

  ///
  /// Returns the number of levels in the tree.
  ///
  int depth() {
    return root.depth();
  }

  ///
  /// Tests whether the index contains any items.
  ///
  /// @return true if the index does not contain any items
  ///
  bool isEmpty() {
    return root.isEmpty();
  }

  ///
  /// Returns the number of items in the tree.
  ///
  /// @return the number of items in the tree
  ///
  int size() {
    return root.size();
  }

  @override
  void insert(Envelope itemEnv, Object item) {
    collectStats(itemEnv);
    Envelope insertEnv = ensureExtent(itemEnv, minExtent);
    root.insert(insertEnv, item);
  }

  void collectStats(Envelope itemEnv) {
    double delX = itemEnv.getWidth();
    if (delX < minExtent && delX > 0.0) minExtent = delX;

    double delY = itemEnv.getHeight();
    if (delY < minExtent && delY > 0.0) minExtent = delY;
  }

  ///
  /// Queries the tree and returns items which may lie in the given search envelope.
  /// Precisely, the items that are returned are all items in the tree
  /// whose envelope <b>may</b> intersect the search Envelope.
  /// Note that some items with non-intersecting envelopes may be returned as well;
  /// the client is responsible for filtering these out.
  /// In most situations there will be many items in the tree which do not
  /// intersect the search envelope and which are not returned - thus
  /// providing improved performance over a simple linear scan.
  ///
  /// @param searchEnv the envelope of the desired query area.
  /// @return a List of items which may intersect the search envelope
  ///
  @override
  List<dynamic> query(Envelope? searchEnv) {
    ///
    /// the items that are matched are the items in quads which
    /// overlap the search envelope
    ///
    ArrayListVisitor visitor = new ArrayListVisitor();
    queryWithVisitor(searchEnv, visitor);
    return visitor.getItems();
  }

  ///
  /// Queries the tree and visits items which may lie in the given search envelope.
  /// Precisely, the items that are visited are all items in the tree
  /// whose envelope <b>may</b> intersect the search Envelope.
  /// Note that some items with non-intersecting envelopes may be visited as well;
  /// the client is responsible for filtering these out.
  /// In most situations there will be many items in the tree which do not
  /// intersect the search envelope and which are not visited - thus
  /// providing improved performance over a simple linear scan.
  ///
  /// @param searchEnv the envelope of the desired query area.
  /// @param visitor a visitor object which is passed the visited items
  ///

  @override
  void queryWithVisitor(Envelope? searchEnv, ItemVisitor visitor) {
    ///
    /// the items that are matched are the items in quads which
    /// overlap the search envelope
    ///
    root.visit(searchEnv, visitor);
  }

  ////
  /// Return a list of all items in the Quadtree
  ///
  List<dynamic> queryAll() {
    List<dynamic> foundItems = List.empty(growable: true);
    root.addAllItems(foundItems);
    return foundItems;
  }

  ///
  /// Removes a single item from the tree.
  ///
  /// @param itemEnv the Envelope of the item to be removed
  /// @param item the item to remove
  /// @return <code>true</code> if the item was found (and thus removed)
  ///

  @override
  bool remove(Envelope itemEnv, Object item) {
    Envelope posEnv = ensureExtent(itemEnv, minExtent);
    return root.remove(posEnv, item);
  }
}

abstract class NodeBase {
  ///
  /// Gets the index of the subquad that wholly contains the given envelope.
  /// If none does, returns -1.
  ///
  /// @return the index of the subquad that wholly contains the given envelope
  /// or -1 if no subquad wholly contains the envelope
  ///
  static int getSubnodeIndex(Envelope? env, double centrex, double centrey) {
    int subnodeIndex = -1;
    if (env!.getMinX() >= centrex) {
      if (env.getMinY() >= centrey) subnodeIndex = 3;
      if (env.getMaxY() <= centrey) subnodeIndex = 1;
    }
    if (env.getMaxX() <= centrex) {
      if (env.getMinY() >= centrey) subnodeIndex = 2;
      if (env.getMaxY() <= centrey) subnodeIndex = 0;
    }
    return subnodeIndex;
  }

  List<dynamic> items = List.empty(growable: true);

  ///
  /// subquads are numbered as follows:
  /// <pre>
  ///  2 | 3
  ///  --+--
  ///  0 | 1
  /// </pre>
  ///
  List<NodeNode?> subnode = List.filled(4, null, growable: false);

  List<dynamic> getItems() {
    return items;
  }

  bool hasItems() {
    return !items.isEmpty;
  }

  void add(Object item) {
    items.add(item);
  }

  ///
  /// Removes a single item from this subtree.
  ///
  /// @param itemEnv the envelope containing the item
  /// @param item the item to remove
  /// @return <code>true</code> if the item was found and removed
  ///
  bool remove(Envelope itemEnv, Object item) {
    // use envelope to restrict nodes scanned
    if (!isSearchMatch(itemEnv)) return false;

    bool found = false;
    for (int i = 0; i < 4; i++) {
      if (subnode[i] != null) {
        found = subnode[i]!.remove(itemEnv, item);
        if (found) {
          // trim subtree if empty
          if (subnode[i]!.isPrunable()) subnode[i] = null;
          break;
        }
      }
    }
    // if item was found lower down, don't need to search for it here
    if (found) return found;
    // otherwise, try and remove the item from the list of items in this node
    found = items.remove(item);
    return found;
  }

  bool isPrunable() {
    return !(hasChildren() || hasItems());
  }

  bool hasChildren() {
    for (int i = 0; i < 4; i++) {
      if (subnode[i] != null) return true;
    }
    return false;
  }

  bool isEmpty() {
    bool isEmpty = true;
    if (!items.isEmpty)
      isEmpty = false;
    else {
      for (int i = 0; i < 4; i++) {
        if (subnode[i] != null) {
          if (!subnode[i]!.isEmpty()) {
            isEmpty = false;
            break;
          }
        }
      }
    }
    return isEmpty;
  }

  //<<TODO:RENAME?>> Sounds like this method adds resultItems to items
  //(like List#addAll). Perhaps it should be renamed to "addAllItemsTo" [Jon Aquino]
  List<dynamic> addAllItems(List<dynamic> resultItems) {
    // this node may have items as well as subnodes (since items may not
    // be wholely contained in any single subnode
    resultItems.addAll(this.items);
    for (int i = 0; i < 4; i++) {
      if (subnode[i] != null) {
        subnode[i]!.addAllItems(resultItems);
      }
    }
    return resultItems;
  }

  bool isSearchMatch(Envelope? searchEnv);

  void addAllItemsFromOverlapping(Envelope searchEnv, List resultItems) {
    if (!isSearchMatch(searchEnv)) return;

    // this node may have items as well as subnodes (since items may not
    // be wholely contained in any single subnode
    resultItems.addAll(items);

    for (int i = 0; i < 4; i++) {
      if (subnode[i] != null) {
        subnode[i]!.addAllItemsFromOverlapping(searchEnv, resultItems);
      }
    }
  }

  void visit(Envelope? searchEnv, ItemVisitor visitor) {
    if (searchEnv == null) {}
    if (!isSearchMatch(searchEnv)) return;

    // this node may have items as well as subnodes (since items may not
    // be wholely contained in any single subnode
    visitItems(searchEnv!, visitor);

    for (int i = 0; i < 4; i++) {
      if (subnode[i] != null) {
        subnode[i]!.visit(searchEnv, visitor);
      }
    }
  }

  void visitItems(Envelope searchEnv, ItemVisitor visitor) {
    // would be nice to filter items based on search envelope, but can't until they contain an envelope
    for (int i = 0; i < items.length; i++) {
      visitor.visitItem(items[i]);
    }
  }

//<<TODO:RENAME?>> In Samet's terminology, I think what we're returning here is
//actually level+1 rather than depth. (See p. 4 of his book) [Jon Aquino]
  int depth() {
    int maxSubDepth = 0;
    for (int i = 0; i < 4; i++) {
      if (subnode[i] != null) {
        int sqd = subnode[i]!.depth();
        if (sqd > maxSubDepth) maxSubDepth = sqd;
      }
    }
    return maxSubDepth + 1;
  }

  int size() {
    int subSize = 0;
    for (int i = 0; i < 4; i++) {
      if (subnode[i] != null) {
        subSize += subnode[i]!.size();
      }
    }
    return subSize + items.length;
  }

  int getNodeCount() {
    int subSize = 0;
    for (int i = 0; i < 4; i++) {
      if (subnode[i] != null) {
        subSize += subnode[i]!.size();
      }
    }
    return subSize + 1;
  }
}

class NodeNode extends NodeBase {
  static NodeNode createNode(Envelope env) {
    Key key = Key.fromEnvelope(env);
    NodeNode node = NodeNode(key.getEnvelope(), key.getLevel());
    return node;
  }

  static NodeNode createExpanded(NodeNode? node, Envelope addEnv) {
    Envelope expandEnv = Envelope.fromEnvelope(addEnv);
    if (node != null) expandEnv.expandToIncludeEnvelope(node.env!);

    NodeNode largerNode = createNode(expandEnv);
    if (node != null) largerNode.insertNode(node);
    return largerNode;
  }

  Envelope? env;
  late double centrex;
  late double centrey;
  int level;

  NodeNode(this.env, this.level) {
    //this.parent = parent;
    this.env = env;
    this.level = level;
    centrex = (env!.getMinX() + env!.getMaxX()) / 2;
    centrey = (env!.getMinY() + env!.getMaxY()) / 2;
  }

  Envelope? getEnvelope() {
    return env;
  }

  bool isSearchMatch(Envelope? searchEnv) {
    if (searchEnv == null) return false;
    return env!.intersectsEnvelope(searchEnv);
  }

  ///
  /// Returns the subquad containing the envelope <tt>searchEnv</tt>.
  /// Creates the subquad if
  /// it does not already exist.
  ///
  /// @return the subquad containing the search envelope
  ///
  NodeNode getNode(Envelope searchEnv) {
    int subnodeIndex = NodeBase.getSubnodeIndex(searchEnv, centrex, centrey);
    // if subquadIndex is -1 searchEnv is not contained in a subquad
    if (subnodeIndex != -1) {
      // create the quad if it does not exist
      NodeNode? node = getSubnode(subnodeIndex);
      // recursively search the found/created quad
      return node!.getNode(searchEnv);
    } else {
      return this;
    }
  }

  ///
  /// Returns the smallest <i>existing</i>
  /// node containing the envelope.
  ///
  NodeBase find(Envelope searchEnv) {
    int subnodeIndex = NodeBase.getSubnodeIndex(searchEnv, centrex, centrey);
    if (subnodeIndex == -1) return this;
    if (subnode[subnodeIndex] != null) {
      // query lies in subquad, so search it
      NodeNode? node = subnode[subnodeIndex];
      return node!.find(searchEnv);
    }
    // no existing subquad, so return this one anyway
    return this;
  }

  void insertNode(NodeNode node) {
    Assert.isTrue(env == null || env!.containsEnvelope(node.env!));
    int index = NodeBase.getSubnodeIndex(node.env, centrex, centrey);
    if (index == -1) {
      return;
    }
    if (node.level == level - 1) {
      subnode[index] = node;
//System.out.println("inserted");
    } else {
      // the quad is not a direct child, so make a new child quad to contain it
      // and recursively insert the quad
      NodeNode childNode = createSubnode(index);
      childNode.insertNode(node);
      subnode[index] = childNode;
    }
  }

  ///
  /// get the subquad for the index.
  /// If it doesn't exist, create it
  ///
  NodeNode? getSubnode(int index) {
    if (subnode[index] == null) {
      subnode[index] = createSubnode(index);
    }
    return subnode[index];
  }

  NodeNode createSubnode(int index) {
    // create a new subquad in the appropriate quadrant

    double minx = 0.0;
    double maxx = 0.0;
    double miny = 0.0;
    double maxy = 0.0;

    switch (index) {
      case 0:
        minx = env!.getMinX();
        maxx = centrex;
        miny = env!.getMinY();
        maxy = centrey;
        break;
      case 1:
        minx = centrex;
        maxx = env!.getMaxX();
        miny = env!.getMinY();
        maxy = centrey;
        break;
      case 2:
        minx = env!.getMinX();
        maxx = centrex;
        miny = centrey;
        maxy = env!.getMaxY();
        break;
      case 3:
        minx = centrex;
        maxx = env!.getMaxX();
        miny = centrey;
        maxy = env!.getMaxY();
        break;
    }
    Envelope sqEnv = new Envelope(minx, maxx, miny, maxy);
    NodeNode node = NodeNode(sqEnv, level - 1);
    return node;
  }

  int getLevel() {
    return level;
  }
}

class Root extends NodeBase {
  // the singleton root quad is centred at the origin.
  static final Coordinate origin = Coordinate(0.0, 0.0);

  Root() {}

  ///
  /// Insert an item into the quadtree this is the root of.
  ////
  void insert(Envelope itemEnv, Object item) {
    int index = NodeBase.getSubnodeIndex(itemEnv, origin.x, origin.y);
    // if index is -1, itemEnv must cross the X or Y axis.
    if (index == -1) {
      this.add(item);
      return;
    }

    ///
    /// the item must be contained in one quadrant, so insert it into the
    /// tree for that quadrant (which may not yet exist)
    ////
    NodeNode? node = subnode[index];

    ///
    ///  If the subquad doesn't exist or this item is not contained in it,
    ///  have to expand the tree upward to contain the item.
    ////

    if (node == null || !node.getEnvelope()!.containsEnvelope(itemEnv)) {
      NodeNode? largerNode = NodeNode.createExpanded(node, itemEnv);
      subnode[index] = largerNode;
    }

    ///
    /// At this point we have a subquad which exists and must contain
    /// contains the env for the item.  Insert the item into the tree.
    ////
    insertContained(subnode[index]!, itemEnv, item);
    //System.out.println("depth = " + root.depth() + " size = " + root.size());
    //System.out.println(" size = " + size());
  }

  ///
  /// insert an item which is known to be contained in the tree rooted at
  /// the given QuadNode root.  Lower levels of the tree will be created
  /// if necessary to hold the item.
  ///
  void insertContained(NodeNode tree, Envelope itemEnv, Object item) {
    ///
    /// Do NOT create a new quad for zero-area envelopes - this would lead
    /// to infinite recursion. Instead, use a heuristic of simply returning
    /// the smallest existing quad containing the query
    ///
    bool isZeroX =
        IntervalSize.isZeroWidth(itemEnv.getMinX(), itemEnv.getMaxX());
    bool isZeroY =
        IntervalSize.isZeroWidth(itemEnv.getMinY(), itemEnv.getMaxY());
    NodeBase? node;
    if (isZeroX || isZeroY) {
      node = tree.find(itemEnv);
    } else {
      node = tree.getNode(itemEnv);
    }
    node.add(item);
  }

  bool isSearchMatch(Envelope? searchEnv) {
    if (searchEnv == null) {
      return false;
    }
    return true;
  }
}

class Key {
  static int computeQuadLevel(Envelope env) {
    double dx = env.getWidth();
    double dy = env.getHeight();
    double dMax = dx > dy ? dx : dy;
    int level = DoubleBits.exponent(dMax) + 1;
    return level;
  }

  // the fields which make up the key
  Coordinate pt = Coordinate.empty2D();
  int level = 0;
  // auxiliary data which is derived from the key for use in computation
  late Envelope env;

  Key.fromEnvelope(Envelope itemEnv) {
    computeKey(itemEnv);
  }

  Coordinate getPoint() {
    return pt;
  }

  int getLevel() {
    return level;
  }

  Envelope getEnvelope() {
    return env;
  }

  Coordinate getCentre() {
    return new Coordinate((env.getMinX() + env.getMaxX()) / 2,
        (env.getMinY() + env.getMaxY()) / 2);
  }

  ///
  /// return a square envelope containing the argument envelope,
  /// whose extent is a power of two and which is based at a power of 2
  ///
  void computeKey(Envelope itemEnv) {
    level = computeQuadLevel(itemEnv);
    env = Envelope.empty();
    computeKeyLevel(level, itemEnv);
    // MD - would be nice to have a non-iterative form of this algorithm
    while (!env.containsEnvelope(itemEnv)) {
      level += 1;
      computeKeyLevel(level, itemEnv);
    }
  }

  void computeKeyLevel(int level, Envelope itemEnv) {
    double quadSize = math.pow(2.0, level) as double;
    pt.x = (itemEnv.getMinX() / quadSize).floorToDouble() * quadSize;
    pt.y = (itemEnv.getMinY() / quadSize).floorToDouble() * quadSize;
    env.init(pt.x, pt.x + quadSize, pt.y, pt.y + quadSize);
  }
}

class DoubleBits {
  static final int EXPONENT_BIAS = 1023;
  int doubleToBits(double value) {
    const pow52 = 4503599627370496.0; // 2^52
    const pow1022 = 4.49423283715579e+307; // 2^1022
    if (value.isNaN) {
      return 0x7FF8000000000000;
    }
    int signbit = 0;
    if (value.isNegative) {
      signbit = 0x8000000000000000;
      value = -value;
    }
    if (value.isInfinite) {
      return signbit | 0x7FF0000000000000;
    } else if (value < 2.2250738585072014e-308) {
// Denormal or zero.
// Multiply by 2^(1022+52) to get the bits into the correct position.
      int bits = (value * pow1022 * pow52).toInt();
      return signbit | bits;
    } else {
// Slow linear search to move bits into correct position for mantissa.
// Use binary search or something even smarter for speed.
      int exponent = 52;
      while (value < pow52) {
        value *= 2;
        exponent -= 1;
      }
      while (value >= pow52 * 2) {
        value /= 2;
        exponent += 1;
      }
      int mantissaBits = (value - pow52).toInt();
      int exponentBits = (exponent + 1023);
      return signbit | (exponentBits << 52) | mantissaBits;
    }
  }

  static double powerOf2(int exp) {
    if (exp > 1023 || exp < -1022)
      throw new ArgumentError("Exponent out of bounds");
    int expBias = exp + EXPONENT_BIAS;
    int bits = expBias << 52;
    return bits.toDouble();
  }

  static int exponent(double d) {
    DoubleBits db = new DoubleBits.fromDouble(d);
    return db.getExponent();
  }

  static double truncateToPowerOfTwo(double d) {
    DoubleBits db = new DoubleBits.fromDouble(d);
    db.zeroLowerBits(52);
    return db.getDouble();
  }

  static String toBinaryString(double d) {
    DoubleBits db = new DoubleBits.fromDouble(d);
    return db.toString();
  }

  static double maximumCommonMantissa(double d1, double d2) {
    if (d1 == 0.0 || d2 == 0.0) return 0.0;

    DoubleBits db1 = new DoubleBits.fromDouble(d1);
    DoubleBits db2 = new DoubleBits.fromDouble(d2);

    if (db1.getExponent() != db2.getExponent()) return 0.0;

    int maxCommon = db1.numCommonMantissaBits(db2);
    db1.zeroLowerBits(64 - (12 + maxCommon));
    return db1.getDouble();
  }

  double x;
  int? xBits;

  DoubleBits.fromDouble(this.x) {
    xBits = doubleToBits(x);
  }

  double getDouble() {
    return xBits!.toDouble();
  }

  ///
  /// Determines the exponent for the number
  ////
  int biasedExponent() {
    int signExp = (xBits! >> 52);
    int exp = signExp & 0x07ff;
    return exp;
  }

  ///
  /// Determines the exponent for the number
  ////
  int getExponent() {
    return biasedExponent() - EXPONENT_BIAS;
  }

  void zeroLowerBits(int nBits) {
    int invMask = (1 << nBits) - 1;
    int mask = ~invMask;
    xBits = (xBits! & mask);
  }

  int getBit(int i) {
    int mask = (1 << i);
    return (xBits! & mask) != 0 ? 1 : 0;
  }

  ///
  /// This computes the number of common most-significant bits in the mantissa.
  /// It does not count the hidden bit, which is always 1.
  /// It does not determine whether the numbers have the same exponent - if they do
  /// not, the value computed by this function is meaningless.
  /// @param db
  /// @return the number of common most-significant mantissa bits
  ///
  int numCommonMantissaBits(DoubleBits db) {
    for (int i = 0; i < 52; i++) {
      if (getBit(i) != db.getBit(i)) return i;
    }
    return 52;
  }

  ///
  /// A representation of the Double bits formatted for easy readability
  ///
  String toString() {
    String numStr = xBits.toString();
    // 64 zeroes!
    String zero64 =
        "0000000000000000000000000000000000000000000000000000000000000000";
    String padStr = zero64 + numStr;
    String bitStr = padStr.substring(padStr.length - 64);
    String str = bitStr.substring(0, 1) +
        "  " +
        bitStr.substring(1, 12) +
        "(" +
        getExponent().toString() +
        ") " +
        bitStr.substring(12) +
        " [ " +
        x.toString() +
        " ]";
    return str;
  }
}

class IntervalSize {
  ///
  /// This value is chosen to be a few powers of 2 less than the
  /// number of bits available in the double representation (i.e. 53).
  /// This should allow enough extra precision for simple computations to be correct,
  /// at least for comparison purposes.
  ///
  static final int MIN_BINARY_EXPONENT = -50;

  ///
  /// Computes whether the interval [min, max] is effectively zero width.
  /// I.e. the width of the interval is so much less than the
  /// location of the interval that the midpoint of the interval cannot be
  /// represented precisely.
  ////
  static bool isZeroWidth(double minValue, double maxValue) {
    double width = maxValue - minValue;
    if (width == 0.0) return true;
    minValue = minValue.abs();
    maxValue = maxValue.abs();
    double maxAbs = maxValue >= minValue ? maxValue : minValue;
    double scaledInterval = width / maxAbs;
    int level = DoubleBits.exponent(scaledInterval);
    return level <= MIN_BINARY_EXPONENT;
  }
}
