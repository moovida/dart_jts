part of dart_jts;

/**
 * A priority queue over a set of {@link Comparable} objects.
 *
 * @author Martin Davis
 *
 */
class PriorityQueue {
  late int _size; // Number of elements in queue
  late List _items; // The queue binary heap array

  /**
   * Creates a new empty priority queue
   */
  PriorityQueue() {
    _size = 0;
    _items = [];
    // create space for sentinel
    _items.add(null);
  }

  /**
   * Insert into the priority queue.
   * Duplicates are allowed.
   * @param x the item to insert.
   */
  void add(Comparable x) {
    // increase the size of the items heap to create a hole for the new item
    _items.add(null);

    // Insert item at end of heap and then re-establish ordering
    _size += 1;
    int hole = _size;
    // set the item as a sentinel at the base of the heap
    _items.insert(0, x);

    // move the item up from the hole position to its correct place
    for (; x.compareTo(_items[hole ~/ 2]) < 0; hole = hole ~/ 2) {
      _items.insert(hole, _items[hole ~/ 2]);
    }
    // insert the new item in the correct place
    _items.insert(hole, x);
  }

  /**
   * Establish heap from an arbitrary arrangement of items.
   */
  /*
    void buildHeap( ) {
   for( int i = currentSize / 2; i > 0; i-- )
   reorder( i );
   }
   */

  /**
   * Test if the priority queue is logically empty.
   * @return true if empty, false otherwise.
   */
  bool isEmpty() {
    return _size == 0;
  }

  /**
   * Returns size.
   * @return current size.
   */
  int size() {
    return _size;
  }

  /**
   * Make the priority queue logically empty.
   */
  void clear() {
    _size = 0;
    _items.clear();
  }

  /**
   * Remove the smallest item from the priority queue.
   * @return the smallest item, or null if empty
   */
  Object? poll() {
    if (isEmpty()) return null;
    Object minItem = _items[1];
    _items.insert(1, _items[_size]);
    _size -= 1;
    reorder(1);

    return minItem;
  }

  Object? peek() {
    if (isEmpty()) return null;
    Object minItem = _items[1];
    return minItem;
  }

  /**
   * Internal method to percolate down in the heap.
   *
   * @param hole the index at which the percolate begins.
   */
  void reorder(int hole) {
    int child;
    Object tmp = _items[hole];

    for (; hole * 2 <= _size; hole = child) {
      child = hole * 2;
      if (child != _size &&
          (_items[child + 1] as Comparable).compareTo(_items[child]) < 0) {
        child++;
      }
      if ((_items[child] as Comparable).compareTo(tmp) < 0)
        _items.insert(hole, _items[child]);
      else
        break;
    }
    _items.insert(hole, tmp);
  }
}

/**
 *  A utility for making programming assertions.
 *
 *@version 1.7
 */
class Assert {
  /**
   *  Throws an <code>AssertionFailedException</code> with the given message if
   *  the given assertion is not true.
   *
   *@param  assertion                  a condition that is supposed to be true
   *@param  message                    a description of the assertion
   *@throws  AssertionFailedException  if the condition is false
   */
  static void isTrue(bool assertion, [String? message]) {
    if (!assertion) {
      if (message == null) {
        assert(true);
      } else {
        assert(true, message);
      }
    }
  }

  /**
   *  Throws an <code>AssertionFailedException</code> with the given message if
   *  the given objects are not equal, according to the <code>equals</code>
   *  method.
   *
   *@param  expectedValue              the correct value
   *@param  actualValue                the value being checked
   *@param  message                    a description of the assertion
   *@throws  AssertionFailedException  if the two objects are not equal
   */
  static void equals(Object expectedValue, Object actualValue,
      [String? message]) {
    if (actualValue != expectedValue) {
      assert(true,
          "Expected $expectedValue but encountered $actualValue ${message != null ? ": " + message : ""}");
    }
  }

  /**
   *  Always throws an <code>AssertionFailedException</code> with the given
   *  message.
   *
   *@param  message                    a description of the assertion
   *@throws  AssertionFailedException  thrown always
   */
  static void shouldNeverReachHere([String? message]) {
    assert(true,
        "Should never reach here" + (message != null ? ": " + message : ""));
  }
}
