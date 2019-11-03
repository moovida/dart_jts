part of dart_sfs;

_require(cond, [msg]) {
  if (!cond) throw new ArgumentError(msg);
}

//abstract class ComparableMixin<T> implements Comparable {
//  bool operator ==(T other) => compareTo(other) == 0;
//  bool operator <(T other) => compareTo(other) == -1;
//  bool operator <=(T other) => compareTo(other) <= 0;
//  bool operator >(T other) => compareTo(other) == 1;
//  bool operator >=(T other) => compareTo(other) > 0;
//  int compareTo(T other);
//}

class FormattingUtils {}

/// Utilities for collections
class CollectionsUtils {
  /// Shift a list of items from a given [index] to the first position.
  static List<T> shiftToFirst<T>(List<T> list, int index) {
    if (list == null || list.isEmpty) return list;
    if (index > list.length - 1) throw ArgumentError("The shift index can't be bigger than the list size.");
    return list.sublist(index)..addAll(list.sublist(0, index));
  }

  /// Removes subsequent equal items from a [list].
  static List<T> removeRepeated<T>(List<T> list) {
    List<T> newList = [];
    for (int i = 1; i < list.length; i++) {
      if (list[i - 1] != list[i]) {
        newList.add(list[i - 1]);
      }
    }
    if (newList.last != list.last) {
      newList.add(list.last);
    }
    return newList;
  }

  /// CHecks if there are subsequent repeating items.
  static bool hasRepeated<T>(List<T> list) {
    for (int i = 1; i < list.length; i++) {
      if (list[i - 1] == list[i]) {
        return true;
      }
    }
    return false;
  }

  static bool addIfNotEqualToLast<T>(List<T> list, T item) {
    if (list.isEmpty || list.last == item) {
      list.add(item);
    }
  }

  static bool areEqual<T>(List<T> listA, List<T> listB) {
    Function eq = const ListEquality().equals;
    return eq(listA, listB);
  }
}

class StringUtils {
  static bool equalsIgnoreCase(String string1, String string2) {
    return string1?.toLowerCase() == string2?.toLowerCase();
  }

  static String replaceCharAt(String oldString, int index, String newChar) {
    return oldString.substring(0, index) + newChar + oldString.substring(index + 1);
  }
}

class MatrixUtils {
  static List<List<int>> createIntMatrix(int rows, int cols) {
    final grid = List<List<int>>.generate(rows, (i) => List<int>.generate(cols, (j) => i * cols + j));
    return grid;
  }
}

class MathUtils {
  /// Clamps a <tt>double</tt> value to a given range.
  /// @param x the value to clamp
  /// @param min the minimum value of the range
  /// @param max the maximum value of the range
  /// @return the clamped value
  static num clamp(num x, num min, num max) {
    if (x < min) return min;
    if (x > max) return max;
    return x;
  }
}

class NumberUtils {
  static bool equalsWithTolerance(double x1, double x2, double tolerance) {
    return (x1 - x2).abs() <= tolerance;
  }
}
