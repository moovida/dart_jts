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

class Stringutils {
  static String replaceCharAt(String oldString, int index, String newChar) {
    return oldString.substring(0, index) +
        newChar +
        oldString.substring(index + 1);
  }
}

class MatrixUtils {
  static List<List<int>> createIntMatrix(int rows, int cols) {
    final grid = List<List<int>>.generate(
        rows, (i) => List<int>.generate(cols, (j) => i * cols + j));
    return grid;
  }
}
