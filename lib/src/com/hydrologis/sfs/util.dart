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