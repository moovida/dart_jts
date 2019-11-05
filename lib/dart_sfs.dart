/// Support for doing something awesome.
///
/// More dartdocs go here.
library dart_sfs;

import "dart:collection";
import "dart:core" ;
import 'package:intl/intl.dart';
import "package:meta/meta.dart";
import "dart:math" as math;
import "dart:convert" as json;
import 'package:collection/collection.dart';

part "src/com/hydrologis/sfs/io.dart";
part "src/com/hydrologis/sfs/annotations.dart";
part 'src/com/hydrologis/sfs/geom/coordinate.dart';
part 'src/com/hydrologis/sfs/geom/envelope.dart';
part 'src/com/hydrologis/sfs/geom/geom.dart';
part 'src/com/hydrologis/sfs/geom/util.dart';
part 'src/com/hydrologis/sfs/algorithm/distance.dart';
part 'src/com/hydrologis/sfs/algorithm/algorithm.dart';
part 'src/com/hydrologis/sfs/algorithm/locate.dart';
part 'src/com/hydrologis/sfs/operation/operation.dart';
part 'src/com/hydrologis/sfs/operation/valid.dart';
part 'src/com/hydrologis/sfs/operation/relate.dart';
part 'src/com/hydrologis/sfs/operation/overlay.dart';
part 'src/com/hydrologis/sfs/operation/predicate.dart';
part "src/com/hydrologis/sfs/ops.dart";
part "src/com/hydrologis/sfs/math/math.dart";
part "src/com/hydrologis/sfs/geomgraph/geomgraph.dart";
part "src/com/hydrologis/sfs/geomgraph/index.dart";
part "src/com/hydrologis/sfs/util.dart";
part 'src/com/hydrologis/sfs/geom/geometry.dart';
part "src/com/hydrologis/sfs/point.dart";
part "src/com/hydrologis/sfs/geometry_collection.dart";
part "src/com/hydrologis/sfs/multipoint.dart";
part "src/com/hydrologis/sfs/linestring.dart";
part "src/com/hydrologis/sfs/multilinestring.dart";
part "src/com/hydrologis/sfs/surface.dart";
part "src/com/hydrologis/sfs/polygon.dart";
part "src/com/hydrologis/sfs/multipolygon.dart";
part "src/com/hydrologis/sfs/polyhedral_surface.dart";
part "src/com/hydrologis/sfs/util/geom_impl.dart";
part "src/com/hydrologis/sfs/util/util.dart";
part "src/com/hydrologis/sfs/index/index.dart";
part "src/com/hydrologis/sfs/index/strtree.dart";
part "src/com/hydrologis/sfs/util/avltree.dart";


