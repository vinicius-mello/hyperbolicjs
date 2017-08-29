// Generated by CoffeeScript 1.10.0
(function() {
  var Complex, Isometry, acosh, asinh, circle, circleSVG, delaunayEmbedding, delaunayEmbeddingRec, distance, is0, isBoundary, isCollinear0, line, lineCenter, midpoint, midpoint0, midpointIdeal, perpBisector, perpBisectorHalfPlaneSVG, perpBisectorSVG, polygonSVG, randomPoint, ray, reflect, regularPolygon, regularPolygonSVG, regularTiling, regularTilingRadius, root, sArea, segmentSVG, svg, toDisc, translate0, translate0Inv;

  Complex = (function() {
    function Complex(x1, y) {
      this.x = x1 != null ? x1 : 0.0;
      this.y = y != null ? y : 0.0;
    }

    Complex.zero = new Complex(0.0, 0.0);

    Complex.one = new Complex(1.0, 0.0);

    Complex.i = new Complex(0.0, 1.0);

    Complex.prototype.plus = function(w) {
      return new Complex(this.x + w.x, this.y + w.y);
    };

    Complex.prototype.minus = function(w) {
      return new Complex(this.x - w.x, this.y - w.y);
    };

    Complex.prototype.times = function(w) {
      return new Complex(this.x * w.x - this.y * w.y, this.x * w.y + this.y * w.x);
    };

    Complex.prototype.divide = function(w) {
      var w2;
      w2 = w.norm2();
      if (w2 === 0.0) {
        throw Error("no inverse");
      }
      return new Complex((this.x * w.x + this.y * w.y) / w2, (this.y * w.x - this.x * w.y) / w2);
    };

    Complex.prototype.negation = function() {
      return new Complex(-1.0 * this.x, -1.0 * this.y);
    };

    Complex.prototype.perp = function() {
      return new Complex(-1.0 * this.y, this.x);
    };

    Complex.prototype.normalized = function() {
      var n2;
      n2 = this.x * this.x + this.y * this.y;
      if (n2 === 0.0) {
        throw Error("no inverse");
      }
      return new Complex(this.x / Math.sqrt(n2), this.y / Math.sqrt(n2));
    };

    Complex.prototype.inverse = function() {
      var n2;
      n2 = this.x * this.x + this.y * this.y;
      if (n2 === 0.0) {
        throw Error("no inverse");
      }
      return new Complex(this.x / n2, -1.0 * this.y / n2);
    };

    Complex.prototype.conjugate = function() {
      return new Complex(this.x, -this.y);
    };

    Complex.prototype.scale = function(s) {
      return new Complex(s * this.x, s * this.y);
    };

    Complex.prototype.dot = function(w) {
      return this.x * w.x + this.y * w.y;
    };

    Complex.prototype.norm = function() {
      var n2;
      n2 = this.x * this.x + this.y * this.y;
      return Math.sqrt(n2);
    };

    Complex.prototype.norm2 = function() {
      return this.x * this.x + this.y * this.y;
    };

    Complex.setPrecision = function(n) {
      Complex.precision = n;
      Complex.precisionStr1 = '0.' + '0'.repeat(n);
      return Complex.precisionStr2 = '-' + Complex.precisionStr1;
    };

    Complex.precision = 6;

    Complex.precisionStr1 = '0.000000';

    Complex.precisionStr2 = '-0.000000';

    Complex.prototype.toString = function() {
      var nx, ny, pre;
      pre = Complex.precision;
      nx = this.x.toFixed(pre);
      if (nx === Complex.precisionStr2) {
        nx = Complex.precisionStr1;
      }
      ny = this.y.toFixed(pre);
      if (ny === Complex.precisionStr2) {
        ny = Complex.precisionStr1;
      }
      return "(" + nx + "," + ny + ")";
    };

    return Complex;

  })();

  Isometry = (function() {
    function Isometry(m) {
      if (typeof m === "object") {
        this.a = m.a;
        this.b = m.b;
        this.c = m.c;
        this.d = m.d;
        this.orientation = m.orientation;
      } else {
        this.a = Complex.one;
        this.b = Complex.zero;
        this.c = Complex.zero;
        this.d = Complex.one;
        this.orientation = 1.0;
      }
    }

    Isometry.prototype.apply = function(z) {
      if (this.orientation < 0) {
        z = z.conjugate();
      }
      return z.times(this.a).plus(this.b).divide(z.times(this.c).plus(this.d));
    };

    Isometry.prototype.compose = function(m) {
      var t;
      if (this.orientation < 0) {
        m = m.conjugate();
      }
      t = new Isometry();
      t.a = this.a.times(m.a).plus(this.b.times(m.c));
      t.b = this.a.times(m.b).plus(this.b.times(m.d));
      t.c = this.c.times(m.a).plus(this.d.times(m.c));
      t.d = this.c.times(m.b).plus(this.d.times(m.d));
      t.orientation = this.orientation * m.orientation;
      return t;
    };

    Isometry.prototype.conjugate = function() {
      var t;
      t = new Isometry();
      t.a = this.a.conjugate();
      t.b = this.b.conjugate();
      t.c = this.c.conjugate();
      t.d = this.d.conjugate();
      t.orientation = this.orientation;
      return t;
    };

    Isometry.translate0 = function(a) {
      var t;
      t = new Isometry();
      t.a = Complex.one;
      t.b = a.negation();
      t.c = t.b.conjugate();
      t.d = t.a;
      return t;
    };

    Isometry.translate0Inv = function(a) {
      return this.translate0(a.negation());
    };

    Isometry.translate = function(a, b) {
      var s, t;
      t = this.translate0(a);
      s = this.translate0(b.negation());
      return s.compose(t);
    };

    Isometry.rotate0 = function(a, b) {
      var e, t;
      e = null;
      if (typeof a === "number") {
        e = new Complex(Math.cos(a), Math.sin(a));
      } else {
        if (b != null) {
          a = a.normalized();
          b = b.normalized();
          e = b.divide(a);
        } else {
          e = a.normalized();
        }
      }
      t = new Isometry();
      t.a = e;
      t.b = Complex.zero;
      t.c = t.b;
      t.d = Complex.one;
      return t;
    };

    Isometry.reflect = function(p, q) {
      var a, b, ce, l, r2, t;
      t = new Isometry();
      t.orientation = -1.0;
      if (isCollinear0(p, q)) {
        if (p.norm2() < q.norm2()) {
          t.a = q.divide(q.conjugate());
          return t;
        } else {
          t.a = p.divide(p.conjugate());
          return t;
        }
      }
      l = line(p, q);
      a = l[0], b = l[1];
      ce = lineCenter(l);
      r2 = ce.minus(a).norm2();
      t.a = ce;
      t.b = new Complex(r2 - ce.norm2(), 0.0);
      t.c = Complex.one;
      t.d = ce.conjugate().negation();
      return t;
    };

    Isometry.prototype.force = function() {
      var s;
      s = 1.0 / this.a.norm();
      this.a = this.a.scale(s);
      this.b = this.b.scale(s);
      this.c = this.c.scale(s);
      this.d = this.d.scale(s);
      return this;
    };

    return Isometry;

  })();

  is0 = function(z) {
    return z.norm2() < 0.000000000001;
  };

  isBoundary = function(z) {
    return Math.abs(1.0 - z.norm2()) < 0.000000000001;
  };

  isCollinear0 = function(p, q) {
    return Math.abs(p.x * q.y - p.y * q.x) < 0.000001;
  };

  toDisc = function(z) {
    var n;
    n = z.norm2();
    if (n > 1) {
      return z.scale(1.0 / Math.sqrt(n));
    } else {
      return z;
    }
  };

  translate0 = function(a, z) {
    var w;
    w = Complex.one;
    w = z.minus(a).divide(w.minus(z.times(a.conjugate())));
    return w;
  };

  translate0Inv = function(a, z) {
    var w;
    w = Complex.one;
    w = z.plus(a).divide(w.plus(z.times(a.conjugate())));
    return w;
  };

  midpoint0 = function(a) {
    var a2, s;
    a2 = a.norm2();
    s = (1.0 - Math.sqrt(1.0 - a2)) / a2;
    return a.scale(s);
  };

  midpointIdeal = function(a, b) {
    var co, s;
    co = a.dot(b);
    s = a.plus(b);
    return s.scale(1.0 / (2.0 * (1.0 + Math.sqrt(Math.abs(1.0 - co) / 2.0))));
  };

  midpoint = function(a, b) {
    return translate0Inv(a, midpoint0(translate0(a, b)));
  };

  circle = function(O, R) {
    var O2, Ox, Oy, p, p2, r;
    p = Math.tanh(R / 2.0);
    p2 = p * p;
    O2 = O.norm2();
    Ox = (O.x * (1.0 - p2)) / (1.0 - p2 * O2);
    Oy = (O.y * (1.0 - p2)) / (1.0 - p2 * O2);
    r = (p * (1.0 - O2)) / (1.0 - p2 * O2);
    return [new Complex(Ox, Oy), r];
  };

  line = function(p, q) {
    var a, b, bp, bq;
    if (isCollinear0(p, q)) {
      if (p.norm2() < q.norm2()) {
        b = q;
        b = b.normalized();
        a = b.negation();
      } else {
        a = p;
        a = a.normalized();
        b = a.negation();
      }
      return [a, b];
    }
    bp = isBoundary(p);
    bq = isBoundary(q);
    if (bp && bq) {
      a = p;
      b = q;
    } else if (!bp) {
      b = translate0(p, q);
      b = b.normalized();
      a = b.negation();
      a = translate0Inv(p, a);
      b = translate0Inv(p, b);
    } else if (!bq) {
      a = translate0(q, p);
      a = a.normalized();
      b = a.negation();
      a = translate0Inv(q, a);
      b = translate0Inv(q, b);
    }
    return [a, b];
  };

  lineCenter = function(l) {
    var a, b, ce, s;
    a = l[0], b = l[1];
    s = a.plus(b);
    ce = s.scale(1.0 / (1.0 + a.dot(b)));
    return ce;
  };

  reflect = function(p, q, z) {
    var a, b, ce, l, r2, v, v2;
    if (isCollinear0(p, q)) {
      if (p.norm2() < q.norm2()) {
        return q.times(z.divide(q).conjugate());
      } else {
        return p.times(z.divide(p).conjugate());
      }
    }
    l = line(p, q);
    a = l[0], b = l[1];
    ce = lineCenter(l);
    r2 = ce.minus(a).norm2();
    v = z.minus(ce);
    v2 = v.norm2();
    return ce.plus(v.scale(r2 / v2));
  };

  perpBisector = function(p, q) {
    var a, b, m, m2, ref;
    ref = line(p, q), a = ref[0], b = ref[1];
    m = midpointIdeal(a, b);
    p = translate0(m, p);
    q = translate0(m, q);
    m2 = midpoint(p, q);
    q = translate0(m2, q);
    b = q.perp();
    b = b.normalized();
    a = b.negation();
    a = translate0Inv(m, translate0Inv(m2, a));
    b = translate0Inv(m, translate0Inv(m2, b));
    return [a, b];
  };

  ray = function(p, q) {
    var a, b, ref;
    ref = line(p, q), a = ref[0], b = ref[1];
    return [p, b];
  };

  distance = function(z, w) {
    var one;
    one = Complex.one;
    return 2.0 * Math.atanh((z.minus(w).divide(one.minus(z.times(w.conjugate())))).norm());
  };

  regularPolygon = function(c, v, n) {
    var a, i, vp, vs;
    v = translate0(c, v);
    vp = v.perp();
    a = 2.0 * Math.PI / n;
    vs = (function() {
      var j, ref, results;
      results = [];
      for (i = j = 0, ref = n; 0 <= ref ? j < ref : j > ref; i = 0 <= ref ? ++j : --j) {
        results.push(v.scale(Math.cos(i * a)).plus(vp.scale(Math.sin(i * a))));
      }
      return results;
    })();
    return (function() {
      var j, len1, results;
      results = [];
      for (j = 0, len1 = vs.length; j < len1; j++) {
        v = vs[j];
        results.push(translate0Inv(c, v));
      }
      return results;
    })();
  };

  sArea = function(a, b, c) {
    return a.x * b.y + b.x * c.y + c.x * a.y - b.x * a.y - c.x * b.y - a.x * c.y;
  };

  asinh = function(x) {
    return Math.log(x + Math.sqrt(x * x + 1.0));
  };

  acosh = function(x) {
    return Math.log(x + Math.sqrt(-1.0 + x * x));
  };

  segmentSVG = function(p, q, move) {
    var ce, flag, l, r, str;
    if (move == null) {
      move = true;
    }
    if (isCollinear0(p, q)) {
      str = "L" + q.x + "," + q.y;
      if (move) {
        str = ("M" + p.x + "," + p.y + " ") + str;
      }
      return str;
    }
    l = line(p, q);
    ce = lineCenter(l);
    r = ce.minus(p).norm();
    flag = sArea(p, q, ce) > 0 ? 1 : 0;
    str = "A" + r + "," + r + " 0 0," + flag + " " + q.x + "," + q.y;
    if (move) {
      str = ("M" + p.x + "," + p.y + " ") + str;
    }
    return str;
  };

  circleSVG = function(O, R) {
    var o, r, ref;
    ref = circle(O, R), o = ref[0], r = ref[1];
    return "M " + o.x + " " + o.y + " m -" + r + ", 0 a " + r + "," + r + " 0 1,0 " + (r * 2) + ",0 a " + r + "," + r + " 0 1,0 " + (-(r * 2)) + ",0";
  };

  perpBisectorSVG = function(p, q) {
    var a, b, ref;
    ref = perpBisector(p, q), a = ref[0], b = ref[1];
    return segmentSVG(a, b);
  };

  perpBisectorHalfPlaneSVG = function(p, q) {
    var a, b, flag, ref, zero;
    ref = perpBisector(p, q), a = ref[0], b = ref[1];
    zero = Complex.zero;
    flag = sArea(a, b, zero) < 0 ? 1 : 0;
    return segmentSVG(a, b) + (" A 1.0,1.0 0 " + flag + ",0 " + a.x + "," + a.y + "Z");
  };

  polygonSVG = function(vs) {
    var i, n, p, path;
    n = vs.length;
    p = vs[0];
    path = (function() {
      var j, ref, results;
      results = [];
      for (i = j = 0, ref = n; 0 <= ref ? j < ref : j > ref; i = 0 <= ref ? ++j : --j) {
        results.push(segmentSVG(vs[i], vs[(i + 1) % n], false));
      }
      return results;
    })();
    return ("M" + p.x + "," + p.y) + path.join("") + "Z";
  };

  regularPolygonSVG = function(c, v, n) {
    var vs;
    vs = regularPolygon(c, v, n);
    return polygonSVG(vs);
  };

  svg = {
    segment: segmentSVG,
    circle: circleSVG,
    perpBisector: perpBisectorSVG,
    perpBisectorHalfPlane: perpBisectorHalfPlaneSVG,
    polygon: polygonSVG,
    regularPolygon: regularPolygonSVG
  };

  randomPoint = function(R, alpha) {
    var r, s, th;
    if (alpha == null) {
      alpha = 1.0;
    }
    th = 2.0 * Math.PI * Math.random();
    s = Math.random();
    r = Math.acosh((Math.cosh(alpha * R) - 1.0) * s + 1.0) / alpha;
    r = Math.tanh(r / 2.0);
    return new Complex(r * Math.cos(th), r * Math.sin(th));
  };

  regularTilingRadius = function(n, m) {
    var alpha, beta, chc, len;
    alpha = Math.PI / n;
    beta = Math.PI / m;
    chc = 1.0 / (Math.tan(alpha) * Math.tan(beta));
    len = Math.sqrt(chc * chc - 1.0) / (1.0 + chc);
    return len;
  };

  regularTiling = function(n, m, stop, adjacency) {
    var i, j, k, len, len1, p, polygons, q, ref, ref1, reflectSide, tovisit, u, visited, vs;
    if (adjacency == null) {
      adjacency = false;
    }
    reflectSide = function(poly, vs, i) {
      var c, hop, nt, p, q, t, v;
      c = poly.c, v = poly.v, t = poly.t, hop = poly.hop;
      p = vs[i];
      q = vs[(i + 1) % n];
      nt = Isometry.reflect(p, q);
      return {
        c: reflect(p, q, c),
        v: p,
        t: nt.compose(t),
        hop: hop + 1
      };
    };
    len = regularTilingRadius(n, m);
    tovisit = [
      {
        c: Complex.zero,
        v: new Complex(len, 0.0),
        t: new Isometry(),
        hop: 0
      }
    ];
    visited = d3.map();
    polygons = [];
    while (tovisit.length > 0) {
      p = tovisit.shift();
      if (!visited.has(p.c)) {
        if (!stop(p)) {
          visited.set(p.c, polygons.length);
          polygons.push(p);
          vs = regularPolygon(p.c, p.v, n);
          for (i = j = 0, ref = n; 0 <= ref ? j < ref : j > ref; i = 0 <= ref ? ++j : --j) {
            tovisit.push(reflectSide(p, vs, i));
          }
        }
      }
    }
    adjacency = true;
    if (adjacency) {
      console.log("1");
      for (k = 0, len1 = polygons.length; k < len1; k++) {
        p = polygons[k];
        console.log("2");
        vs = regularPolygon(p.c, p.v, n);
        p.adjacency = [];
        for (i = u = 0, ref1 = n; 0 <= ref1 ? u < ref1 : u > ref1; i = 0 <= ref1 ? ++u : --u) {
          q = reflectSide(p, vs, i);
          if (visited.has(q.c)) {
            p.adjacency.push(visited.get(q.c));
          }
        }
      }
    }
    return polygons;
  };

  delaunayEmbedding = function(root) {
    return delaunayEmbeddingRec(root, 0, new Isometry());
  };

  delaunayEmbeddingRec = function(node, alpha, m) {
    var calpha, child, childSpan, i, j, mm, nchildren, r, ref, results, rot, span, theta, tr, z, zero;
    zero = Complex.zero;
    z = m.apply(zero);
    node.x = z.x;
    node.y = z.y;
    if (node.children == null) {
      return;
    }
    nchildren = node.children.length;
    if (nchildren === 0) {
      return;
    }
    span = 2.0 * (Math.PI - alpha);
    childSpan = span / nchildren;
    calpha = Math.min(Math.PI / 3.0, childSpan / 2.0);
    r = Math.cos(calpha);
    tr = Isometry.translate0Inv(new Complex(r, 0));
    results = [];
    for (i = j = 0, ref = nchildren; 0 <= ref ? j < ref : j > ref; i = 0 <= ref ? ++j : --j) {
      theta = -span / 2.0 + calpha + i * childSpan;
      rot = Isometry.rotate0(theta);
      mm = m.compose(rot.compose(tr));
      child = node.children[i];
      results.push(delaunayEmbeddingRec(child, calpha, mm));
    }
    return results;
  };

  root = self.hyperbolic != null ? self.hyperbolic : self.hyperbolic = {};

  root.Complex = Complex;

  root.Isometry = Isometry;

  root.segmentSVG = segmentSVG;

  root.circleSVG = circleSVG;

  root.reflect = reflect;

  root.circle = circle;

  root.perpBisectorSVG = perpBisectorSVG;

  root.perpBisectorHalfPlaneSVG = perpBisectorHalfPlaneSVG;

  root.randomPoint = randomPoint;

  root.regularPolygon = regularPolygon;

  root.regularPolygonSVG = regularPolygonSVG;

  root.polygonSVG = polygonSVG;

  root.svg = svg;

  root.distance = distance;

  root.regularTiling = regularTiling;

  root.regularTilingRadius = regularTilingRadius;

  root.delaunayEmbedding = delaunayEmbedding;

  Math.asinh = asinh;

}).call(this);
