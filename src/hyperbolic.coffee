class Complex
    constructor: (@x=0.0, @y=0.0) ->

    plus: (w) ->
        new Complex(
            @x + w.x,
            @y + w.y
        )

    minus: (w) ->
        new Complex(
            @x - w.x,
            @y - w.y
        )

    times: (w) ->
        new Complex(
            @x*w.x - @y*w.y,
            @x*w.y + @y*w.x
        )

    divide: (w) ->
        w2 = w.norm2()
        throw Error "no inverse" if w2 is 0.0
        new Complex(
            (@x*w.x + @y*w.y)/w2,
            (@y*w.x - @x*w.y)/w2
        )

    negation: ->
        new Complex(
            -1.0 * @x,
            -1.0 * @y
        )

    perp: ->
        new Complex(
            -1.0 * @y,
            @x
        )

    normalized: ->
        n2 = @x*@x + @y*@y
        throw Error "no inverse" if n2 is 0.0
        new Complex(
            @x/Math.sqrt(n2),
            @y/Math.sqrt(n2)
        )

    inverse: ->
        n2 = @x*@x + @y*@y
        throw Error "no inverse" if n2 is 0.0
        new Complex(
            @x / n2,
            -1.0 * @y / n2
        )

    conjugate: ->
        new Complex(
            @x ,
            -@y
        )

    scale: (s) ->
        new Complex(s*@x,s*@y)

    dot: (w) ->
        @x*w.x+@y*w.y

    norm: ->
        n2 = @x*@x + @y*@y
        Math.sqrt(n2)

    norm2: ->
        @x*@x + @y*@y

    toString: ->
        return "#{@x.toFixed(6)}" if @y == 0.0
        return "#{@y.toFixed(6)}i" if @x == 0.0
        if @y > 0
            "#{@x.toFixed(6)}+#{@y.toFixed(6)}i"
        else
            "#{@x.toFixed(6)}-#{(-1.0 * @y).toFixed(6)}i"

class Moebius
    constructor: (m) ->
        if typeof m is "object"
            @a = m.a
            @b = m.b
            @c = m.c
            @d = m.d
        else
            @a = new Complex(1.0,0.0)
            @b = new Complex(0.0,0.0)
            @c = new Complex(0.0,0.0)
            @d = new Complex(1.0,0.0)

    apply: (z) ->
        z.times(@a).plus(@b).divide(z.times(@c).plus(@d))

    applyInv: (z) ->
        z.times(@d).minus(@b).divide(@a.minus(z.times(@c)))

    compose: (m) ->
        t = new Moebius()
        t.a = @a.times(m.a).plus(@b.times(m.c))
        t.b = @a.times(m.b).plus(@b.times(m.d))
        t.c = @c.times(m.a).plus(@d.times(m.c))
        t.d = @c.times(m.b).plus(@d.times(m.d))
        return t

    conjugate: ->
        t = new Moebius()
        t.a = @a.conjugate()
        t.b = @b.conjugate()
        t.c = @c.conjugate()
        t.d = @d.conjugate()
        return t

class Isometry extends Moebius
    orientation: 1.0

    constructor: ->
        super()

    apply: (z) ->
        if @orientation<0
            z = z.conjugate()
        return super(z)

    compose: (m) ->
        if @orientation<0
            m = m.conjugate()
        t = new Isometry()
        t.a = @a.times(m.a).plus(@b.times(m.c))
        t.b = @a.times(m.b).plus(@b.times(m.d))
        t.c = @c.times(m.a).plus(@d.times(m.c))
        t.d = @c.times(m.b).plus(@d.times(m.d))
        t.orientation = @orientation * m.orientation
        return t

    @translate0: (a) ->
        t = new Isometry()
        t.a = new Complex(1.0,0.0)
        t.b = a.negation()
        t.c = t.b.conjugate()
        t.d = t.a
        return t

    @translate0Inv: (a) ->
        return @translate0(a.negation())

    @translate: (a, b) ->
        t = @translate0(a)
        s = @translate0(b.negation())
        return s.compose(t)

    @rotate0: (a,b) ->
        e = null
        if typeof a is "number"
            e =  new Complex(Math.cos(a),Math.sin(a))
        else
            if b?
                a = a.normalized()
                b = b.normalized()
                e = b.divide(a)
            else
                e = a.normalized()
        t = new Isometry()
        t.a = e
        t.b = new Complex(0.0,0.0)
        t.c = t.b
        t.d = new Complex(1.0,0.0)
        return t

    @reflect: (a,b) ->
        t = new Isometry()
        t.orientation = -1.0
        if isCollinear0(p,q)
            if p.norm2()<q.norm2()
                t.a = q.divide(q.conjugate())
                return t
            else
                t.a = p.divide(p.conjugate())
                return t
        l = line(a,b)
        ce = lineCenter(l)
        r2 = ce.minus(a).norm2()
        [a, b] = l
        t.a = ce
        t.b = new Complex(r2-ce.norm2(),0.0)
        t.c = new Complex(1.0,0.0)
        t.d = ce.negation()
        return t.force()

    force: ->
        s= 1.0/@a.norm()
        @a=@a.scale(s)
        @b=@b.scale(s)
        @c=@b.divide(@a).conjugate()
        @d=new Complex(1.0,0.0)
        return this


is0 = (z) ->
    return z.norm2()<0.000000000001

isBoundary = (z) ->
    return Math.abs(1.0-z.norm2())<0.000000000001

isCollinear0 = (p,q) ->
    return Math.abs(p.x*q.y-p.y*q.x)<0.000001

toDisc = (z) ->
    n = z.norm2()
    if n>1
        return z.scale(1.0/Math.sqrt(n))
    else
        return z

translate0 = (a,z) ->
    w = new Complex(1.0,0.0)
    w = z.minus(a).divide(w.minus(z.times(a.conjugate())))
    return w

translate0Inv = (a,z) ->
    w = new Complex(1.0,0.0)
    w = z.plus(a).divide(w.plus(z.times(a.conjugate())))
    return w

midpoint0 = (a) ->
    a2 = a.norm2()
    s = (1.0-Math.sqrt(1.0-a2))/a2
    return a.scale(s)

midpointIdeal = (a,b) ->
    co = a.dot(b)
    s = a.plus(b)
    return s.scale(1.0/(2.0*(1.0+Math.sqrt(Math.abs(1.0-co)/2.0))))

midpoint = (a,b) ->
    return translate0Inv(a,midpoint0(translate0(a,b)))

circle = (O,R) ->
    p = Math.tanh(R/2.0)
    p2 = p*p
    O2 = O.norm2()
    Ox = (O.x*(1.0-p2))/(1.0-p2*O2)
    Oy = (O.y*(1.0-p2))/(1.0-p2*O2)
    r = (p*(1.0-O2))/(1.0-p2*O2)
    return [new Complex(Ox,Oy), r]

line = (p,q) ->
    if isCollinear0(p,q)
        if p.norm2()<q.norm2()
            b = q
            b = b.normalized()
            a = b.negation()
        else
            a = p
            a = a.normalized()
            b = a.negation()
        return [a,b]
    bp = isBoundary(p)
    bq = isBoundary(q)
    if bp && bq
        a = p
        b = q
    else if not bp
        b = translate0(p,q)
        b = b.normalized()
        a = b.negation()
        a = translate0Inv(p,a)
        b = translate0Inv(p,b)
    else if not bq
        a = translate0(q,p)
        a = a.normalized()
        b = a.negation()
        a = translate0Inv(q,a)
        b = translate0Inv(q,b)
    return [a,b]

lineCenter = (l) ->
    [a, b] = l
    s = a.plus(b)
    ce = s.scale(1.0/(1.0+a.dot(b)))
    return ce

reflect = (p,q,z) ->
    if isCollinear0(p,q)
        if p.norm2()<q.norm2()
            return q.times(z.divide(q).conjugate())
        else
            return p.times(z.divide(p).conjugate())
    l = line(p,q)
    [a, b] = l
    ce = lineCenter(l)
    r2 = ce.minus(a).norm2()
    v = z.minus(ce)
    v2 = v.norm2()
    return ce.plus(v.scale(r2/v2))

perpBisector = (p,q) ->
    [a, b] = line(p,q)
    m = midpointIdeal(a,b)
    p = translate0(m,p)
    q = translate0(m,q)
    m2 = midpoint(p,q)
    q = translate0(m2,q)
    b = q.perp()
    b = b.normalized()
    a = b.negation()
    a = translate0Inv(m,translate0Inv(m2,a))
    b = translate0Inv(m,translate0Inv(m2,b))
    return [a, b]

ray = (p,q) ->
    [a, b] = line(p,q)
    return [p, b]

distance = (z,w) ->
    one = new Complex(1.0,0.0)
    return 2.0*Math.atanh((z.minus(w).divide(one.minus(z.times(w.conjugate())))).norm())

regularPolygon = (c,v,n) ->
    v = translate0(c,v)
    vp = v.perp()
    a = 2.0*Math.PI/n
    vs = (v.scale(Math.cos(i*a)).plus(vp.scale(Math.sin(i*a))) for i in [0...n])
    return (translate0Inv(c,v) for v in vs)

sArea = (a,b,c) ->
    return a.x*b.y+b.x*c.y+c.x*a.y - b.x*a.y-c.x*b.y-a.x*c.y

asinh = (x) ->
    return Math.log(x+Math.sqrt(x*x + 1.0))

acosh = (x) ->
    return Math.log(x+Math.sqrt(-1.0+x*x))

segmentSVG = (p,q,move=true) ->
    if isCollinear0(p,q)
        return "M#{p.x},#{p.y} L#{q.x},#{q.y}"
    l = line(p,q)
    ce = lineCenter(l)
    r = ce.minus(p).norm()
    flag = if sArea(p,q,ce)>0 then 1 else 0
    str = "A#{r},#{r} 0 0,#{flag} #{q.x},#{q.y}"
    if move
        str = "M#{p.x},#{p.y} "+str
    return str


perpBisectorSVG = (p,q) ->
    [a, b] = perpBisector(p,q)
    return segmentSVG(a,b)

perpBisectorHalfPlaneSVG = (p,q) ->
    [a, b] = perpBisector(p,q)
    zero = new Complex(0.0,0.0)
    flag = if sArea(a,b,zero)<0 then 1 else 0
    return segmentSVG(a,b)+" A 1.0,1.0 0 #{flag},0 #{a.x},#{a.y}Z"

regularPolygonSVG = (c,v,n) ->
    vs = regularPolygon(c,v,n)
    p = vs[0]
    path = (segmentSVG(vs[i],vs[(i+1)%n],false) for i in [0...n])
    return "M#{p.x},#{p.y}"+path.join("")+"Z"

randomPoint = (R) ->
    th = 2.0 * Math.PI * Math.random()
    s = Math.random()
    r = Math.acosh((Math.cosh(R)-1.0)*s + 1.0)
    r = Math.tanh(r/2.0)
    return new Complex(r*Math.cos(th), r*Math.sin(th))

delaunayEmbedding = (root) ->
    delaunayEmbeddingRec(root,0,new Isometry())

delaunayEmbeddingRec = (node, alpha, m) ->
    zero = new Complex(0.0,0.0)
    z = m.apply(zero)
    node.x = z.x
    node.y = z.y
    if not node.children?
        return
    nchildren = node.children.length
    if nchildren==0
        return
    span = 2.0*(Math.PI - alpha)
    childSpan = span/nchildren
    calpha = Math.min(Math.PI/3.0,childSpan/2.0)
    r = Math.cos(calpha)
    tr = Isometry.translate0Inv(new Complex(r,0))
    for i in [0...nchildren]
        theta = -span/2.0+calpha+i*childSpan
        rot = Isometry.rotate0(theta)
        mm = m.compose(rot.compose(tr))
        child = node.children[i]
        delaunayEmbeddingRec(child,calpha,mm)


root = self.hyperbolic ?= {}
root.Complex = Complex
root.Moebius = Moebius
root.Isometry = Isometry
root.segmentSVG = segmentSVG
root.circle = circle
root.perpBisectorSVG = perpBisectorSVG
root.perpBisectorHalfPlaneSVG = perpBisectorHalfPlaneSVG
root.randomPoint = randomPoint
root.regularPolygon = regularPolygon
root.regularPolygonSVG = regularPolygonSVG
root.distance = distance
Math.asinh = asinh
root.delaunayEmbedding = delaunayEmbedding
