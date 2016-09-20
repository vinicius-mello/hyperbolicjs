class Complex

    constructor: (@x=0.0, @y=0.0) ->

    @zero: new Complex(0.0,0.0)
    @one: new Complex(1.0,0.0)
    @i: new Complex(0.0,1.0)

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

    @setPrecision: (n)->
        Complex.precision = n
        Complex.precisionStr1 = '0.'+'0'.repeat(n)
        Complex.precisionStr2 = '-'+Complex.precisionStr1
    @precision: 6
    @precisionStr1 : '0.000000'
    @precisionStr2 : '-0.000000'
    toString: ->
        pre = Complex.precision
        nx = (@x).toFixed(pre)
        if nx == Complex.precisionStr2
            nx = Complex.precisionStr1
        ny = (@y).toFixed(pre)
        if ny == Complex.precisionStr2
            ny = Complex.precisionStr1
        "(#{nx},#{ny})"


class Isometry

    constructor: (m) ->
        if typeof m is "object"
            @a = m.a
            @b = m.b
            @c = m.c
            @d = m.d
            @orientation = m.orientation
        else
            @a = Complex.one
            @b = Complex.zero
            @c = Complex.zero
            @d = Complex.one
            @orientation = 1.0

    apply: (z) ->
        if @orientation<0
            z = z.conjugate()
        return z.times(@a).plus(@b).divide(z.times(@c).plus(@d))

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

    conjugate: ->
        t = new Isometry()
        t.a = @a.conjugate()
        t.b = @b.conjugate()
        t.c = @c.conjugate()
        t.d = @d.conjugate()
        t.orientation = @orientation
        return t

    @translate0: (a) ->
        t = new Isometry()
        t.a = Complex.one
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
        t.b = Complex.zero
        t.c = t.b
        t.d = Complex.one
        return t

    @reflect: (p,q) ->
        t = new Isometry()
        t.orientation = -1.0
        if isCollinear0(p,q)
            if p.norm2()<q.norm2()
                t.a = q.divide(q.conjugate())
                return t
            else
                t.a = p.divide(p.conjugate())
                return t
        l = line(p,q)
        [a, b] = l
        ce = lineCenter(l)
        r2 = ce.minus(a).norm2()
        t.a = ce
        t.b = new Complex(r2-ce.norm2(),0.0)
        t.c = Complex.one
        t.d = ce.conjugate().negation()
        return t

    force: ->
        s= 1.0/@a.norm()
        @a=@a.scale(s)
        @b=@b.scale(s)
        @c=@c.scale(s)
        @d=@d.scale(s)
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
    w = Complex.one
    w = z.minus(a).divide(w.minus(z.times(a.conjugate())))
    return w

translate0Inv = (a,z) ->
    w = Complex.one
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
    one = Complex.one
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
        str = "L#{q.x},#{q.y}"
        if move
            str = "M#{p.x},#{p.y} "+str
        return str
    l = line(p,q)
    ce = lineCenter(l)
    r = ce.minus(p).norm()
    flag = if sArea(p,q,ce)>0 then 1 else 0
    str = "A#{r},#{r} 0 0,#{flag} #{q.x},#{q.y}"
    if move
        str = "M#{p.x},#{p.y} "+str
    return str

circleSVG = (O,R)->
    [o,r] = circle(O,R)
    return "M #{o.x} #{o.y} m -#{r}, 0 a #{r},#{r} 0 1,0 #{r * 2},0 a #{r},#{r} 0 1,0 #{-(r * 2)},0"

perpBisectorSVG = (p,q) ->
    [a, b] = perpBisector(p,q)
    return segmentSVG(a,b)

perpBisectorHalfPlaneSVG = (p,q) ->
    [a, b] = perpBisector(p,q)
    zero = Complex.zero
    flag = if sArea(a,b,zero)<0 then 1 else 0
    return segmentSVG(a,b)+" A 1.0,1.0 0 #{flag},0 #{a.x},#{a.y}Z"

polygonSVG = (vs) ->
    n = vs.length
    p = vs[0]
    path = (segmentSVG(vs[i],vs[(i+1)%n],false) for i in [0...n])
    return "M#{p.x},#{p.y}"+path.join("")+"Z"

regularPolygonSVG = (c,v,n) ->
    vs = regularPolygon(c,v,n)
    return polygonSVG(vs)

randomPoint = (R) ->
    th = 2.0 * Math.PI * Math.random()
    s = Math.random()
    r = Math.acosh((Math.cosh(R)-1.0)*s + 1.0)
    r = Math.tanh(r/2.0)
    return new Complex(r*Math.cos(th), r*Math.sin(th))

regularTilingRadius = (n,m)->
    alpha = Math.PI/n
    beta = Math.PI/m
    chc = 1.0/(Math.tan(alpha)*Math.tan(beta))
    len = Math.sqrt(chc*chc-1.0)/(1.0+chc)
    return len

regularTiling = (n,m,stop)->

    reflectSide = (poly,vs,i)->
        {c: c, v: v, t: t} = poly
        p = vs[i]
        q = vs[(i+1)%n]
        nt = Isometry.reflect(p,q)
        { c: reflect(p,q,c), v: p, t: nt.compose(t) }

    len = regularTilingRadius(n,m)
    tovisit = [ { c: Complex.zero, v: new Complex(len,0.0), t: new Isometry() } ]

    visited = d3.set()
    polygons = []

    while tovisit.length>0
        p = tovisit.shift()
        if !visited.has(p.c)
            visited.add(p.c)
            if !stop(p)
                polygons.push(p)
                vs = regularPolygon(p.c,p.v,n)
                for i in [0...n]
                    tovisit.push(reflectSide(p,vs,i))

    return polygons


delaunayEmbedding = (root) ->
    delaunayEmbeddingRec(root,0,new Isometry())

delaunayEmbeddingRec = (node, alpha, m) ->
    zero = Complex.zero
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

simulatedAnnealing = (s, kmax, driver) ->
    k = 0
    Emin = E = driver.energy(s)
    driver.save(s)
    saIter =  ->
        t = k/kmax
        if t>=1.0
            return 1.0
        k = k+1
        T = driver.temperature(t)
        [dE, c] = driver.change(s)
        if driver.accept(dE,T)>Math.random()
            driver.move(c,s)
            E = E+dE
            if E<Emin
                Emin = E
                driver.save(s)
        return t

root = self.hyperbolic ?= {}
root.Complex = Complex
root.Isometry = Isometry
root.segmentSVG = segmentSVG
root.circleSVG = circleSVG
root.reflect = reflect
root.circle = circle
root.perpBisectorSVG = perpBisectorSVG
root.perpBisectorHalfPlaneSVG = perpBisectorHalfPlaneSVG
root.randomPoint = randomPoint
root.regularPolygon = regularPolygon
root.regularPolygonSVG = regularPolygonSVG
root.polygonSVG = polygonSVG
root.distance = distance
root.regularTiling = regularTiling
root.regularTilingRadius = regularTilingRadius
root.delaunayEmbedding = delaunayEmbedding
root.simulatedAnnealing = simulatedAnnealing
Math.asinh = asinh
