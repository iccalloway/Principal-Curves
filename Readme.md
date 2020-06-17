Projecting onto tongue surface data with principal curves
================
Ian Calloway
June 2, 2020

  - [What are Principal Curves?](#what-are-principal-curves)
  - [Pakcage Info](#pakcage-info)
  - [Defining Coordinates](#defining-coordinates)
      - [PCA](#pca)
      - [Principal Curves](#principal-curves)
  - [Tongue Surface Data](#tongue-surface-data)

<style>
body {
text-align: justify}
</style>

``` r
knitr::opts_chunk<img src="https://rawgit.com/iccalloway/Principal-Curves/None/svgs/20e0cbdbc877259346b2c7169abba350.svg?invert_in_darkmode" align=middle width=1531.0806670499999pt height=2436.7478057999997pt/>lambda,
      c=orth_proj(rbind(as.matrix(x),as.matrix(y)), test_proj$s, test_proj$ord)
      )
    )[1:nrow(x),]
  
  ylam <- range(test_proj<img src="https://rawgit.com/iccalloway/Principal-Curves/None/svgs/56a44f4e8ac2b5f2290edf3104b5d637.svg?invert_in_darkmode" align=middle width=270.68596499999995pt height=24.65753399999998pt/>lambda)],na.rm=T)

  ## Option: remove points projected to endpoints
  if(endpoint){ 
    ind <- which(x_f$lambda <= (ylam[1] + tol) | x_f$lambda >= (ylam[2]-tol))
    if(length(ind)>0){
      x_f[ind,]$c <- NA
    }
  }

  ## Option: remove points that are uniformly projecting to a single point on the curve
  if(trim){
    excess <- as.numeric(names(table(x_f$lambda))[table(x_f$lambda) > threshold*nrow(x_f)])
    if(length(excess)>0){
      ind <- sapply(x_f$lambda, function(x) any(abs(x - excess) < tol))
      if(length(ind)>0){
        x_f[ind,]$c <- NA  
      }
    }
  }
  return(x_f)
}
```

-----

## Defining Coordinates

### PCA

It is possible to project the data onto its principal components. For
example, for a vector \(\vec{u} \in \mathbb{R}^2\) in the data set and
**PC1** as the first principal component:

\[\mathbf{P}_{PC1}(\vec{u})=\vec{u}_{PC1}=\lambda\vec{w}_1\]
\[\vec{u} - \mathbf{P}_{PC1}(\vec{u})=\vec{u}_{PC1^\perp}=c\vec{w}_2\]
where \(\mathbf{P}_{PC1}\) is a projection matrix onto the first
principal component, \(\lambda,c \in \mathbb{R}\) are scalars,
\(\vec{w}_1\) is a unit length vector parallel to the first principal
component, and \(\vec{w}_2\) is a unit length vector orthogonal to the
first principal component and parallel to \(\vec{y}\). \(\vec{w}_1\) and
\(\vec{w}_2\) add up to \(\vec{u}\) and are orthogonal to one another.

\(W=\{\vec{w}_1, \vec{w}_2\}\) form an orthogonal basis in
\(\mathbb{R}^2\), and the coordinate vector \([\vec{u}]_W=(\lambda,c)\).

#### Plot Generating Code

The code below can be used to help get a feel for what effect projection
onto an arbitrary principal component (or principal curve) has on the
data. It takes as input:

1.  A range (the same range is used for both x and y axes)

2.  A sampling step

3.  An arbitrary function (\(\mathbb{R} \mapsto \mathbb{R}\))

It returns three plots:

1.  A grid displaying the function and each arrows from each point to
    its projection on the principal curve

2.  A grid displaying the finite gradient magnitude for each point on
    the grid

3.  A plot of each point on the grid in the coordinate vector
    \((\lambda, c)\).

<!-- end list -->

``` r
example1 <- function(range,step, mapping){
  tol <- 1e-10
  
  grange <- seq(range[1],range[2],step)
  grange_narrow <- seq(range[1],range[2],step/5)
  
  grid <- expand.grid(grange,grange)
  grid_narrow <- expand.grid(grange_narrow,grange_narrow)
  
  x <- grange_narrow
  data <- as.matrix(cbind(x,y=mapping(x)))
  projections <- project_to_curve(as.matrix(grid),as.matrix(data))
  grid<img src="https://rawgit.com/iccalloway/Principal-Curves/None/svgs/ecc508266e26c840ca493a0ed6c4e11b.svg?invert_in_darkmode" align=middle width=203.80215239999998pt height=24.65753399999998pt/>s - grid,1, function(x){
    if (abs(x[1]) < tol & abs(x[2]) < tol){
      return(NA)
    } else {
      return(atan2(x[2],x[1]))
    }
  })
  data_f<-as.data.frame(data)
  
  
  p1<-ggplot() +
  geom_spoke(data=grid, aes(x=Var1, y=Var2, angle=angle, radius=0.4/2, colour=Var2), arrow = arrow(length = unit(.05, 'inches')),size=2)+
  geom_line(data=data_f[data_f<img src="https://rawgit.com/iccalloway/Principal-Curves/None/svgs/e693e7b7712e40e7c8095b77087f7173.svg?invert_in_darkmode" align=middle width=142.78764554999998pt height=24.65753399999998pt/>y <= range[2],], aes(x=x,y=y),colour="black",size=1) +
  scale_colour_viridis()+
  labs(x="x",y="y", colour="y")
  
  projections_narrow <- project_to_curve(as.matrix(grid_narrow),as.matrix(data))

  
  diffs <- cbind(grid_narrow,projections_narrow<img src="https://rawgit.com/iccalloway/Principal-Curves/None/svgs/2d867fa6d95764e8ead478bc3436deb5.svg?invert_in_darkmode" align=middle width=596.2681609499999pt height=47.671232400000015pt/>"Var1.1")
  ys_m <- daply(ys, .(Var1,Var2), function(x) x<img src="https://rawgit.com/iccalloway/Principal-Curves/None/svgs/5bbd98cb5bae80ab9eeb28fab858d9f8.svg?invert_in_darkmode" align=middle width=851.2126507499999pt height=325.02278490000003pt/>y >= range[1] & data_f<img src="https://rawgit.com/iccalloway/Principal-Curves/None/svgs/a82893adf831b9df93e3c5b8ebf1cb74.svg?invert_in_darkmode" align=middle width=700.2744424499999pt height=85.29680940000001pt/>type <- "grid"
  
  projections <- grid[,c(1,2,4)]
  
  projections_f <- cbind(projections,proj_simple(as.matrix(projections[,c(1,2)]),as.matrix(data)))
   
   p3 <- ggplot() +
     geom_point(data=projections_f[projections_f$type == 'grid',], aes(x=lambda, y=c, colour=Var2, group=Var1),size=2)+
     geom_text(data=projections_f[projections_f$type == 'grid',], aes(x=lambda, y=c, colour=Var2, label=Var1),size=4,vjust=2)+
     scale_colour_viridis(option="viridis")+
     labs(x="\u03BB",y="c", title="Local Coordinates",colour="y")
   return(list(p1,p2,p3))
   
  
}
```

-----

#### Example 1: Projecting onto \(y=x\)

The figure below is a visualization of what it would look like if all
points were projected onto a principal component spanned by \([1,1]\).
For each point in the grid, the arrow is pointed in the direction of the
place on the principal curve that it would project to. All points below
the line point in the direction of \([-1,1]\), and those above the line
point in the exact opposite direction.

``` r
result <- example1(c(-5,5),0.5, function(x) x)
print(result[[1]])
```

![](PrincipalCurves_files/figure-gfm/example_pca_1-1.png)<!-- -->

If we use the coordinate vectors \((\lambda,c)\) to plot each point on
the grid, then we get the figure below.

  - The color of the point corresponds to its original y-value

  - The value below each point corresponds to its original x-value.

In this new coordinate system, the data have essentially been rotated 45
degrees clockwise.

``` r
print(result[[3]])
```

![](PrincipalCurves_files/figure-gfm/example_pca_2-1.png)<!-- -->

-----

### Principal Curves

There is not a guaranteed linear transformation for a projection onto a
principal curve. However, one can define *local coordinates* based on
projections onto the curve.

Let \(f: \Lambda \subset \mathbb{R}^1 \mapsto A\) be a smooth function
that maps real values to points on the principal curve -
\(f(\lambda \in \Lambda)=a \in A\).

As noted earlier, all \(\vec{u} \in \mathbb{R}^p\) can be projected onto
the principal curve.

\[\phi(\vec{u})= \sup_\lambda\{\lambda: ||\vec{u}- f(\lambda)|| = \inf_\mu ||\vec{u}-f(\mu)||\}\]
If there is only one point on the curve closest to \(\vec{u}\), then
that point is the projection. If there are multiple points, however,
then the point corresponding to the largest \(\lambda\) is chosen. While
this guarantees that every point in \(\mathbb{R}^2\) will project to the
principal curve, it will also create some headaches.

The difference of the vector and its projection onto the curve is:

\[\vec{u}-\phi(\vec{u}) = c\vec{w}_3\] Where \(\vec{w_3}\) is a unit
vector parallel to \(\vec{u}-\phi(\vec{u})\). While \(\vec{w_2}\) is the
same for different points along the principal component line, \(w_3\) is
not guaranteed to be the same for different points along the principal
curve.

The projection and the orthogonal projection still sum to \(\vec{u}\).
If \(f\) is differentiable around \(\lambda'\), where
\(f(\lambda')=\phi(\vec{u})\), then:

\[(\vec{u}-\phi(\vec{u})) \perp \nabla f(\lambda')\]

The line tangent to the curve where where \(\vec{u}\) projects onto the
principal curve is orthogonal to \(\vec{u} - \phi(\vec{u})\).

-----

#### A Coordinate Chart

We can use a principal curve to define local coordinates according to a
point’s position relative to the curve. Unfortunately, unlike PCA, we do
not have a guaranteed ability to describe every point in
\(\mathbb{R}^2\) this way. The next section describes \(\mathbf{M}\),
the subset of \(\mathbb{R}^2\) that would allow us to move between
standard coordinates and position relative to the curve.

Let \(\mathbf{M} \subset \mathbb{R}^2\) be a 2-manifold. A function to
map from \(\mathbb{R}^2\) to \(\mathbf{M}\) would be:

\[g(\vec{u})=[\lambda_f(\vec{u}),\vec{u}-\phi(\vec{u})]^T\]

where \(\lambda_f(\vec{u})\) is the value of \(\lambda\) where
\(\vec{u}\) projects to the principal curve.

For this mapping to useful, it should be bijective - it should be
possible to go from \(\mathbf{M}\) to \(\mathbb{R}^2\) as well as
backwards. If \(g\) is not, then the coordinate system will have some
undesirable properties. \(\mathbf{M}\) includes all the points in
\(\mathbb{R}^2\) where:

-----

**1. There is only one point on the principal curve closest to
\(\vec{u}\)**

Although every point in \(\mathbb{R}^2\) is projected to the principal
curve, points closest to multiple parts of a curve do not behave as
well.

If there are two points \(f(\lambda_1)\) and \(f(\lambda_2)\) on the
principal curve that are closest to \(\vec{u}\), and we assume
\(g^{-1}\) exists, then:

1.  \(\lambda_f(\vec{u})\) must equal equal \(\lambda_1\) or
    \(\lambda_2\) but not both. Assume it equals \(\lambda_1\).
2.  \(g^{-1}(\lambda_1,c)=g^{-1}(\lambda_2,c)=\vec{u}\) for some \(c\).
3.  \(g(g^{-1}(\lambda_2,c))=g(\vec{u})=[\lambda_f(u),u-\phi(u)^T]=[\lambda_1,c]\)

This is contradiction (\(f(f^{-1}(x))=x\)), so the inverse cannot exist
if there are multiple points on the principal curve equidistant to
\(\vec{u}\).

-----

**2. The point does not project to either endpoint of the curve.**

If we restrict the domain of \(g\) to points \(\vec{u}\) closest to only
one point on the principal curve, then we can define an inverse function
\(h=g^{-1}\) as follows:

\[h(\lambda, c)=f(\lambda) + c(T_{f(\lambda)}M)^\perp\]

where \((T_{f(\lambda)}M)^\perp\) is a unit vector orthogonal to \(f\)
at \(\lambda\).

This function requires the principal curve to be differentiable at
\(\lambda\), which means that we must restrict the domain of \(g\)
further to those regions that have that property. Principal curves are
not smooth at their endpoints, so these must be excluded. *(Although you
can get around this issue by extrapolating past the defined endpoints of
the curve. The princurve package does offer this functionality, but I
don’t discuss it here)*

If the point does not have a unique nearest point on the curve, and it
does not project to an endpoint, then the local coordinates
\((\lambda, c)\) should not behave in a wacky manner.

-----

#### Example 2: Projecting onto \(y=x^2\)

The figure below is a visualization of what it would look like if all
points were projected onto a principal curve that happens to defined as
\(y=x^2\). Once again, the arrows are pointed in the direction of the
place on the curve it would project to.

For the ray extending upward from \([0,0.5]\), there is a sudden change
in the direction of the points on the grid even though we are not
crossing the boundary of the curve. Over this region, each point is
equidistant from *two* points on the curve. The principal curve
algorithm is able to decide which of the two choices to project to, but
there do seem to be some random fluctuations in which side it prefers.

``` r
result <- example1(c(-5,5),0.5, function(x) x^2)
print(result[[1]])
```

![](PrincipalCurves_files/figure-gfm/example_squared_1-1.png)<!-- -->

We can track these sudden changes in projection direction. For the
vectors extending from each point to the point on the curve it projects
to, we can compute the finite gradient magnitude, shown in the plot
below. While the x and y axes are the same as before, the heat map
colors correspond to a differences in gradient magnitude moving
diagonally to the right and upwards along the plot surface. We see that
at or around \([0,0.5]\) upward (with some jitteriness apparently due to
variability in the projection function), there is a sudden shift in the
direction of projection.

``` r
print(result[[2]])
```

![](PrincipalCurves_files/figure-gfm/example_squared_2-1.png)<!-- -->

As discussed above, data that sit along this region or cross it may not
behave in the way we like. What does our grid look like when we project
it onto this curve?

``` r
print(result[[3]])
```

![](PrincipalCurves_files/figure-gfm/example_squared_3-1.png)<!-- -->

Unlike the first example, it is not possible to describe the
transformation to the data as a simple rotation. It seems to be
stretched in some fashion. A few items worth noting:

1.  If you identify the points that sit along \(c=0\), these are the
    same as the points that would lie along \(y=x^2\).

2.  At above \(y=0.5\), the grid seems to be split in two\! If you trace
    the line originally at \(y=5\) (the brightest yellow), for example,
    you will reach a discontinuity and the line will then continue
    elsewhere. This is one of the dangers of having data that lie on or
    cross regions with non-unique nearest principal curve points -
    crossing these regions results in a discontinuity in the local
    coordinate system.

3.  For larger negative y-values, the data points seem to be clustered
    together more closely. The line originally at \(y=-5\) (the deepest
    purple) appears to take the shape of a tight parabola.

If you plan to do analysis that involves this kind of projection onto a
principal curve, it would be useful to identify the problematic regions
beforehand.

-----

## Tongue Surface Data

For this demo, I am using ultrasound video frames of rhymes from /CʌC/,
/CʌlC/, and /CVoʊC/ words produced by a speaker of American English. The
figures below plot a single video frame in the speaker’s production of
“cult” (on the left) and “soak” (on the right).

  - Higher X values correspond to a position further forward in the
    vocal tract
  - Higher Y values correspond to a higher position in the vocal tract

In the speaker’s production of “cult”, there appear to be two tongue
constrictions - one at the tongue tip and another back in the tongue
body. In contrast, for their production of “soak”, the tongue dorsum
seems to be bunched.

``` r
merged <- fread('~/p4.csv')[,-c(1,2)] ## Load Data

ggplot(merged[frame %in% c(113,291)], aes(x=x_scaled, y=y_scaled)) +
  geom_point() +
  facet_wrap(~word) +
  labs(x="x",y="y", title="Tongue Surfaces")+
  coord_fixed(ratio=0.5)
```

![](PrincipalCurves_files/figure-gfm/plot_raw-1.png)<!-- -->

What if we wanted to project either of these tongue surface curves onto
the other? Analogous to the problematic region for \(y=x^2\), it is
possible for certain tongue shapes to have regions that result in
discontinuities if another tongue is projected onto it.

This second bit of code gives projection-related information for a grid
of points centered around an arbitrary curve.

``` r
example2 <- function(curve, step){
  tol <- 1e-10
  
  x_range <- range(curve[,1],na.rm=T)
  y_range <- range(curve[,2],na.rm=T)
  
  x_mid <- mean(x_range,na.rm=T)
  y_mid <- mean(y_range,na.rm=T)

  xspan <- seq(x_mid - 2*(x_mid-x_range[1]),x_mid+2*(x_range[2]-x_mid),step)
  yspan <- seq(y_mid - 2*(y_mid-y_range[1]),y_mid+2*(y_range[2]-y_mid),step)
  
  xspan_narrow <- seq(x_mid - 2*(x_mid-x_range[1]),x_mid+2*(x_range[2]-x_mid),step/8)
  yspan_narrow <- seq(y_mid - 2*(y_mid-y_range[1]),y_mid+2*(y_range[2]-y_mid),step/8)
  
  grid <- expand.grid(xspan,yspan)
  grid_narrow <- expand.grid(xspan_narrow,yspan_narrow)
  
  
  projections <- project_to_curve(as.matrix(grid),as.matrix(curve))
  grid<img src="https://rawgit.com/iccalloway/Principal-Curves/None/svgs/ecc508266e26c840ca493a0ed6c4e11b.svg?invert_in_darkmode" align=middle width=203.80215239999998pt height=24.65753399999998pt/>s - grid,1, function(x){
    if (abs(x[1]) < tol & abs(x[2]) < tol){
      return(NA)
    } else {
      return(atan2(x[2],x[1]))
    }
  })
  data_f<-as.data.frame(curve)
  
  
  p1<-ggplot() +
  geom_spoke(data=grid, aes(x=Var1, y=Var2, angle=angle, radius=step/3, colour=Var2), arrow = arrow(length = unit(.05, 'inches')),size=1)+
  geom_line(data=data_f, aes(x=x_scaled,y=y_scaled),colour="black",size=1) +
  scale_colour_viridis()+
  labs(x="x",y="y", colour="y")

  projections_narrow <- project_to_curve(as.matrix(grid_narrow),as.matrix(curve))
  
  diffs <- cbind(grid_narrow,projections_narrow<img src="https://rawgit.com/iccalloway/Principal-Curves/None/svgs/2d867fa6d95764e8ead478bc3436deb5.svg?invert_in_darkmode" align=middle width=596.2681609499999pt height=47.671232400000015pt/>"Var1.1")
  ys_m <- daply(ys, .(Var1,Var2), function(x) x<img src="https://rawgit.com/iccalloway/Principal-Curves/None/svgs/03a25acf4d057fd1fba360f8fa6a6b75.svg?invert_in_darkmode" align=middle width=851.2126507499999pt height=440.36530559999994pt/>type <- "grid"
  projections <- grid[,c(1,2,4)]
  
  projections_f <- cbind(projections,proj_simple(as.matrix(projections[,c(1,2)]),as.matrix(curve),T))
   
   p3 <- ggplot() +
     geom_point(data=projections_f[projections_f$type == 'grid',], aes(x=lambda, y=c, colour=Var2, group=Var1),size=2)+
     geom_text(data=projections_f[projections_f$type == 'grid',], aes(x=lambda, y=c, colour=Var2, label=round(Var1,2)),size=4,vjust=1.5)+
     scale_colour_viridis(option="viridis")+
     labs(x="\u03BB",y="c", title="Local Coordinates",colour="y")
   return(list(p1,p2,p3))
}
```

-----

Looking at the tongue surface for ‘soak’, we notice that there seems to
be a sudden change in projection direction directly below the tongue
dorsum but not really elsewhere.

``` r
soak <- merged[frame==291,c('x_scaled','y_scaled')]
soak_pc <- principal_curve(as.matrix(soak))
test <-project_to_curve(as.matrix(soak),soak_pc<img src="https://rawgit.com/iccalloway/Principal-Curves/None/svgs/6c906b64d1eb34d99085f3b559977e87.svg?invert_in_darkmode" align=middle width=59.90444294999998pt height=24.65753399999998pt/>ord,])
soak_srt <- soak[test<img src="https://rawgit.com/iccalloway/Principal-Curves/None/svgs/09f183b33e6ddbbde2c7a487c9a63f5a.svg?invert_in_darkmode" align=middle width=1827.23121075pt height=745.8447050999999pt/>s[cult_pc<img src="https://rawgit.com/iccalloway/Principal-Curves/None/svgs/372f5f38f1ffd9e0596bec654d181418.svg?invert_in_darkmode" align=middle width=185.37526755pt height=24.65753399999998pt/>ord,]
result <- example2(cult_srt,10)
print(result[[2]])
```

![](PrincipalCurves_files/figure-gfm/cult-1.png)<!-- -->

As you might guess, it would be very difficult to project a tongue
surface onto this particular curve without running into discontinuities
unless it already happened to lie close to the xcurve.

``` r
print(result[[3]])
```

![](PrincipalCurves_files/figure-gfm/cult_3-1.png)<!-- -->

-----

I’ve chosen a production that is suitable for the data set I’m working
with. This a vowel-initial ultrasound production of ‘suck’.

``` r
suck <- merged[word=='suck' & rep==1 & nucpoint==0,c('x_scaled','y_scaled')]
suck_pc <- principal_curve(as.matrix(suck))
test <-project_to_curve(as.matrix(suck),suck_pc<img src="https://rawgit.com/iccalloway/Principal-Curves/None/svgs/1fca65165aa5aacb0503d270470939ee.svg?invert_in_darkmode" align=middle width=59.771315999999985pt height=24.65753399999998pt/>ord,])
suck_srt <- suck[test$ord,]
result <- example2(suck_srt,10)
print(result[[2]])
```

![](PrincipalCurves_files/figure-gfm/suck-1.png)<!-- --> These are the
two frames from earlier projected onto this production of suck. A few
things to note:

1.  These new curves do not appear to have discontinuities\!

2.  The y-position seems to match with constriction degree relative to
    ‘suck’ at that location in the tongue. ‘Soak’ shows greater
    constriction at the tongue dorsum and not much constriction
    elsewhere, while ‘cult’ shows more constriction closer to the tongue
    root and relatively low constriction elsewhere.

The main utility of these projections will be to identify changes in
constriction degree and constriction location over time.

``` r
news <- cbind(merged[frame %in% c(113,291),],proj_simple(as.matrix(merged[frame %in% c(113,291),c('x_scaled','y_scaled')]),suck_srt))
ggplot(news, aes(x=200-lambda,y=c))+
geom_point() +
facet_wrap(~frame)+
labs(x="\u03BB", y="c")
```

![](PrincipalCurves_files/figure-gfm/tongue_project-1.png)<!-- -->

``` r
news <- proj_simple(as.matrix(merged[,c('x_scaled','y_scaled')]),suck_srt)
```

More sections to come soon…
