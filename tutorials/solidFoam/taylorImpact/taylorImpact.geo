// cc = Circle circumference
// cr = Circle radius
// cy = Cylinder axis

//cc = 2; cr = 2; cy = 40;  // 480 cells
//cc = 3; cr = 3; cy = 50;  // 1350 cells
cc = 5; cr = 5; cy = 100;   // 7500 cells

// Domain limits
r = 0.0032; 		// Radius of circle
d = 0.0014; 		// Square length
ymax = 0.0324;		// Bar height
p = r*Sqrt(2)/2; 	// Coordinate of point on circle

// Nodes
nc = cc + 1;
nr = cr + 1;
ny = cy + 1;


// Square section
Point(1) = {0, 0, 0};
Point(2) = {d, 0, 0};
Point(3) = {d, 0, d};
Point(4) = {0, 0, d};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Transfinite Line{1,2,3,4} = nc;
Line Loop(1) = {1,2,3,4} ;
Plane Surface(1) = {1};
Transfinite Surface{1} = {1,2,3,4};
Recombine Surface{1};

// Circular section 1
Point(5) = {r, 0, 0};
Point(6) = {p, 0, p};
Line(5) = {2,5};
Circle(6) = {5,1,6};
Line(7) = {6,3};
Transfinite Line{6} = nc;
Transfinite Line{5,7} = nr;
Line Loop(2) = {5,6,7,-2};
Plane Surface(2) = {2};
Transfinite Surface{2} = {2,5,6,3};
Recombine Surface{2};

// Circular section 2
Point(8) = {0, 0, r};
Circle(8) = {6,1,8};
Line(9) = {8,4};
Transfinite Line{8} = nc;
Transfinite Line{9} = nr;
Line Loop(3) = {-7,8,9,-3};
Plane Surface(3) = {3};
Transfinite Surface{3} = {3,6,8,4};
Recombine Surface{3};

// Extrusion
Extrude{0,ymax,0}{
Surface{1,2,3};
Layers{cy};Recombine;
}

// Definition of surfaces for boundary conditions
Physical Surface("free") 	   = {31,53,75,44,66};
Physical Surface("symmetricY") = {1,2,3};
Physical Surface("symmetricX") = {30,70};
Physical Surface("symmetricZ") = {40,18};

// Definition of a volume
Physical Volume("volume") = {1,2,3};
