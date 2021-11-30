// Unit: cm
dx = 100.0;
dy = 100.0;
dz = 100.0;

// Plate size
Plate_length = 10.0;
Plate_width = 10.0;

// Fine region
Fine_length = 2.04;
Fine_width = 2.04;

// Number of meshes
Num_fine = 34;
Num_trans = 16;
Param_trans = 1.1;

// z-coordinate of the top plane
Z_top = -0.001;

// Points in the top plane
Point(1) = {dx,                  dy,                 Z_top+dz, 1.0};
Point(2) = {dx,                  Fine_width/2.0+dy,  Z_top+dz, 1.0};
Point(3) = {Fine_length/2.0+dx,  Fine_width/2.0+dy,  Z_top+dz, 1.0};
Point(4) = {Fine_length/2.0+dx,  dy,                 Z_top+dz, 1.0};
Point(5) = {dx,                  Plate_width/2.0+dy, Z_top+dz, 1.0};
Point(6) = {Plate_length/2.0+dx, Plate_width/2.0+dy, Z_top+dz, 1.0};
Point(7) = {Plate_length/2.0+dx, dy,                 Z_top+dz, 1.0};

// Lines in the top plane
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {2, 5};
Line(6) = {3, 6};
Line(7) = {4, 7};
Line(8) = {5, 6};
Line(9) = {6, 7};

// Planes in the top plane
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Curve Loop(2) = {5, 8, -6, -2};
Plane Surface(2) = {2};
Curve Loop(3) = {6, 9, -7, -3};
Plane Surface(3) = {3};

// z-coordinate of the middle plane
Z_mid = -Plate_thick1 - 0.001;

// Points in the middle plane +7
Point(8)  = {dx,                  dy,                 Z_mid+dz, 1.0};
Point(9)  = {dx,                  Fine_width/2.0+dy,  Z_mid+dz, 1.0};
Point(10) = {Fine_length/2.0+dx,  Fine_width/2.0+dy,  Z_mid+dz, 1.0};
Point(11) = {Fine_length/2.0+dx,  dy,                 Z_mid+dz, 1.0};
Point(12) = {dx,                  Plate_width/2.0+dy, Z_mid+dz, 1.0};
Point(13) = {Plate_length/2.0+dx, Plate_width/2.0+dy, Z_mid+dz, 1.0};
Point(14) = {Plate_length/2.0+dx, dy,                 Z_mid+dz, 1.0};

// Lines in the middle plane +9
Line(10) = {8, 9};
Line(11) = {9, 10};
Line(12) = {10, 11};
Line(13) = {11, 8};
Line(14) = {9, 12};
Line(15) = {10, 13};
Line(16) = {11, 14};
Line(17) = {12, 13};
Line(18) = {13, 14};

// Planes in the middle plane +3
Curve Loop(4) = {10, 11, 12, 13};
Plane Surface(4) = {4};
Curve Loop(5) = {14, 17, -15, -11};
Plane Surface(5) = {5};
Curve Loop(6) = {15, 18, -16, -12};
Plane Surface(6) = {6};

// z-coordinate of the bottom plane
Z_bot = -Plate_thick2 - Plate_thick1 - 0.001;

// Points in the bottom plane +7
Point(15) = {dx,                  dy,                 Z_bot+dz, 1.0};
Point(16) = {dx,                  Fine_width/2.0+dy,  Z_bot+dz, 1.0};
Point(17) = {Fine_length/2.0+dx,  Fine_width/2.0+dy,  Z_bot+dz, 1.0};
Point(18) = {Fine_length/2.0+dx,  dy,                 Z_bot+dz, 1.0};
Point(19) = {dx,                  Plate_width/2.0+dy, Z_bot+dz, 1.0};
Point(20) = {Plate_length/2.0+dx, Plate_width/2.0+dy, Z_bot+dz, 1.0};
Point(21) = {Plate_length/2.0+dx, dy,                 Z_bot+dz, 1.0};

// Lines in the bottom plane +9
Line(19) = {15, 16};
Line(20) = {16, 17};
Line(21) = {17, 18};
Line(22) = {18, 15};
Line(23) = {16, 19};
Line(24) = {17, 20};
Line(25) = {18, 21};
Line(26) = {19, 20};
Line(27) = {20, 21};

// Planes in the bottom plane +3
Curve Loop(7) = {19, 20, 21, 22};
Plane Surface(7) = {7};
Curve Loop(8) = {23, 26, -24, -20};
Plane Surface(8) = {8};
Curve Loop(9) = {24, 27, -25, -21};
Plane Surface(9) = {9};

// Lines between different planes
Line(28) = {5, 12};
Line(29) = {2, 9};
Line(30) = {1, 8};
Line(31) = {4, 11};
Line(32) = {7, 14};
Line(33) = {6, 13};
Line(34) = {3, 10};
Line(35) = {12, 19};
Line(36) = {9, 16};
Line(37) = {8, 15};
Line(38) = {11, 18};
Line(39) = {14, 21};
Line(40) = {13, 20};
Line(41) = {10, 17};

// Vertical planes
Curve Loop(10) = {28, -14, -29, 5};
Plane Surface(10) = {10};
Curve Loop(11) = {29, -10, -30, 1};
Plane Surface(11) = {11};
Curve Loop(12) = {14, 35, -23, -36};
Plane Surface(12) = {12};
Curve Loop(13) = {10, 36, -19, -37};
Plane Surface(13) = {13};
Curve Loop(14) = {4, 30, -13, -31};
Plane Surface(14) = {14};
Curve Loop(15) = {13, 37, -22, -38};
Plane Surface(15) = {15};
Curve Loop(16) = {7, 32, -16, -31};
Plane Surface(16) = {16};
Curve Loop(17) = {16, 39, -25, -38};
Plane Surface(17) = {17};
Curve Loop(18) = {32, -18, -33, 9};
Plane Surface(18) = {18};
Curve Loop(19) = {18, 39, -27, -40};
Plane Surface(19) = {19};
Curve Loop(20) = {8, 33, -17, -28};
Plane Surface(20) = {20};
Curve Loop(21) = {17, 40, -26, -35};
Plane Surface(21) = {21};
Curve Loop(22) = {2, 34, -11, -29};
Plane Surface(22) = {22};
Curve Loop(23) = {11, 41, -20, -36};
Plane Surface(23) = {23};
Curve Loop(24) = {3, 31, -12, -34};
Plane Surface(24) = {24};
Curve Loop(25) = {12, 38, -21, -41};
Plane Surface(25) = {25};
Curve Loop(26) = {6, 33, -15, -34};
Plane Surface(26) = {26};
Curve Loop(27) = {15, 40, -24, -41};
Plane Surface(27) = {27};

// Volumes
Surface Loop(1) = {3, 18, 16, 26, 24, 6};
Volume(1) = {1};
Surface Loop(2) = {20, 2, 10, 5, 26, 22};
Volume(2) = {2};
Surface Loop(3) = {1, 11, 14, 22, 24, 4};
Volume(3) = {3};
Surface Loop(4) = {5, 8, 12, 21, 27, 23};
Volume(4) = {4};
Surface Loop(5) = {9, 19, 17, 27, 25, 6};
Volume(5) = {5};
Surface Loop(6) = {7, 13, 15, 25, 23, 4};
Volume(6) = {6};

Physical Surface(1) = {10, 11, 12, 13};
Physical Surface(2) = {14, 15, 16, 17};

Physical Volume(1) = {4, 5, 6};
Physical Volume(2) = {1, 2, 3};

// Number of nodes on each line, for generating meshes
Transfinite Curve {1, 2, 3, 4, 10, 11, 12, 13, 19, 20, 21, 22} = Num_fine + 1 Using Progression 1;
Transfinite Curve {8, 9, 17, 18, 26, 27} = Num_fine + 1 Using Progression 1;
Transfinite Curve {5, 6, 7, 14, 15, 16, 23, 24, 25} = Num_trans + 1 Using Progression Param_trans;
Transfinite Curve {32, 31, 30, 29, 34, 28, 33} = Num_mesh1 + 1 Using Progression 1;
Transfinite Curve {39, 38, 37, 36, 41, 35, 40} = Num_mesh2 + 1 Using Progression 1;

Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";

Mesh.SaveGroupsOfNodes = 1;
Mesh.SaveAll = 32;
