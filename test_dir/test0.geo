cl = 0.015625;
Point(1) = { 0, 0, 0, cl };
Point(2) = { 1, 0, 0, cl };
Point(3) = { 1, 1, 0, cl };
Point(4) = { 0, 1, 0, cl };

Line(1) = { 1, 2 };
Line(2) = { 2, 3 };
Line(3) = { 3, 4 };
Line(4) = { 4, 1 };

Line Loop(5) = { 1, 2, 3, 4 };

Plane Surface(1) = { 5 };
Transfinite Surface { 1 };

Physical Surface(1) = { 1 };

