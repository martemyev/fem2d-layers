cl = 0.1;
Point(1) = { 0, 0, 0, cl };
Point(2) = { 1, 0, 0, cl };
Point(3) = { 1, 1, 0, cl };
Point(4) = { 0, 1, 0, cl };

Point(5) = { 0.3, 0.8, 0, cl };
Point(6) = { 0.4, 0.8, 0, cl };
Point(7) = { 0.3, 0.9, 0, cl };
Point(8) = { 0.2, 0.8, 0, cl };
Point(9) = { 0.3, 0.7, 0, cl };

Line(1) = { 1, 2 };
Line(2) = { 2, 3 };
Line(3) = { 3, 4 };
Line(4) = { 4, 1 };

Circle(6) = { 7, 5, 6 };
Circle(7) = { 7, 5, 8 };
Circle(8) = { 9, 5, 6 };
Circle(9) = { 9, 5, 8 };

Line Loop(5) = { 1, 2, 3, 4 };
Line Loop(10) = { 7, -9, 8, -6 };

Plane Surface(11) = {10};
Physical Surface(11) = { 11 }; // circle

Plane Surface(1) = { 5, 10 };
Physical Surface(1) = { 1 }; // square without circle

