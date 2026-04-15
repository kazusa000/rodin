SetFactory("OpenCASCADE");

Mesh.MeshSizeMin = 0.015;
Mesh.MeshSizeMax = 0.04;

Rectangle(1) = {0, 0, 0, 1, 1};
Disk(2) = {0.5, 0.5, 0, 0.18, 0.18};
Disk(3) = {0.5, 0.5, 0, 0.08, 0.08};

ring() = BooleanDifference{ Surface{2}; Delete; }{ Surface{3}; Delete; };
v() = BooleanFragments{ Surface{1}; Delete; }{ Surface{ring()}; Delete; };

eps = 1e-6;
obstacle[] = Surface In BoundingBox{0.31, 0.31, -eps, 0.69, 0.69, eps};
fluid[] = Surface{:};
fluid() -= obstacle();

interface[] = Curve In BoundingBox{0.31 - eps, 0.31 - eps, -eps, 0.69 + eps, 0.69 + eps, eps};
outer[] = Curve In BoundingBox{-eps, -eps, -eps, 1 + eps, 1 + eps, eps};
outer() -= interface();

Physical Surface(2) = {obstacle()};
Physical Surface(3) = {fluid()};
Physical Curve(13) = {interface()};
Physical Curve(7) = {outer()};
