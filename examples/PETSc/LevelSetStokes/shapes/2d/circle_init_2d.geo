SetFactory("OpenCASCADE");

Mesh.MeshSizeMin = 0.02;
Mesh.MeshSizeMax = 0.05;

Rectangle(1) = {0, 0, 0, 1, 1};
Disk(2) = {0.5, 0.5, 0, 0.12, 0.12};

v() = BooleanFragments{ Surface{1}; Delete; }{ Surface{2}; Delete; };

eps = 1e-6;
obstacle[] = Surface In BoundingBox{0.38, 0.38, -eps, 0.62, 0.62, eps};
fluid[] = Surface{:};
fluid() -= obstacle();

interface[] = Curve In BoundingBox{0.38 - eps, 0.38 - eps, -eps, 0.62 + eps, 0.62 + eps, eps};
outer[] = Curve In BoundingBox{-eps, -eps, -eps, 1 + eps, 1 + eps, eps};
outer() -= interface();

Physical Surface(2) = {obstacle()};
Physical Surface(3) = {fluid()};
Physical Curve(13) = {interface()};
Physical Curve(7) = {outer()};
