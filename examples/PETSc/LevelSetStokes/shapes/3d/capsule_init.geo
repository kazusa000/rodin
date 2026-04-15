SetFactory("OpenCASCADE");

Mesh.MeshSizeMin = 0.06;
Mesh.MeshSizeMax = 0.12;

Box(1) = {0, 0, 0, 1, 1, 1};
Cylinder(2) = {0.34, 0.5, 0.5, 0.32, 0, 0, 0.12, 2*Pi};
Sphere(3) = {0.34, 0.5, 0.5, 0.12};
Sphere(4) = {0.66, 0.5, 0.5, 0.12};
ov() = BooleanUnion{ Volume{2}; Delete; }{ Volume{3,4}; Delete; };

v() = BooleanFragments{ Volume{1}; Delete; }{ Volume{ov()}; Delete; };

eps = 1e-6;
obstacle[] = Volume In BoundingBox{0.22, 0.37, 0.37, 0.78, 0.63, 0.63};
fluid[] = Volume{:};
fluid() -= obstacle();

interface[] = Surface In BoundingBox{0.22-eps, 0.37-eps, 0.37-eps, 0.78+eps, 0.63+eps, 0.63+eps};
outer[] = Surface In BoundingBox{-eps, -eps, -eps, 1+eps, 1+eps, 1+eps};
outer() -= interface();

Physical Volume(2) = {obstacle()};
Physical Volume(3) = {fluid()};
Physical Surface(13) = {interface()};
Physical Surface(7) = {outer()};

