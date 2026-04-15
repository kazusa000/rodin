SetFactory("OpenCASCADE");

Mesh.MeshSizeMin = 0.055;
Mesh.MeshSizeMax = 0.11;

Box(1) = {0, 0, 0, 1, 1, 1};
Sphere(2) = {0.52, 0.5, 0.5, 0.11};
Sphere(3) = {0.62, 0.5, 0.5, 0.14};
Sphere(4) = {0.45, 0.60, 0.5, 0.14};
Sphere(5) = {0.45, 0.40, 0.5, 0.14};
ov() = BooleanUnion{ Volume{2}; Delete; }{ Volume{3,4,5}; Delete; };

v() = BooleanFragments{ Volume{1}; Delete; }{ Volume{ov()}; Delete; };

eps = 1e-6;
obstacle[] = Volume In BoundingBox{0.28, 0.24, 0.35, 0.77, 0.76, 0.65};
fluid[] = Volume{:};
fluid() -= obstacle();

interface[] = Surface In BoundingBox{0.28-eps, 0.24-eps, 0.35-eps, 0.77+eps, 0.76+eps, 0.65+eps};
outer[] = Surface In BoundingBox{-eps, -eps, -eps, 1+eps, 1+eps, 1+eps};
outer() -= interface();

Physical Volume(2) = {obstacle()};
Physical Volume(3) = {fluid()};
Physical Surface(13) = {interface()};
Physical Surface(7) = {outer()};
