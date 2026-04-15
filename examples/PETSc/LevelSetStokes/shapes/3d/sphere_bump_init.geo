SetFactory("OpenCASCADE");

Mesh.MeshSizeMin = 0.06;
Mesh.MeshSizeMax = 0.12;

Box(1) = {0, 0, 0, 1, 1, 1};
Sphere(2) = {0.5, 0.5, 0.5, 0.19};
Sphere(3) = {0.62, 0.54, 0.5, 0.05};
ov() = BooleanUnion{ Volume{2}; Delete; }{ Volume{3}; Delete; };

v() = BooleanFragments{ Volume{1}; Delete; }{ Volume{ov()}; Delete; };

eps = 1e-6;
obstacle[] = Volume In BoundingBox{0.30, 0.30, 0.30, 0.70, 0.70, 0.70};
fluid[] = Volume{:};
fluid() -= obstacle();

interface[] = Surface In BoundingBox{0.30-eps, 0.30-eps, 0.30-eps, 0.70+eps, 0.70+eps, 0.70+eps};
outer[] = Surface In BoundingBox{-eps, -eps, -eps, 1+eps, 1+eps, 1+eps};
outer() -= interface();

Physical Volume(2) = {obstacle()};
Physical Volume(3) = {fluid()};
Physical Surface(13) = {interface()};
Physical Surface(7) = {outer()};
