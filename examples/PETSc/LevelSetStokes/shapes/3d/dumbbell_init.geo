SetFactory("OpenCASCADE");

Mesh.MeshSizeMin = 0.055;
Mesh.MeshSizeMax = 0.11;

Box(1) = {-0.5, -0.5, -0.5, 2, 2, 2};
Cylinder(2) = {0.38, 0.5, 0.5, 0.24, 0, 0, 0.07, 2*Pi};
Sphere(3) = {0.38, 0.5, 0.5, 0.15};
Sphere(4) = {0.62, 0.5, 0.5, 0.15};
ov() = BooleanUnion{ Volume{2}; Delete; }{ Volume{3,4}; Delete; };

v() = BooleanFragments{ Volume{1}; Delete; }{ Volume{ov()}; Delete; };

eps = 1e-6;
obstacle[] = Volume In BoundingBox{0.22, 0.34, 0.34, 0.78, 0.66, 0.66};
fluid[] = Volume{:};
fluid() -= obstacle();

interface[] = Surface In BoundingBox{0.22-eps, 0.34-eps, 0.34-eps, 0.78+eps, 0.66+eps, 0.66+eps};
outer[] = Surface In BoundingBox{-0.5-eps, -0.5-eps, -0.5-eps, 1.5+eps, 1.5+eps, 1.5+eps};
outer() -= interface();

Physical Volume(2) = {obstacle()};
Physical Volume(3) = {fluid()};
Physical Surface(13) = {interface()};
Physical Surface(7) = {outer()};
