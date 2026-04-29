SetFactory("OpenCASCADE");

Mesh.MeshSizeMin = 0.06;
Mesh.MeshSizeMax = 0.12;

Box(1) = {-0.5, -0.5, -0.5, 2, 2, 2};
Sphere(2) = {0.42, 0.5, 0.5, 0.18};
Sphere(3) = {0.58, 0.5, 0.5, 0.18};
ov() = BooleanUnion{ Volume{2}; Delete; }{ Volume{3}; Delete; };

v() = BooleanFragments{ Volume{1}; Delete; }{ Volume{ov()}; Delete; };

eps = 1e-6;
obstacle[] = Volume In BoundingBox{0.24, 0.31, 0.31, 0.76, 0.69, 0.69};
fluid[] = Volume{:};
fluid() -= obstacle();

interface[] = Surface In BoundingBox{0.24-eps, 0.31-eps, 0.31-eps, 0.76+eps, 0.69+eps, 0.69+eps};
outer[] = Surface In BoundingBox{-0.5-eps, -0.5-eps, -0.5-eps, 1.5+eps, 1.5+eps, 1.5+eps};
outer() -= interface();

Physical Volume(2) = {obstacle()};
Physical Volume(3) = {fluid()};
Physical Surface(13) = {interface()};
Physical Surface(7) = {outer()};
