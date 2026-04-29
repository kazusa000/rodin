SetFactory("OpenCASCADE");

Mesh.MeshSizeMin = 0.06;
Mesh.MeshSizeMax = 0.12;

Box(1) = {-0.5, -0.5, -0.5, 2, 2, 2};
Sphere(2) = {0.5, 0.5, 0.5, 0.22};
Dilate {{0.5, 0.5, 0.5}, {1.10, 1.10, 0.55}} { Volume{2}; }

v() = BooleanFragments{ Volume{1}; Delete; }{ Volume{2}; Delete; };

eps = 1e-6;
obstacle[] = Volume In BoundingBox{0.25, 0.25, 0.38, 0.75, 0.75, 0.62};
fluid[] = Volume{:};
fluid() -= obstacle();

interface[] = Surface In BoundingBox{0.25-eps, 0.25-eps, 0.38-eps, 0.75+eps, 0.75+eps, 0.62+eps};
outer[] = Surface In BoundingBox{-0.5-eps, -0.5-eps, -0.5-eps, 1.5+eps, 1.5+eps, 1.5+eps};
outer() -= interface();

Physical Volume(2) = {obstacle()};
Physical Volume(3) = {fluid()};
Physical Surface(13) = {interface()};
Physical Surface(7) = {outer()};
