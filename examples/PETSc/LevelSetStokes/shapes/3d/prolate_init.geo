SetFactory("OpenCASCADE");

Mesh.MeshSizeMin = 0.06;
Mesh.MeshSizeMax = 0.12;

Box(1) = {0, 0, 0, 1, 1, 1};
Sphere(2) = {0.5, 0.5, 0.5, 0.22};
Dilate {{0.5, 0.5, 0.5}, {0.85, 0.85, 1.35}} { Volume{2}; }

v() = BooleanFragments{ Volume{1}; Delete; }{ Volume{2}; Delete; };

eps = 1e-6;
obstacle[] = Volume In BoundingBox{0.31, 0.31, 0.20, 0.69, 0.69, 0.80};
fluid[] = Volume{:};
fluid() -= obstacle();

interface[] = Surface In BoundingBox{0.31-eps, 0.31-eps, 0.20-eps, 0.69+eps, 0.69+eps, 0.80+eps};
outer[] = Surface In BoundingBox{-eps, -eps, -eps, 1+eps, 1+eps, 1+eps};
outer() -= interface();

Physical Volume(2) = {obstacle()};
Physical Volume(3) = {fluid()};
Physical Surface(13) = {interface()};
Physical Surface(7) = {outer()};

