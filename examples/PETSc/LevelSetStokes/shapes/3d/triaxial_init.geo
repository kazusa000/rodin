SetFactory("OpenCASCADE");

Mesh.MeshSizeMin = 0.06;
Mesh.MeshSizeMax = 0.12;

Box(1) = {0, 0, 0, 1, 1, 1};
Sphere(2) = {0.5, 0.5, 0.5, 0.22};
Dilate {{0.5, 0.5, 0.5}, {1.18, 0.92, 0.72}} { Volume{2}; }

v() = BooleanFragments{ Volume{1}; Delete; }{ Volume{2}; Delete; };

eps = 1e-6;
obstacle[] = Volume In BoundingBox{0.23, 0.29, 0.33, 0.77, 0.71, 0.67};
fluid[] = Volume{:};
fluid() -= obstacle();

interface[] = Surface In BoundingBox{0.23-eps, 0.29-eps, 0.33-eps, 0.77+eps, 0.71+eps, 0.67+eps};
outer[] = Surface In BoundingBox{-eps, -eps, -eps, 1+eps, 1+eps, 1+eps};
outer() -= interface();

Physical Volume(2) = {obstacle()};
Physical Volume(3) = {fluid()};
Physical Surface(13) = {interface()};
Physical Surface(7) = {outer()};
