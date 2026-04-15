SetFactory("OpenCASCADE");

Mesh.MeshSizeMin = 0.05;
Mesh.MeshSizeMax = 0.09;

Box(1) = {0, 0, 0, 1, 1, 1};
Torus(2) = {0.5, 0.5, 0.5, 0.14, 0.08};

v() = BooleanFragments{ Volume{1}; Delete; }{ Volume{2}; Delete; };

eps = 1e-6;
obstacle[] = Volume In BoundingBox{0.28, 0.28, 0.41, 0.72, 0.72, 0.59};
fluid[] = Volume{:};
fluid() -= obstacle();

interface[] = Surface In BoundingBox{0.28-eps, 0.28-eps, 0.41-eps, 0.72+eps, 0.72+eps, 0.59+eps};
outer[] = Surface In BoundingBox{-eps, -eps, -eps, 1+eps, 1+eps, 1+eps};
outer() -= interface();

Physical Volume(2) = {obstacle()};
Physical Volume(3) = {fluid()};
Physical Surface(13) = {interface()};
Physical Surface(7) = {outer()};
