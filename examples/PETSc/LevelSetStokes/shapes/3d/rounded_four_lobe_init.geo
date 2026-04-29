SetFactory("OpenCASCADE");

Mesh.MeshSizeMin = 0.05;
Mesh.MeshSizeMax = 0.10;

Box(1) = {-0.5, -0.5, -0.5, 2, 2, 2};

// Rounded four-lobed seed close to Fig. (a3): a smooth core with
// four overlapping bulbs arranged along the cardinal directions.
Sphere(2) = {0.50, 0.50, 0.50, 0.11};
Sphere(3) = {0.50, 0.66, 0.50, 0.12};
Sphere(4) = {0.50, 0.34, 0.50, 0.12};
Sphere(5) = {0.66, 0.50, 0.50, 0.12};
Sphere(6) = {0.34, 0.50, 0.50, 0.12};

ov() = BooleanUnion{ Volume{2}; Delete; }{ Volume{3,4,5,6}; Delete; };
v() = BooleanFragments{ Volume{1}; Delete; }{ Volume{ov()}; Delete; };

eps = 1e-6;
obstacle[] = Volume In BoundingBox{0.20, 0.20, 0.36, 0.80, 0.80, 0.64};
fluid[] = Volume{:};
fluid() -= obstacle();

interface[] = Surface In BoundingBox{0.20-eps, 0.20-eps, 0.36-eps, 0.80+eps, 0.80+eps, 0.64+eps};
outer[] = Surface In BoundingBox{-0.5-eps, -0.5-eps, -0.5-eps, 1.5+eps, 1.5+eps, 1.5+eps};
outer() -= interface();

Physical Volume(2) = {obstacle()};
Physical Volume(3) = {fluid()};
Physical Surface(13) = {interface()};
Physical Surface(7) = {outer()};
