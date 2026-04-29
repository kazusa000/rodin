SetFactory("OpenCASCADE");

Mesh.MeshSizeMin = 0.045;
Mesh.MeshSizeMax = 0.09;

Box(1) = {-0.5, -0.5, -0.5, 2, 2, 2};

// Rounded many-lobed seed close to Fig. (a4): a compact center
// with seven smooth outer bulbs.
Sphere(2) = {0.50, 0.50, 0.50, 0.095};
Sphere(3) = {0.50, 0.69, 0.50, 0.095};
Sphere(4) = {0.65, 0.62, 0.50, 0.095};
Sphere(5) = {0.69, 0.47, 0.50, 0.095};
Sphere(6) = {0.60, 0.33, 0.50, 0.095};
Sphere(7) = {0.43, 0.30, 0.50, 0.095};
Sphere(8) = {0.31, 0.41, 0.50, 0.095};
Sphere(9) = {0.34, 0.59, 0.50, 0.095};

ov() = BooleanUnion{ Volume{2}; Delete; }{ Volume{3,4,5,6,7,8,9}; Delete; };
v() = BooleanFragments{ Volume{1}; Delete; }{ Volume{ov()}; Delete; };

eps = 1e-6;
obstacle[] = Volume In BoundingBox{0.18, 0.16, 0.38, 0.82, 0.84, 0.62};
fluid[] = Volume{:};
fluid() -= obstacle();

interface[] = Surface In BoundingBox{0.18-eps, 0.16-eps, 0.38-eps, 0.82+eps, 0.84+eps, 0.62+eps};
outer[] = Surface In BoundingBox{-0.5-eps, -0.5-eps, -0.5-eps, 1.5+eps, 1.5+eps, 1.5+eps};
outer() -= interface();

Physical Volume(2) = {obstacle()};
Physical Volume(3) = {fluid()};
Physical Surface(13) = {interface()};
Physical Surface(7) = {outer()};
