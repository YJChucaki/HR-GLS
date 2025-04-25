/// saddle_surface_fixed.geo
SetFactory("OpenCASCADE");

// 定义参数时使用变量避免硬编码
a = 2.0;
b = 2.0;

// 创建参数化曲面（自动分配标签）
surf = newv;
ParametricSurface(surf) = {
    {"u", "v", "u*v"},
    {-a, a, -b, b},
    {"", "", ""},
    {20, 20}
};

// 必须显式同步才能正确获取几何信息
Synchronize;

// 自动获取曲面边界
boundary[] = CombinedBoundary{ Surface{surf}; };

// 物理组使用变量引用
Physical Curve("BoundaryEdges") = boundary[];
Physical Surface("Saddle") = {surf};

// 设置网格参数
Mesh.CharacteristicLengthMin = 0.2;
Mesh.CharacteristicLengthMax = 0.2;

// 生成网格
Mesh 2;
