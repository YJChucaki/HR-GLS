// Gmsh脚本：绘制半个阴阳鱼（太极图的上半部分）

// 设置基本参数
r = 1.0;  // 大圆的半径
h = 0.01; // 网格精度

// 定义大圆的圆心和边界点
Point(1) = {0, 0, 0, h};    // 圆心
Point(2) = {r, 0, 0, h};    // 大圆右端
Point(3) = {0, r, 0, h};    // 大圆上端
Point(4) = {-r, 0, 0, h};   // 大圆左端

// 定义小圆的圆心和边界点（小圆半径为大圆的一半）
Point(5) = {0, r/2, 0, h};  // 小圆圆心（上半部分）
Point(6) = {0, r, 0, h};    // 小圆上端
Point(7) = {r/2, r/2, 0, h}; // 小圆右端
Point(8) = {-r/2, r/2, 0, h}; // 小圆左端

// 定义大圆的圆弧（上半部分）
Circle(1) = {2, 1, 3};  // 从右到上
Circle(2) = {3, 1, 4};  // 从上到左（完整大圆的上半部分用这两个弧）

// 定义小圆的圆弧（上半部分的小圆）
Circle(3) = {6, 5, 7};  // 小圆的上半部分弧
Circle(4) = {7, 5, 8};  // 小圆的右到左弧
Circle(5) = {8, 5, 6};  // 小圆的下半部分弧（完整小圆）

// 定义曲线循环并创建面（上半部分）
Line Loop(1) = {1, -3, -4, -5}; // 大圆上半弧减去小圆
Plane Surface(1) = {1};         // 创建面

// Set color to blue (optional)
Color Blue { Surface{1}; }

//Mesh.Algorithm = 1;
//Mesh.MshFileVersion = 2;
//Mesh 2;
//Mesh.SecondOrderIncomplete = 1;
//SetOrder 1;
//RecombineMesh;
