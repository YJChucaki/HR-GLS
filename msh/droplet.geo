// Gmsh script to draw a teardrop shape with a larger front (bottom circle) and smaller tail (top point)

// Define points
Point(1) = {0, 0, 0};      // 半圆中心
Point(2) = {-1.5, 0, 0};   // 半圆左侧端点 (A)，半径增大到 1.5
Point(3) = {0, -1.5, 0};   // 半圆底部点 (B)，半径增大到 1.5
Point(4) = {1.5, 0, 0};    // 半圆右侧端点 (C)，半径增大到 1.5
Point(5) = {0, 2, 0};      // 泪滴顶部尖点 (P)，高度减少到 2，使尾部更小
Point(6) = {-0.75, 1, 0};  // 左侧样条中间控制点 (D)，调整以匹配新比例
Point(7) = {0.75, 1, 0};   // 右侧样条中间控制点 (E)，调整以匹配新比例

// Define arcs for the bottom semicircle (larger radius)
Circle(1) = {2, 1, 3};     // 从点 2 到点 3 的弧，中心在点 1
Circle(2) = {3, 1, 4};     // 从点 3 到点 4 的弧，中心在点 1

// Define splines for the sides (steeper for smaller tail)
Spline(3) = {4, 7, 5};     // 从点 4 经点 7 到点 5 的样条曲线（右侧）
Spline(4) = {5, 6, 2};     // 从点 5 经点 6 到点 2 的样条曲线（左侧）

// Define curve loop (顺时针方向)
Curve Loop(1) = {1, 2, 3, 4};

// Define plane surface
Plane Surface(1) = {1};

// Set color to blue (optional)
Color Blue { Surface{1}; }

//Mesh.Algorithm = 1;
//Mesh.MshFileVersion = 2;
//Mesh 2;
//Mesh.SecondOrderIncomplete = 1;
//SetOrder 1;
//RecombineMesh;
