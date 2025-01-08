R = 0.3;  // 定义圆弧半径为0.3
L = 1.0;  // 定义正方形边长的一半为1.0

Point(1) = {L, -L, 0};  // 创建点1，坐标为(1.0, -1.0, 0)
Point(2) = {L, L, 0};   // 创建点2，坐标为(1.0, 1.0, 0)
Point(3) = {-L, L, 0};  // 创建点3，坐标为(-1.0, 1.0, 0)
Point(4) = {-L, -L, 0}; // 创建点4，坐标为(-1.0, -1.0, 0)
Point(5) = {-L + R, -L, 0}; // 创建点5，坐标为(-1.0 + 0.3, -1.0, 0)
Point(6) = {-L, -L + R, 0}; // 创建点6，坐标为(-1.0, -1.0 + 0.3, 0)
Point(7) = {-L + Cos(Pi/4) * R, -L + Sin(Pi/4) * R, 0}; // 创建点7，坐标为(-1.0 + 0.3*sqrt(2)/2, -1.0 + 0.3*sqrt(2)/2, 0)

Circle(1) = {5, 4, 7};  // 创建圆弧1，连接点5、4和7
Circle(2) = {7, 4, 6};  // 创建圆弧2，连接点7、4和6

Line(3) = {6, 3};  // 创建线段3，连接点6和3
Line(4) = {3, 2};  // 创建线段4，连接点3和2
Line(5) = {2, 1};  // 创建线段5，连接点2和1
Line(6) = {1, 5};  // 创建线段6，连接点1和5
Line(7) = {2, 7};  // 创建线段7，连接点2和7

Curve Loop(1) = {4, 7, 2, 3}; // 创建曲线环1，包含线段4、7、2和3
Plane Surface(1) = {1};       // 创建平面表面1，使用曲线环1

Curve Loop(2) = {7, -1, -6, -5}; // 创建曲线环2，包含线段7和负数索引的线段
Plane Surface(2) = {2};          // 创建平面表面2，使用曲线环2

Transfinite Line{3, 4, 5, 6} = 15; // 设置线段1到7的跨越点数为10
Transfinite Line{7}=15;
Transfinite Line{1,2}=5;

Transfinite Surface{1}; // 对表面1进行超有限元网格划分
Transfinite Surface{2}; // 对表面2进行超有限元网格划分

Mesh.ElementOrder = 1; // 设置网格元素的阶数为1
Mesh.Algorithm = 8;    // 设置网格生成算法为8

// EOF  // 文件结束标记//+
Physical Curve("left") = {3};
//+
Physical Curve("right") = {5};
//+
Physical Curve("up") = {4};
//+
Physical Curve("bottom") = {6};
//+
Physical Curve("circle") = {2, 1};
//+
Physical Surface("surface") = {2, 1};
