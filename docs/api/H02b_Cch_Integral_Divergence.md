# H02b / CCH 唯一积分散度

状态：2026-10-03，**H02b / H02 PASS**；[任务书](../task/H02b_Integral_Divergence_TODO.md)、[验收报告](../task/H02b_Integral_Divergence_Report.md)。范围为冻结 CCH §6.1 的原生积分散度；其后 [H03](H03_Cch_Pressure_Load.md)已从此算子生成压力载荷并通过伴随验收，[H04](H04_Cch_Single_Phase_Saddle.md)已完成封闭域鞍点参考求解。[H05](H05_Cch_Physical_Divergence.md)已消费同一D_I验证真实物理散度与轴规则，[H06](H06_Cch_Energy_Viscosity.md)已 PASS，下一项 H09，尚无 CCH 流体 App 或完整时间推进。

## 入口与所有权

公开头为 `include/Modules/Hydrodynamics/CchIntegralDivergence.h`，类型位于 `Quantum::Hydrodynamics`，消费者链接 [lamp_hydrodynamics_cch](H02a_Cch_Build_Boundary.md)。

```cpp
const auto layout = CchDofLayout::build(mesh, domain);
const auto divergence = CchIntegralDivergence::build(mesh, layout);
const CchNativeValues values(layout, nativeNormalVelocity, cellPressure);
const auto netOutflow = divergence.apply(mesh, values);
const auto& integralMatrix = divergence.matrix();
```

构造先对 actual mesh 验证 retained layout，再逐个 participating native face 查询同一 Core 几何，装配 compressed sparse matrix。只保存稀疏系数和不可变 layout handle，不借用 mesh 指针、不保存第二份 face geometry、B 或散度 stencil。复制算子共享同一只读矩阵与 layout identity；原 mesh/domain 对象析构后仍可读取矩阵，应用时必须提供身份完全等价的实际 mesh。

## 数学与编号

\[
D_I=B A_f,\qquad (D_Iu)_c=\sum_{f\subset\partial c}s_{cf}A_fu_f.
\]

矩阵为 Np×Nu，行/列直接采用 [H01](H01_Cch_Native_Dof_Layout.md) 的 pressure/velocity 编号。面速度取正 x/r 或 y/z 方向：negative cell 写入 +A_f，positive cell 写入 -A_f。同一内部面只有一个速度列，其两项系数精确相反；selected 域边界只有所选 cell 的贡献，包括未选物理邻居形成的边界。

粗侧面由两个真实 fine subfaces 求和，无重叠 coarse-face unknown 或 virtual-child 插值。XY 的 A_f 为边长（单位厚度）；RZ 的 A_f 为 full 2πr 旋转测度。轴上径向面 A_f=0，保留 H01 velocity 列，不存显式零项；这不代替后续轴正则和边界消元。

返回的是各 cell 的积分净外流量：XY 为单位厚度体积流率，RZ 为三维体积流率。不除以 cell volume，不乘 dt、密度、相分数或压力；pressure 数据改变不影响结果。nnz≤2Nu 且≤8Np，矩阵乘法随 nnz 线性增长；mesh identity 在每次批次入口验证一次。

## 只读消费与拒绝

| 接口 | 契约 |
| --- | --- |
| build(mesh, layout) | 校验同一 actual mesh，准备唯一积分矩阵 |
| layout() | 返回 retained CchDofLayout 的 const 引用 |
| matrix() | 返回同一 `Eigen::SparseMatrix<double, RowMajor, Eigen::Index>` 的 const 引用；供 H03 直接负转置 |
| apply(mesh, values) | 验证 mesh 内容及 values 的同一 layout identity，使用上述矩阵乘 velocity，返回长度 Np 的有限结果 |

同尺寸但几何/连接不同的 stale mesh、独立 rebuild 的 layout、XY/RZ 或 selection 错来源均抛出 invalid_argument。CchNativeValues 原有尺寸/有限输入检查继续有效；有限输入造成输出溢出也拒绝，矩阵和值不变，后续合法调用可继续。

matrix() 是 prepared operator 的只读代数端口，直接做外部矩阵乘法不执行来源验证；携带状态的公开计算应使用 apply。它不是额外可写状态或另一套算子装配入口。

75 个场景和 2,400 组随机/零/常量速度通过独立边积分、解析多项式积分、流函数局部闭合、真实 AMR 子面求和和身份/寿命验证，所有归一化缺陷≤1e-13。算子重建测试只验证重新准备的矩阵，不证明物理速度 transfer、压力求解或 Gate H。[H03](H03_Cch_Pressure_Load.md)已实现 C=-D_I^T 并通过伴随验收，保留同一权威散度 handle；[H04](H04_Cch_Single_Phase_Saddle.md)单独验证封闭域参考求解及均匀 XY pressure/velocity MMS，[H06](H06_Cch_Energy_Viscosity.md)已 PASS，下一项 H09。
