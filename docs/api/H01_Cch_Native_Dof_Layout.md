# H01 / CCH 原生自由度与数据 contract

状态：H01 PASS，2026-10-03；验收见 [报告](../task/H01_Native_Face_Cell_Dof_Layout_Report.md)。其后 [H02a](H02a_Cch_Build_Boundary.md)已建立真实编译入口，[H02b](H02b_Cch_Integral_Divergence.md)实现唯一积分散度，[H03](H03_Cch_Pressure_Load.md)实现同源负转置压力载荷并通过伴随验收，[H04](H04_Cch_Single_Phase_Saddle.md)实现封闭域单相鞍点参考求解。[H05](H05_Cch_Physical_Divergence.md)已消费同一D_I验证真实物理散度与轴规则，[H06](H06_Cch_Energy_Viscosity.md)已 PASS，下一项 H09，尚无可运行的 CCH 流体 App。

范围：只有原生编号、数据所有权和身份验证；不装配流体算子，不注册可运行的 CCH App。任务与验收见 [H01 TODO](../task/H01_Native_Face_Cell_Dof_Layout_TODO.md)。冻结 CCH 方法 §3.1 与 H-R1 规定面速度/单元压力；本 contract 直接复用 Core 的 NativeLeafDomain、NativeMeshIdentity 和真实 2:1 subfaces。

## 入口与 owner

公开头为 `include/Modules/Hydrodynamics/CchDofLayout.h`，类型位于 `Quantum::Hydrodynamics`。这是紧凑容器 contract；build 已原样移入对应 .cpp，消费者链接真实 `lamp_hydrodynamics_cch`，其余短容器操作保持头内实现。测试经该库依赖 Core，既有 Hydro/Interface Compatibility 库没有成为其依赖。

```cpp
const auto layout = CchDofLayout::build(mesh, domain);
const CchNativeValues values(layout, nativeNormalVelocity, cellPressure);
values.validateFor(layout, mesh);
```

`domain` 必须先由已验证的 actual mesh/FaceSegmentGraph 构造。layout 拥有其 compact selected-cell metadata 副本并保留同一个 Core mesh identity，不保留 mesh 指针、不复制面几何。正常查询编号和读取字段时，原 domain/mesh 的寿命不影响布局和数据；复用实际几何前必须对所用 mesh 调用 validateMesh / validateFor。

每个 `build` 产生独立不可变 layout identity；复制 layout 共享 identity。Core 验证精确几何/拓扑内容，因此内容完全等价的独立 mesh 对象可验证；同尺寸但坐标/连接不同的对象必须拒绝。即使两次 build 的尺寸和几何相同，也不能用新 layout 静默解释旧 values，须由后续显式 transfer 建立新值。

## 编号与测度

| 接口 | 含义 |
| --- | --- |
| velocityCount / pressureCount | 所选域参与的实际面数 / pressure cell 数 |
| primaryScalarCount | Nu + Np；不包含派生重构或任何高阶全局状态 |
| velocityDof(nativeFace) / nativeFace(velocityDof) | native face 与紧凑速度编号双向映射 |
| pressureDof(nativeCell) / nativeCell(pressureDof) | 直接复用 retained domain 的唯一 pressure 编号 |
| normalAxis(velocityDof) | First 是 x/r 分量，Second 是 y/z 分量，均取正坐标方向 |
| domain | 同一 selected NativeLeafDomain 的只读引用，测度继续由 Core 给出 |
| sameIdentity / validateMesh | 布局同源检查 / actual mesh 完整身份验证 |

native row 沿用现有 vertical edges 在前、horizontal edges 在后的编号；不新增另一套 physical face ordering。有效但完全不参与所选域的 native face/cell 返回 -1；负索引和超范围索引明确拒绝。内部面只占一个速度数据位置，两侧 cell 共享；coarse–fine 接口只有实际 fine subfaces，不能再增加重叠 coarse-face unknown。压力不放在顶点、face 或 hanging 点。

XY/RZ 在本项只影响坐标名称与 retained Core domain。RZ 轴面即使 area=0，仍保留原生速度数据位置；它不意味着后续可独立求解的轴上自由速度。边界固定/消元、pressure gauge、轴正则和物理边界分类在后续对应任务实现。

均匀 n×n 全域有 Np=n²、Nu=2n(n+1)。任意 selected 2:1 leaf 至多八个侧段，因此 Nu≤8Np；全域/部分域编号、metadata 和 primary values 均随 native mesh 线性增长。layout 保存索引/axis 标签与 Core domain，不保存 C1/Q8、face moments 或几何重构系数。

## 字段隔离与拒绝

`CchNativeValues` 只拥有长度 Nu 的 `normalVelocity` 和长度 Np 的 `pressure`。构造复制调用者输入，返回 const vector 引用；没有可写别名、gradient、库存或 ghost primary state。值可以精确为零、任意有限正负数；不引入物理幅值限制。

错尺寸、NaN/Inf、不同 layout identity 或 stale mesh 抛出 invalid_argument，不改调用者输入、不裁剪或替换数值。values 没有主时钟、candidate/accepted 或自动 commit；H14 之后由 App/Coupling 决定时间层与接受事务。

H02 已从这份 layout 和 Core face geometry 生成[唯一 integral divergence](H02b_Cch_Integral_Divergence.md)，[H03](H03_Cch_Pressure_Load.md)已从同一散度生成负转置压力载荷，[H04](H04_Cch_Single_Phase_Saddle.md)消费两者并完成封闭域 normal velocity 消元与单体参考求解。[H05](H05_Cch_Physical_Divergence.md)已通过，[H06](H06_Cch_Energy_Viscosity.md)已 PASS，下一项 H09；AMR Hodge、质量/动量通量、完整 production 边界和时间推进仍待对应任务，不能新装配第二套散度/gradient。
