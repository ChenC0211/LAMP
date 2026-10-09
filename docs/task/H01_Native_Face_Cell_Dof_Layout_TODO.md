# H01 / 原生面速度与单元压力自由度布局

任务依据：[现行分模块路线图](../roadmap/05_LAMP_Module_Task_Roadmap.md) H01。执行日期：2026-10-03。前置 Gate I 已通过；本项不改变 Interface 数学实现。完成后按顺序进入 H02，不在本任务实现流体推进。

## 目标

只解决一个架构事实：CCH 主数据为真实 native face 上的单个正坐标方向法向速度与 selected pressure cell 上的单个压力；自由度随 native mesh 线性增长，布局与输入数据有明确的不可变身份和寿命，不依赖全局 C1/Q8 速度状态。

## 允许修改

- `include/Modules/Hydrodynamics/CchDofLayout.h`：紧凑不可变编号和同布局的只读 native values；容器型 contract 直接实现，不为此增加算法翻译单元。
- `tests/contracts/hydrodynamics/cch_dof_layout_contract.cpp`、`tests/CMakeLists.txt`：独立编号、规模、AMR、选择域和身份/寿命拒绝契约；测试只链接 Core，不链接现有 Hydro/Interface Compatibility。
- 本 TODO、对应报告、`docs/api/H01_Cch_Native_Dof_Layout.md`、路线图和必要文档索引。
- 构建与验收证据存入 `runs/baselines/H01_20261003`，临时读取/渲染存入 `runs/checks/H01/<时间戳>`。

预算：1 个 production 文件、2 个测试文件。实现按已有 NativeLeafDomain 和 NativeMeshIdentity 复用实际几何，不复制网格、face geometry 或另设几何 identity。

## 禁止修改

- Core/Solvers、冻结的 Hydro CompatibleMoment/Fitted/legacy 实现、Interface/Electricity/Plasma、App、字典 reader、case 物理参数和现有接受容限。
- 不装配 H02 散度、H03 压力载荷、H04 saddle solve、H05 RZ 测度/轴正则处理、H06 黏性、H07–H15 耦合/输运/受力/推进。
- 不增加高阶全局速度、coarse-face 额外未知量、ghost primary data、旧路线 fallback 或可运行的占位 cch 方法选项。

## 数学约束

- 每个参与所选域的实际 native face 恰有一个速度 DOF；2:1 coarse–fine 处只有实际 fine subfaces，没有重叠 coarse-face DOF。内部面由两侧共享。
- 每个 selected native leaf cell 恰有一个压力 DOF；不在顶点、面或 hanging 点增设压力。布局不做 pressure gauge 或边界消元。
- 速度方向固定为第一/第二坐标的正方向，即 XY 的 x/y、RZ 的 r/z；outward sign 仍属于 Core face geometry，不能变成第二套面定向。
- 全域均匀 n×n 网格满足 Np=n²、Nu=2n(n+1)。任意 selected 2:1 域满足 Nu≤8Np；所有存储随 native mesh 线性增长。
- RZ 轴面保留一个原生速度数据位置，即使真实面测度为零；本项不生成伪测度，也不在未接入边界条件前宣称它是独立可求解自由度。
- layout 保留不可变 selected-domain 数据与 Core mesh identity；相同维数不等于相同 mesh/selection/measure/layout。已保存字段不能用于不同布局或过期网格。

## 实现要求

- [x] 核对冻结 CCH 方法中原生未知量定义、H-R1 和 Core native face 编号，保存开工前文件/程序身份。
- [x] layout 构造只使用已验证 NativeLeafDomain；velocity 双向编号紧凑，pressure 编号复用该 domain 的唯一映射。
- [x] 只保留索引、轴类型和已存在的 domain；不复制 face geometry/mesh 或引入物理边界/材料数据。
- [x] native values 拥有 velocity/pressure 的副本，向外只读；尺寸、非有限数据、错误 layout 和 stale mesh 均拒绝，不修正输入。
- [x] 验证全域/部分域、XY/RZ、均匀/2:1 AMR、真实 coarse–fine 面共享、输入及 mesh 销毁后的寿命、相同尺寸不同来源的拒绝。
- [x] 完成标准 CMake 构建、T0/T1/T2，记录文件预算、未实现边界与验收结论。

## 验收

- **T0：** 标准 preset 配置；新公开头可独立作为首个 include；契约构建无新增 warning；architecture、自测、repository_layout、lamp_help 通过。
- **T1：** 三档以上均匀网格的解析计数；AMR/selected 域与独立 FaceSegmentGraph 的唯一面集合完全一致；轴分量/pressure 编号双射；primary 标量数恰为 Nu+Np，无高阶数据；输入复制和只读约束；非法 native/dof 索引、空 domain、错维数、NaN/Inf、不同 layout、改动后的同维数 mesh 及 AMR 改变均拒绝。
- **T2：** 现有短 Interface App/输出契约及代表性 Hydro native 几何契约通过；旧 production 源码和原输入哈希不变。H01 没有新的流体方程，因此不宣称流体 MMS、低黏稳定性或耦合 gate。

## 失败规则

若验收暴露本项编号、寿命或输入隔离错误，记 FAIL 并在范围内修复。若必须更改既有 Core 几何或旧物理算法，或超过预算才能继续，记 BLOCKED 并给出接口缺口；先拆分任务，不加 fallback。最终结论只用 PASS / FAIL / BLOCKED；未实现的 H02–H15 不能计作本项缺失或通过。

## 交付

公开 contract、独立短测试、API、TODO/报告、构建/检查日志和来源哈希。范围明确为布局与数据所有权；下一单任务为 H02。

**最终状态：PASS。** 23 个网格场景、29 个拒绝，7/7 定向 CTest；1 production + 2 测试文件，旧源码/程序/方法和 120 份原字典哈希保持。详见 [H01 报告](H01_Native_Face_Cell_Dof_Layout_Report.md)。
