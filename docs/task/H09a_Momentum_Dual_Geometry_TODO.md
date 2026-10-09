# H09a / 原生交错动量控制体与同源几何质量

状态：**PASS**，2026-10-03；验收见[报告](H09a_Momentum_Dual_Geometry_Report.md)。前置 H01–H08 PASS；依据[路线图](../roadmap/05_LAMP_Module_Task_Roadmap.md)、冻结 CCH §5.2 与[完全非守恒材料设计](H07_Pure_Gals_Hydro_Coupling_Design.md)。H09 拆为本项几何/质量和 H09b 共享流量/非零缺陷/动量一致性；两项均通过才能关闭 H09。

## 目标

只解决一个数学事实：每个 native face 速度拥有由相邻实际压力单元半体积组成的动量控制体；在 XY/RZ、真实 2:1 和选中子域上，使用唯一 FirstJet Q4 对该控制体及各真实交集分相积分，得到只读几何材料质量。

## 允许修改

- 新增 `include/Modules/Hydrodynamics/CchMomentumGeometry.h`、`src/Modules/Hydrodynamics/CchMomentumGeometry.cpp`，同一文件对内包含几何和材料质量两个紧密相关的只读 contract。
- root `CMakeLists.txt` 唯一 source 注册；`tests/contracts/hydrodynamics/cch_momentum_geometry_contract.cpp` 与 `tests/CMakeLists.txt`。
- 本 TODO、API、报告、现行路线图/文档索引/耦合设计；证据 `runs/baselines/H09a_20261003`，临时工作 `runs/checks/H09a/<时间戳>`。

预算：3 production/build +2 tests。超过预算先拆任务。

## 禁止修改

H01–H08 实现/旧测试、Interface 状态/输运/几何实现、Core/Solvers、旧 Hydro/Compatibility、电场/App、cases、冻结 PDF、历史 TODO/报告。不链接旧 momentum-dual owner。不新增全局速度自由度、材料库存、第二 scalar 或第二 primal 散度；不按中心判相、不用平均 cell density 代替 cut dual 积分、不修复体积/质量。

## 数学约束

- 对速度方向 a，各实际 cell 的两侧 native faces 分别拥有从实际 cell 中心到该侧的半矩形；2:1 的双子面再沿切向分开。每方向的这些交集精确覆盖该 cell，且不重叠。
- 同一 native face 的一或两个 owner 交集合成实际矩形 dual；选中域边界只有内侧半体积。共同边段仅保存一次，内段有两个同速度方向的 owner，外段一个；各 dual 四边端点精确覆盖，RZ 轴边面积精确为零。
- XY 用 dx*dy，RZ 用 pi*(r_hi+r_lo)*dr*dz；RZ 径向半体积一般不相等。几何转移权重为实际交集 V_fc/V_c，不能用统一 1/2。
- phase 取同一 retained FirstJet publication 的 Q4 裁剪区域，复用正式 `HermiteGalsPhaseRegionSource`；开放端点声明从 retained curve metadata 读取。禁止重拟合、复制 scalar 或初始化 M。
- `m_f=rho_l*V_l,f+rho_g*V_g,f`；液气体积分别保留，小相不得由两个大数相减消失。密度须有限正数，无 floor/真空。质量是当前层几何派生量，不是守恒更新时间状态。
- 每方向 dual 液气总积分等于同层 cell 液气总积分至测度误差；真实交集相积分的 cell 分区误差显式报告。该空间积分自洽检查不要求时间连续性缺陷为零。
- 固定 layout/D_I、publication identity 和完整 time layer 保留；AMR 后旧对象拒绝，须显式重建。本项不实现跨网格速度转移或 H15 Hodge。

## 实现要求

- [x] TODO 先于源码；核验 H08 交付与冻结源码/输入/PDF/程序。
- [x] 不可变几何/质量对象，真实交集、共享边段、来源校验、finite/underflow/overflow 拒绝。
- [x] 独立几何枚举/解析相积分参考，覆盖 XY/RZ、uniform/2:1、多层AMR、选中半域/孔洞/不连通域；纯相、平面、圆/球、细小液气相及密度比1/1e3/1e6。
- [x] 等效 mesh/原 mesh 析构寿命、身份/完整时间层/材料错误、实际 refine 后拒绝与重建；不给旧数据重新贴标签。
- [x] 标准构建、短契约和前置回归；来源/链接/数值/文档哈希审计归档。

## 验收

**T0：** local Release 无新增 warning/error；新 cpp 唯一属 CCH、无旧 Hydro/Compatibility/App 依赖；公开头自足；help/architecture/self_test/repository_layout/FirstJet App 通过。

**T1：** 独立从实际 cell/face 枚举各 dual owner box、方向、体积、V_fc/V_c 和边段覆盖；归一化几何/共享关联误差≤1e-12。独立解析纯相/平面/圆球裁剪积分与质量误差≤1e-10；逐 cell 分区及每方向总质量误差≤1e-10。包括在同一 cell 内两侧 dual 相含量显著不同，明确检出平均 density 替代。轴边零面积、小相非零、复制寿命、真实 AMR 失效/重建及误用拒绝通过。

**T2：** 非多项式圆/球 geometry 的 n=4/8/16 uniform/2:1 dual 材料质量独立解析误差递减且末档阶≥2。每档分别重建来源，不能插值旧 mass；测出 rebuilt AMR 几何误差。本项是短算子收敛检查，不是论文级完整流体 PDE benchmark。

## 失败规则

同源、覆盖、积分或预定收敛未通过即 FAIL，在允许范围修复，不放宽门槛。若需要修改 Interface 或超预算，先拆缺失 contract；H09 保持未关闭。只有确无可推进的外部依赖才 BLOCKED。

## 交付

真实 dual 几何/交集/共享边段、同源分相材料质量、短契约、API/TODO/报告及原始验证/哈希证据。H09b 负责唯一 primal→dual 材料流量、实际非零 dual 连续性缺陷、uR 一致性项与常速度保持；H10/H14/H15 和完整两相/低黏/EHD 仍未完成。总目标保持 active。
