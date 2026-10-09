# H09a / 原生动量控制体与同源几何质量验收

日期：2026-10-03。**H09a PASS；下一项 H09b，H09 未关闭。** 按[任务书](H09a_Momentum_Dual_Geometry_TODO.md)实现当前层真实 dual 几何/材料质量；总目标、Gate H/IH/E/C 保持未完成。

## 实现与数学边界

[生产接口](../api/H09a_Cch_Momentum_Geometry.md)在正式 CCH 库提供 `CchMomentumGeometry/CchMomentumMass`。前者从实际压力 cell 中心和 native faces 建立一或两个半矩形交集，共同边段只存一次；2:1 双子面沿切向分开，selected 边界只取内侧。每速度方向精确分区同一 cell。RZ 用真实物理测度及 `V_fc/V_c`，不把径向半体积统一设为1/2。layout、原 D_I 及来源均保留，无新全局自由度或第二 primal 散度。

质量直接裁剪 retained FirstJet publication 的原 Q4，复用正式 `HermiteGalsPhaseRegionSource`；开放 contacts/symmetry 从已有 curve metadata 读取。液气体积和每个实际交集的相积分独立保存，`m_f=rho_l*V_l,f+rho_g*V_g,f` 为有限正 double。没有 M、守恒库存、中心判相、几何/材料修复或平均 cell density 替代；旧 Hydro owner 未被依赖或改写。

现行[耦合设计](H07_Pure_Gals_Hydro_Coupling_Design.md)进一步明确：真实 cut dual 质量变化通常不等于 `R_m*delta_m_cell`。后续须核算 dual 自身非零缺陷及该几何转移差异，不能只插值 cell 缺陷。H09a 不生成 dual 材料流量或 uR 项，不推进动量。

## 独立数值验收

| 检查 | 实测 |
| --- | --- |
| 几何 / 材料场景 | 26 / 117；XY/RZ、uniform/2:1 depth1/2、全域/半域/孔洞/不连通域 |
| 真实交集 / 独立共享矩形边 / segment检查 | 3,216 / 3,412 / 6,304 |
| 几何测度/权重最大归一化误差 | 3.2167e−15，门槛1e−12 |
| 解析分相体积 / 材料质量最大误差 | 1.4701e−15 / 1.3629e−15，门槛1e−10 |
| 逐cell/方向分区最大相对误差 | 8.8349e−12，门槛1e−10；包含非多项式相积分 |
| 开放 contact 三角形独立参考误差 | 8.8818e−16；非对称开放端点声明保持 |
| 与平均cell密度替代的最大差异 | 0.43999912，按 `(rho_l+rho_g)*V_dual` 归一化；真实积分检出该替代 |
| tiny phase / coarse–fine cut dual | 4组2^−40液气小相、36个非零质量贡献；XY/RZ分别68/55次实际粗细交界cut dual检查 |
| 明确拒绝 / 实际AMR失效与重建 | 37 / 2；非法密度/身份/完整clock/layout/D_I、overflow/underflow、坐标/拓扑变更 |

几何参考独立从 native face 的相邻实际 cell 枚举 owner box；共享边参考由同方向两矩形真实交集得到，核验精确端点覆盖及 reciprocal incidence。相体积参考使用独立圆/球裁剪原函数或平面/三角形积分，覆盖密度比1/1e3/1e6和纯相；未调用被测 region 积分作参考。复制共享只读数据，原 mesh 析构后的等效 mesh 可消费；真实 refine 后旧对象拒绝，新来源重建后准确。

非多项式 `exp(q)-1` 圆/球的独立解析 dual mass 误差，n=4/8/16；AMR 球实际穿过粗细交界，未让 cut 全部留在未细化区：

| 序列 | n=16归一化质量误差 | 末档阶 |
| --- | ---: | ---: |
| XY uniform | 2.0671e−9 | 4.6998 |
| XY 2:1 | 3.6023e−10 | 4.8530 |
| RZ uniform | 9.2242e−10 | 4.9990 |
| RZ 2:1 | 2.7417e−10 | 5.2025 |

四组均误差递减并超过预定阶≥2。每档从新 publication/mesh 重建材料，无旧质量插值。以上是算子空间验证，不是流体 PDE 收敛、完整物理benchmark或时间阶证明。

## 工程验证与证据

预算3 production/build +2 tests：新 header/cpp/root CMake、新短契约/tests CMake。前置1,031文件中1,029不变，仅两个CMake更新、新增3文件；312 case资产和31 fixture资产、H01–H08/Interface/Core/旧Hydro/Compatibility/App/PDF保持。LAMP.exe哈希与前置相同，新source唯一归CCH，依赖闭包仍是正式Core/Interface/Eigen/线代，无旧Hydro或Compatibility。

- local Release定向构建通过，无新增warning/error。
- 新契约 **1/1 PASS，18.14秒**。曾有非法clock fixture、1×1网格不满足已有mixed-closure五点支持的fixture失败，已修正输入并保存原日志；没有放宽门槛。通过后又加强实际AMR cut覆盖并重跑新契约。
- 前置/受影响回归 **16/16 PASS，206.13秒**：H01–H08、正式phase-region、FirstJet App、旧pressure-gradient load、help/architecture/self_test/repository_layout。H01–H08所有选定数值行与H08冻结证据逐字一致，App仍12smokes/49reject。
- 文档收尾后repository_layout与本地链接、来源/交付哈希另行核验。未跑全部80项CTest、Debug或完整流体论文benchmark。

原始证据：[verification](../../runs/baselines/H09a_20261003/verification.json)、[operator detail](../../runs/baselines/H09a_20261003/operator_detail.log)、[regression detail](../../runs/baselines/H09a_20261003/regression_detail.log)、[source audit](../../runs/baselines/H09a_20261003/source_audit.json)、[交付清单](../../runs/baselines/H09a_20261003/delivery_manifest.json)。TODO-before、失败fixture和AMR覆盖加强前的日志保留；临时审计脚本归档到baseline/engineering_scripts。

## 未解决范围

H09b的唯一primal→dual材料流量、真实非零dual缺陷、uR项及高密度比常速度保持仍待实现。H10一般KEP对流/动能材料项、H14完整corrector、H15兼容AMR Hodge/速度转移/生产边界、CCH App及低黏两相/EHD均未通过。区域积分逐查询验证mesh的规模成本需要后续完整算例测量；本项未新增unchecked查询或改变Interface。
