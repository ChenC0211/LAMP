# H08a2r13 / 公共来源观察器迁移设计

实施前设计，2026-10-05。仅研究驱动一文件及其独立 Python 参考一文件；production/test/CMake/input 不改。任务入口见 [TODO](H08a2r13_Benchmark_Source_Observer_TODO.md)。

审计发现 driver 的 `spatialAdvectorAudit` 导出旧 patch rank/deficit/support；`spatialClosureAudit` 从 patch 权重重算旧 slope。当前 cases/tests/tools 没有这些 CSV 的其他活动读者。历史 H10j0/H10l1 PowerShell 参考依赖三系数、旧 rank 库存及 trace/lift slope，只保留在历史证据中，不再作为当前来源参考。独立 `spatial_source_patches.csv` 是 H10i 几何研究输出，不访问 production patches，保持原样。momentum 短契约的 patch 指针比较归 r11 部署迁移，不在本任务改动。

`advector_faces.csv` 删除末尾 rank/deficit/support；其他物理值、身份、完整时钟保持。准备摘要继续原样导出现有公开 report，其计数只是实际 owner 的准备工作量，不能用于跳过物理断言。原 n=4 hole/disconnected 八个不适用解析精度布局按原 fixture 固定列出，不再由 patch deficit 动态选择；原56适用布局保持验收。

`closure_faces.csv` 删除 mean_x/mean_y/scale/slope，保留原208单位面通量基输入、物理边界、native cell 邻接和面积。cells/status 原样保留。新增 closure_atoms（实际完整双端五系数/几何/积分）、closure_views（真实 primal/dual atom 映射及双侧索引）、closure_sources（原身份/full clock/真实规模）和 closure_samples（每 atom 三点物理速度查询）。primal atom 的采样取所属真实 native cell；dual atom 使用原 owner。只读公共 field，不复制拟合、重构或生成伪旧权重。

独立参考直接用 GL3 对完整四次多项式乘 XY/RZ 测度积分，双端图表互验并与物理查询比较，再按实际共享 atom 组合 native 面、内部反号/真实 AMR 子面及 cell 边界积分。零 gross cell 仍严格检查原 D_I 为零；浮点积分按原1e-12无常数地板的完整场工作量预算，不强制消去有符号舍入。参考保留原56适用布局的解析场/平均散度检查，通用来源积分覆盖全部64布局。边界模式以相率、真实共享 atom、native 总量验证；旧 cell 体积坐标割线只作为历史诊断，不能当作新场点散度。

T0：标准研究构建、3项架构/目录。T1：同一当前生产库运行 before/after 的 advector、closure、boundary、trace-contract；未改表逐字比较，仅准备/更新时间单列。所有新 schema 严格读回、行/ID/覆盖计数；三/五系数及双端图表的高项截断负控。T2：新独立参考验收、原 phase/material/dual/momentum/KEP/characteristic 契约；另用冻结 r11 公共来源候选隔离编译当前研究驱动，验证没有残留 patch 构建依赖及原208basis物理参考。候选仍隔离，不作为部署或 PDE 通过。

临时工作全部位于 `runs/checks/H08a2r13/20261005_184026_443`，冻结前保存原 driver/runtime/输入哈希。完整 full benchmark、原 affine/native 及28窗口84gate在r11实际部署后执行，本任务不改其代码与门槛。

验收复核补充（实施五点导出前）：GL3的三个值节点不能唯一确定四次法向迹；在原GL3外增加两个不同内点，`closure_samples` 每atom五行，以五点唯一约束完整四次物理查询。积分仍为GL3。增加“在GL3节点为零的三/四次扰动”负控，避免积分和三点查询同时漏检高项。原三点过渡运行保留为中间证据，交付绑定final五点运行。
