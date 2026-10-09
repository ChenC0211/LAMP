# I12a11r14 / 精确 Q4 细化来源 API

入口：`LevelSetHermiteGeometry::Snapshot::refineNativeQuartic(targetMesh, source)`。这是 Core 标量算子；Interface 的守恒库存转移与受检 regrid publication 尚未接入。

输入须为完整 active native Q4 Snapshot 和同一物理域上的纯细化网格。每个目标叶必须完全落在一个旧父叶内；允许不变叶、uniform→uniform、uniform→AMR 和 AMR→更细 AMR，允许细化界面单元。拒绝 coarsening、跨父叶重拟合、partial-active、corner-sector、非 Q4、域变化及 moved-from 来源。网格仍须满足现有 native 2:1 约束。

新 master 的 value/first/mixed jet、edge value/normal 和 center 保留旧来源的精确齐次比值；已有 cardinal datum 直接复用其精确表示。所有 hanging/fine-side 数据仍由原共享约束引擎产生。每个子叶通过原 Q4 限制内核共用一个父叶分母，避免把逐点清分母的重复因子相乘。没有浮点拟合、删小项或新的基函数。

返回前，每个叶的全部 25 个 normalized Hermite 数据都必须与旧多项式的限制精确相等，共享约束的 physical datum 也逐项核验。比较采用保留全部尾项的交叉相乘；必要的二进制幂缩放仅调整齐次表示。Q4 的唯一性保证整个叶多项式相同。几何查询的浮点容限不代替该证书。

结果拥有自己的数据和查询结构，可在旧 Snapshot、旧 mesh、目标 mesh 销毁后独立查询。构造期间的只读借用在返回前清除；异常不修改旧来源。继续遵守现有 exact expansion 的有限指数窗口和浮点环境检查，不承诺任意尺度、任意深度都可表示；不能精确表示时抛出异常，不返回近似替代结果。

该接口不生成 `HermiteGalsPrimaryState`、不更新 `M_K`、不执行 I10、不构造 `HermiteGalsPublication`，也不是 accepted-step commit。后续推进必须补齐独立守恒库存和 I11 的受检发布，不可通过低层 `bindGeometry/finish` 将这个 Core 结果冒充通过验收的 Interface 状态。

任务及执行证据见 [TODO](../task/I12a11r14_Exact_Quartic_Refinement_Source_TODO.md) 和 [报告](../task/I12a11r14_Exact_Quartic_Refinement_Source_Report.md)。父任务 I12 仍为 FAIL。
