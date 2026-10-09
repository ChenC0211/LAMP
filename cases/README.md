# LAMP 完整算例与 benchmark

这里放完整验证、应用算例和研究配置。小型契约和最小 fixture 放 [tests](../tests/README.md)，结果放 `runs/<名称>/<时间戳_标识>`。除 `benchmarks` 的专属研究驱动外，本目录只放静态配置、输入、参考和说明。

```powershell
.\tools\build.cmd
.\tools\run_case.cmd -List
.\tools\run_case.cmd hermite_gals_vortex_xy_dynamic_amr
```

| 家族 | 目标与状态 |
| --- | --- |
| `hermite_gals_*` | 当前纯 FirstJet 通用配置；六类模板 × 四种固定/动态网格，历史 Gate I 已通过；解析验收由独立 benchmark 执行 |
| [hermite_gals_vortex_long_time](hermite_gals_vortex_long_time/README.md) | 2,000 步可逆差速旋涡，动态/固定 AMR 对照；独立解析曲率、位置和面积验证 |
| [benchmarks](benchmarks/README.md) | 完整时间窗、空间/时间收敛与 AMR 研究；显式构建/运行 |
| `hydro_static_drop`、`hydro_drop_oscillation`、`hydro_hysing_bubble` | 完整 Hydro 基准；参考数据随对应研究保存 |
| `hydro_moment_*`、`hydro_column_*`、`hydro_conservative_rz` | 冻结/实验耦合配置；各 README 保留已知能力及限制 |
| `electrospray` | 电喷雾工作流配置；新的 FirstJet–Hydro/EHD 完整耦合尚未验收 |

全部名称、App 和能力状态见 [catalog.json](catalog.json)。`run_case.cmd -List` 从实际配置和 benchmark 描述生成清单。无配置的旧 smoke 注册及已删除路线的配置不再出现在入口中。

编辑模板参数后重新运行；不在模板下执行物理求解。旧结果位置见 [G01 迁移索引](../runs/MIGRATION_G01.md)。
