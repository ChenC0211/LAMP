# 已保存结果的后处理

`analyze_hydro_standard.ps1` 调用 `HydroStandardAnalysis.cs` 读取已有 Hydro 输出与参考，支持 staticXY、staticRZ、oscillation、hysingTC1、hysingTC2。它不推进时间或修改已接受物理解；输入必须是 runs 下的终态运行目录。

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File tools/postprocess/analyze_hydro_standard.ps1 -Benchmark staticXY -Run runs/hydro_static_drop__xy/<时间戳_标识>
```

Hysing 参考在 cases/hydro_hysing_bubble/reference。生产 FirstJet 只报告几何体积漂移和同源曲率，没有专用解析误差；独立 benchmark 或下述参考验证器给出解析精度。ParaView 打开 solution.pvd。完整研究驱动的 CSV 记录在 run.log。

`analyze_interface_reference.ps1` 调用同目录 Python 验证器，只读终态运行和其冻结 reference/verification.json；当前支持 radialTwistCircle。NumPy/Matplotlib 由 ParaView pvpython 或显式 -Python 提供。用 -SelfTest 交叉检查逆映射、独立参数曲线和 finite differences；-AcceptedPrefix 仅允许失败运行的接受前缀诊断，不能通过完整窗口 gate。

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File tools/postprocess/analyze_interface_reference.ps1 -Run runs/<算例>/<时间戳>
```

输出在原运行 output/analysis/reference_<UTC时间戳>：来源 SHA-256、逐帧误差、验收 JSON、报告和 PNG/SVG 图。0=PASS、1=gate FAIL、2=来源/输入错误。所有 VTP 指标注明采样/折线近似，不回写 accepted geometry，也不构造另一套数值解。模板与方法见 [长旋涡研究](../../cases/hermite_gals_vortex_long_time/README.md)。

## Tecplot 新输出与历史离线导出

当前 VOF 新 run 的体场、界面与 series 只写 DAT；`output/meta/frames.csv` 是唯一帧索引，volume 直接指向 `.dat`。不生成新 VTU/PVD、`meta/tecplot_frames.csv` 或 series CSV；`meta/configuration_trace.csv` 等来源记录保留。体场为原节点和连接的 FEQUADRILATERAL/BLOCK，16 项实际单元中心字段使用 CELLCENTERED，不插值或重算。manifest 明示 `volume_format=tecplotDat`、`series_format=tecplotDat`、`frame_index=meta/frames.csv`。

原生 series DAT 是宽 POINT 表：先按原顺序保留数值列名，再按相同顺序追加每列 `<列名>Available`。原 NA、空值或非有限值写为 `0` 且对应 Available=0，这个 0 只是占位；真实 0/false 的 Available=1。画 XY 曲线时须同时用所选 X、Y 各自的 Available 列做 blanking/筛选。表中保留原记录，缺测不会自动断成独立曲线，不能直接把占位 0 当作测量值。

stage/status/reason 等分类保存在 zone AUXDATA 中，含可读值和原字段字节的可逆十六进制；分类变化或真实时钟非递增时另起 zone。性能按固定阶段写入 `performance_<stage>.dat`，例如 `performance_VOF_AC.dat`，没有重复的 `performance.dat`。没有实际数据时不创建文件，包括无重试/失败的 events 和 σ=0 的 vof_csf。

旧文件不会因新格式被删除。此前 CSV→DAT 的伴生曲线是三列（真实 X、Value、Step/Record）和每指标 zone；新原生宽表保留完整数值列，二者不能按同一 VARIABLES 假设读取。旧分析器须遵循其支持的原 schema，不能将旧 VTU/CSV 专用读取器当成新宽表解析器。

`export_tecplot.py` 只处理已保存运行，不推进数值、不插值、不修改既有物理输出；仅用 Python 标准库，工作文件放 D:。对 `frames.csv` 声明的原生体场 DAT，工具只验证格式、字段/位置、真实维度、连接范围及帧步数/时钟/几何身份后记为 native-existing/skip；不从 DAT 反造 VTU、不补第二份索引、不覆盖 DAT-only 运行。原生 series 保持原文件，工具不将其当 CSV 再转换。一般文件格式验证不等于几何/物理精度验收。

```powershell
python -B tools/postprocess/export_tecplot.py D:\Chen\LAMP\runs\hydro_rayleigh_taylor\<运行目录>
python -B tools/postprocess/export_tecplot.py D:\Chen\LAMP\runs\hydro_static_drop --retrofit
python -B tools/postprocess/export_tecplot.py D:\Chen\LAMP\runs\hydro_inviscid_drop_oscillation --retrofit --series-only
```

`--series-only` 只为历史 `output/series/*.csv` 补曲线，不扫描 VTU、体场或帧索引；配合 `--retrofit` 发现终态 Hydro/Interface 历史运行。运行身份仍须与原始/effective Run.app、launcher run.json、manifest.app 一致；排除 benchmarks、Electric/EHD/PIC 和含 electric.dict 的运行，身份不足/冲突明确 SKIP。只有原生 DAT 而无 CSV 时没有待补历史 CSV，保持文件不变。

`--csv-encoding` 只为编码已明确的旧 CSV 显式指定文本编码，例如 `--series-only --csv-encoding gb18030`。默认严格 `utf-8-sig`，不猜编码或替换错误字符；未知/非文本 codec 拒绝。DAT 为 UTF-8，非默认来源编码记入 CsvEncoding AUXDATA；原 CSV 不修改。

直接传单次运行目录；`--retrofit` 才递归发现。默认保留已有同名 DAT，`--overwrite` 仅能替换 Generator 和源路径均匹配的本工具/生产 CSV 伴生文件，拒绝陌生来源。成功退出 0，输入/格式/写入失败退出 2；此前完整写出的派生文件保留，可重跑补齐。

- 历史 VTU 路径：有 `output/meta/frames.csv` 时只读取索引 volume 列已发布的 VTU，没有索引才扫描 `output/fields/**/*.vtu`。支持本项目 ASCII UnstructuredGrid 单 Piece 的 VTK_LINE/VTK_QUAD，生成 FELINESEG/FEQUADRILATERAL、BLOCK、DOUBLE；保留节点、连接和原文本数值，CellData 仍是 CELLCENTERED。三分量拆 `_X/_Y/_Z`，使用 XY/RZ 两坐标（第三坐标须为零）。非有限字段仅在派生 DAT 中用 0+Available=0 明示，原 FieldData 时间/步数/身份等保存在 AUXDATA。
- 历史 CSV 路径：每文件生成同名 DAT，仅读取完整文本快照；末条记录无完整换行时 SKIP，其他坏行严格报错。三个变量为真实横轴字段名、Value、Step（有实际步数）或 Record（从 1 起的源逻辑记录号）；重名时第三列用 SourceStep/SourceRecord。X 优先 time、accepted_time、time_s、time_code、t_si_s、time_si_s、physical_time、t，再取实际 step/accepted_step/timestep/attempted_step，无 clock 才用 Record 并明确 AUX 来源。每指标/分类组合按有限连续段写 POINT zone，NA/空/非有限断段，X 非增另起段；全无有效数值不造空 DAT。true/false→1/0，分类保存在 AUX，原精度文本不窄化。

历史曲线已补齐：55 个已结束的非电 Hydro run、364 份 CSV 中，326 份有有效曲线、38 份仅表头；新增 294 个 DAT，既有 32 个 DAT 和原 CSV 哈希均保持。GB18030 历史 events 通过显式编码读取，19 个补转小样通过本机 preplot。这些结果属于历史格式；当前原生宽表的独立短合同和 preplot 结果见修复报告。

2026-10-18 的 CY/CX 源码已发布并集中构建，五项合同通过。两步 RT 与故意 Courant 失败前缀的 31 个实际 DAT（4 体场、4 界面、23 series）由本机 preplot 全部验读通过，无 warning/error；二进制 PLT v112 头中的变量、位置、zone 点数/单元数及每段真实行数与 ASCII 一致。另对原生 run 实际执行 `--overwrite` 得到 `written=0, native_existing=3`，33 个原输出文件哈希全部保持。证据为 `D:/Chen/Codex/temp/lamp_vof_cases_20261017/takeover_20261018/preplot/results.json` 和同目录 `native_export_verify.json`；这是输出协议检查，气泡完整窗口另判。
