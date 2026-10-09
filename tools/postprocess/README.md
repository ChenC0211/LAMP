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
