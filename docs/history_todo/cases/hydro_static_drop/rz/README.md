此页为 G02 整理时保留的历史算例说明。数值、限制和既有结果供追溯；其中旧脚本、旧程序名和输出位置不再是当前操作入口。当前说明见 cases/hydro_static_drop/rz/README.md。

运行入口：tools\run_case.cmd hydro_static_drop/rz 。结果：runs/hydro_static_drop__rz/<时间戳_标识>/。程序：LAMP。

# B2-RZ：轴对称静止球滴

无旋轴对称两相不可压缩 NS；球心 (0,0)、R=0.4，域 r∈[0,1]、z∈[−1,1]。r=0 为规则轴，其余三边自由滑移。两相 ρ=1、μ=0.0816496580927726，σ=1，La=120；初速零、压力由方程建立。采用真实旋转测度，总曲率包括子午和环向两部分，平衡压差为 5。

32×64 均匀网格、32 个子午弧段，dt=0.004898979485566356，T=7.83836717690617，共 1600 步，每 40 步写场。这是项目轴对称扩展；参考、模型及判据见[统一说明](../../../../../cases/hydro_static_drop/README.md)。

输入：[子午线公式](../../../../../cases/hydro_static_drop/rz/system/geometry.dict)、[物性与边界](../../../../../cases/hydro_static_drop/rz/system/hydro.dict)、[界面](../../../../../cases/hydro_static_drop/rz/system/interface.dict)、[网格](../../../../../cases/hydro_static_drop/rz/system/mesh.dict)、[运行](../../../../../cases/hydro_static_drop/rz/system/run.dict)。

```powershell
.\tools\run_case.cmd hydro_static_drop/rz
```

现行初值为隐式 r²+z²−R²=0，从下轴极点 (0,−0.4) 追踪，traceTolerance=1e−10；提取后仍使用同一 Hermite 子午线。结果在本目录 output/<时间戳>/。当前只有[隐式输入预检](output/20260928_233052_458_implicit_input_fe6204/run.log)，没有本版物理运行结果；[原参数输入预检](output/20260928_230936_057_analytic_input_c22b7d/run.log)保留追溯。
