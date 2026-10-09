# RZ 双液滴毛细运动

完整物理基准的冻结/实验耦合配置。新 FirstJet–Hydro/EHD 耦合尚未验收，本目录的存在不代表低黏/无黏或破碎能力已完成。

```powershell
.\tools\build.cmd
.\tools\run_case.cmd hydro_moment_two_drops_rz
```

输入为对应 system 字典。修改模板后重新运行，结果自动归 runs/<名称，斜杠转双下划线>/<时间戳_标识>；预检加 -CheckOnly。

[详细参数与历史验收说明](../../docs/history_todo/cases/hydro_moment_two_drops_rz/README.md)。其中的旧脚本仅作追溯。

统一操作见 [运行与调试](../../docs/guides/Run_And_Debug.md)；配置与能力状态见 [算例清单](../catalog.json) 。
