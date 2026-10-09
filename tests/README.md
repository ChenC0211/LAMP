# 短检查

`contracts/<模块>` 放局部算子、状态/所有权、字典、输出等短检查；`support` 是测试辅助代码；`fixtures` 是最小静态输入；`architecture` 和 `repository_layout.cmake` 检查依赖与目录规则。每项由 CMake 显式注册。

```powershell
.\tools\build.cmd -Test
```

定向检查可用标准 `ctest --preset local -R interface_first_jet --output-on-failure`。CMake/CTest 不在 PATH 时，使用 Visual Studio 附带的二进制；通用构建入口会自动寻找。

完整 FirstJet 窗口与收敛研究已移入 [cases/benchmarks](../cases/benchmarks/README.md)。旧失败路线的专属 harness 和无配置的 App 对比入口已退役，快照及清单在 `runs/baselines/G02_20261002`。

两个依赖原始二进制 fixture 的历史检查 `legacy_nodal_state_contract`、`pic_checkpoint_contract` 仍显示为 `missing_fixture` / Disabled；不计作通过。恢复所需原始文件后才启用。默认检查输出在 `runs/checks/<检查>`，临时文件不写进项目根目录或 cases。
