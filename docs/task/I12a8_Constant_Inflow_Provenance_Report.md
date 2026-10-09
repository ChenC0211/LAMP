# I12a8 / 常相入流历史与区间身份报告

**PASS。** 常值单轴平移的实际入流面现在消费显式、不可变的纯相边界历史；半步、整步及其组合保留同一个来源。完整 first-jet driver 尚未接线，I12 仍 FAIL。

## 实现与范围

符合 [TODO](I12a8_Constant_Inflow_Provenance_TODO.md) 的 4 production + 2 test 文件预算。BC source 持有未绑定 BC 的 base interval，验证完整有序入流面、时间范围及旧 state/velocity/characteristic owner。interval、appended、整面库存装配均检查 BC 身份。独立创建的相同数值历史不能混用。

translation-face 入口对入流消费显式常相，对内部/出流保留 I12a7 扫掠；二者共享 strict Core total-volume 度量。没有更改 first jet、B/J、库存容限或边界液量。详见 [API](../api/I12a8_Constant_Inflow_Provenance_API.md)。

## 验证

- 主程序与 26 个测试目标 Release 构建通过，278.01 秒，无 warning；25/25 CTest 通过，160.12 秒。
- 新合同 **74 组检查、22 类拒绝**通过：XY 两方向/RZ 轴向、正反流、uniform/AMR、旧纯液/纯气、边界纯液/纯气/逐面交替、真实时间切分、共享 owner 生命周期。完整面集合实际经过 I12a3/I09，验证严格容量、纯相保持及边界 ledger。
- 内部/出流结果与 I12a7 的四个体积和全部 sweep diagnostics 逐字段相同。24 项旧测试完整输出逐字不变，architecture 保持 12 targets / 380 TUs / 899 files / 3085 includes / 5 legacy edges。
- 缺失/混源 BC、非法相态/面/时间/mesh、错误旧源或速度/特征 owner、opaque 或未支持速度/覆盖、moved-from source 均明确拒绝。原 M 未改。
- 81 冻结/参考、3 PDF 和 I12a7 archive 哈希不变。旧 crossing probe 仍退出 1，完整输出不变；没有将代数链路验证当作完整 benchmark 通过。

原始构建、CTest、旧输出比较、scope/冻结审计、runtime 与 archive 校验保存在 `runs/baselines/I12a8_20261001/`。运行 `quantum_interface_constant_inflow_contract.exe` 复现完整合同，`--guards` 仅检查拒绝与身份。

## 后续

接入由同一旧状态和同一速度/特征/BC owner 驱动的半步、整步，经过既有 I10/I11 后再接受状态。先验证常速度完整平移；一般速度、时变边界和完整 I12 仍待实际验收。
