# I12a10 / 完整平面跨面发布报告

**FAIL：48 组中 36 组通过，12 组 AMR 正向平移在中点外边界对称性检查拒绝。** 原 q=1/4、3/4 的时间跨面库存缺口已在通过的完整链路中消除，但不能以这些通过结果替代剩余失败。

## 实现与成功范围

遵守 [TODO](I12a10_Planar_Crossing_Publication_TODO.md) 的 0 production + 2 test 文件预算。制造解 Case 增加明确平面与 XY 两个常平移方向；初始 I07/I10/I11 后完全复用原 advanceTranslation。外部 transverse symmetry 声明一致传给初始化、region、bubble 和 publication；实际入流纯液/纯气历史对完整窗口验证同号，并共享原 immutable BC owner。所有 interior first jets 和 center 仍由原旧 scalar/B/J 获得。

`--crossing-publication` 原 n=4、dt=.02、速度幅值 1，覆盖 XY x/y、RZ z、uniform/AMR、正反速度/相及 q=1/4、3/4。36 组成功发布真实半步与整步，严格库存有界，穿面四个 gross volumes 与独立解析事件液量一致。最大逐单元相对库存误差 1.15e-15、穿面相对误差 3.10e-15、first-jet/法向误差 5.14e-16、曲率绝对误差 3.74e-14 以下。

成功组合为全部 24 个 uniform 案例，以及 AMR 的全部 12 个负向平移案例；不是只验证单面代数或解析 end skeleton。原 accepted identity/M 在成功和失败调用之后均逐项不变，只有成功返回的 publication 可供上层接受。

## 拒绝与后续

12 个失败均为 AMR、正速度，覆盖三个方向、两种相位、两个跨面时刻。stage 为 midpoint_closure，原错误是 `symmetry section requires identically zero scalar normal derivative`，发生在 predictor measure 的完整边界检查。未提高容限、删除 symmetry 声明、修改 M 或将解析解写入 interior；原状态保持。

错误中的 cell=27 是 bubble 初始化后的诊断上下文，尚未证明它就是实际失去对称性的边界单元。下一项需要逐项检查 transported boundary first jets、I06 mixed 数据及同源 Q4 法向导数，定位 exact symmetry 在哪一层丢失，然后只修复已证实的事实。不能把该失败归因于已修复的时间中点求积，也不能直接改用不同 BC 避开检查。

q=1/2 的 midpoint 整面零迹仍未由本项支持，继续独立验收。非恒定速度与其它 I12 gate 也保持未完成。

## 验证与交付

首次 API 拼写编译错误已修正，最终 Release 构建通过，3.48 秒，无 warning；原始错误日志保留。10 项 CTest 中 9 项通过，新增 crossing publication 失败，退出 8，约 29 秒。**9 项旧测试完整输出逐字不变**；architecture 保持 12 targets / 380 TUs / 899 files / 3085 includes / 5 legacy edges。

所有 production、其它 27 个可执行、81 冻结/参考、3 PDF 与 I12a9r2 archive 哈希不变。旧 midpoint crossing probe 退出 1 且输出不变。第一次 CTest 已运行完，但 PowerShell 将预期失败的 stderr 当终止异常，未写结果 JSON；其原始日志另存，修正 runner 的退出码捕获后重新验证，未将该数值失败设为 WILL_FAIL。证据及逐项交付校验位于 `runs/baselines/I12a10_20261001/`。证据审计 PASS 不改变本项数值 FAIL；I12/Gate I 未通过，不进入 Hydro。
