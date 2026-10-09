# I12a11r2 / 闭合工作精度与库存预算衔接报告

**FAIL：30 项回归通过，但原 XY AMR n=32 配置在首步末新增纯相残差拒绝。** 按 [TODO](I12a11r2_Closure_Working_Precision_TODO.md) 的失败规则保留代码和证据，继续独立根因诊断；不宣称原第 4 步库存问题已解决，I12a11/I12/Gate I 仍 FAIL。

## 实现与已通过验证

4 production + 2 tests。原库存 64 eps 常量移至共享私有 header，值与装配运算顺序不变。I10 同一 bubble 方程的工作残差改为 `min(relativeVolumeTolerance,64 eps/4)`，默认 16 eps；初始跳过、纯相、括区、求根及最终核验一致使用。Report 记录实际精度，M、first jets、B/J、几何证书和工作上限未改。详见 [API](../api/I10_Hermite_Gals_Bubble_API.md)。该工作精度分配不是整体求积误差界证明。

首轮 bubble 测试 70.12 秒 PASS。原 XY/RZ、uniform/AMR、定向扩区制造解的最终残差均在实际工作精度内；定向最坏残差 3.5527e-15。Newton 步数 388→503，二分仍 128，endpoint 2→0，未出现求根停滞。纯相预算内 4 eps 尾量原样保留，原 1e-12 尾量现在拒绝；M 与 identity 隔离断言通过。

完整 Release 重建 32.16 秒，无 warning；30/30 CTest 通过，207.91 秒。26 项输出逐字不变。其余变化为：

- architecture：900 文件、3088 includes；比原来新增 1 header、2 includes，12 targets/380 translation units/5 legacy edges 不变。
- bubble：上述更低残差、额外迭代及新增纯相拒绝验证。
- publication：只改变 TRANSACTION 残差，4.6045e-11→1.1102e-15，其余整段输出不变。
- translation smoke：8 组全部完成；4 组静止和 RZ uniform 平移输出不变，另 3 组几何误差变化。XY AMR 曲率误差由 1.1013e-13 增至 1.6354e-10；volume/ledger/reflux/drift/cell-M 等列不变。更严格体积求解不等于曲率误差必然改善。

## 原配置的新失败

`--translation-inventory-audit` 现要求真实完成原 4 步。105.43 秒后退出 2，首次接受时间仍为 t=0：

```text
end_closure
[Interface/HermiteGalsBubble primary=8 cell=1513 stage=topology/amplitude guard]
pure-cell mismatch would create an unmarked component
```

半步与全步库存装配都已返回，尚未接受第 1 步末状态；不能拿这种更早拒绝作为原 cell 1472 第 4 步修复成功。accepted identity/M 未变。下一项需量化 cell 1513 的目标、纯相证书量、旧 M/几何量和同源扫掠差额，区分工作精度分配、体积测量与空间重构的影响；目前尚无根因结论。

按先决验证失败即保留证据并拆项的规则，本次未继续 40 步多项式矩阵或 12 组非多项式矩阵，明确记为 NOT_RUN。旧 midpoint probe 仍退出 1，输出不变。81 冻结/参考文件、3 PDF、前序 archive 哈希保持。源码、28 可执行与运行库、原始日志及逐项审计归档于 `runs/baselines/I12a11r2_20261001/`；归档完整性 PASS 与本项数值 FAIL 分开记录。
