# Characteristic 回归输入

`near_tangent_native.csv` 是均匀8×8 native闭域速度（144面），来自H13真实载荷负谱分支；只作为短契约输入，不是生产算例或解析参考。行号绑定原native布局，测试严格检查完整行序。

原程序在透明C1内部连接的近切向穿越中，返回Jacobian偏离同一点映射的有限差分；该错误直接破坏H08a必经FirstJet更新。契约在三个速度倍率上沿用原FD 5e-7门槛。原输入、程序、独立复算和完整失败见[反例说明](../../../../docs/task/H13_H08a_Transparent_Picard_Jacobian_Counterexample.md)。

入口：`cch_characteristic_contract`；单独定位使用同一测试程序参数 `--transparent-picard-regression`。不读取cases名称，不向生产算法传递参考答案。
