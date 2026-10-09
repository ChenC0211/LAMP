# H13 / 已授权新source的可信参考消费者同步

状态：2026-10-09 / 081136，**已获联合方案授权，两组reference永久同步，标准重测待完成**。仍属H13，未修改生产source/characteristic或原门槛。以下独立预算请求为历史，现纳入[联合10份已有契约预算](H13_Phase_Cut_Kinetic_Hodge_Method_Change.md)。

## 唯一必要的数据同步

已授权native面矩准备改变了实际共同source。现有 `tests/contracts/hydrodynamics/cch_characteristic_contract.cpp` 的两组固定reference仍属于旧P2顶点source，导致合法的新场被拿去比较旧流映射。原测试分组、时间序列、FD-J、C1及精度断言均有效，须保留；只同步这两行reference数据与来源注释。

当前[H13联合提案](H13_Native_Flux_Port_Method_Change.md)将新增验证限制为regular/traction两份已有契约，现有工程默认任务预算也限两份测试文件。因此本同步请求增加第三份已有测试文件的**reference数据维护**，不扩展测试分组/生产模型，不由此授权修改characteristic/Jacobian公式。该预算例外独立于第五production文件的[KEP补充提案](H13_Kep_Preparation_Accuracy_Amendment.md)，未获明确范围补充前只保留隔离候选。

## 参考如何取得及误差定义

程序 `runs/checks/H13/20261009_003903/characteristic_coefficients.cpp` 导出原native输入所生成的实际Q4 tensor数据、实际subface bubble系数、分片、query jet与40/80/160/2560步map。独立 `reference_taylor.py` 从二进制系数解释的真实分片多项式重组Ux/Uy；没有从production tracer生成真值或修改接受状态。

对Ux采用100位多精度30阶Taylor和解析一阶变分；逐真实y分片推进，不跨错x support piece。复Picard圆盘满足 `r<=min(Ry/2,Rx/(2*M),0.25/Lx)`，真实步长 `h<=r/4`；Cauchy余项按 `B*(h/r)^(p+1)/(1-h/r)` 累计。实际Uy-1、其梯度与速度Hessian由分片多项式系数范数界定，用Gronwall界传播到完整二维point/J。多精度算术误差另列，原double转换4e-16 allowance仍保留。

原始多项式对actual query jet核验最大绝对差2.95347e-14/5.52395e-14，单独记录数值evaluation差。参考自身的不确定度是数学多项式连续流的误差界；actual离散map与它的差包含原tracer时间/evaluation误差，并继续按原断言验收。

| 原witness | reference point界 | reference J界 | 2560步point误差 | 2560步J误差 | 40/80/160末档阶 |
| --- | --- | --- | --- | --- | --- |
| uniform | 1.23104e-16 | 1.42503e-14 | 2.72621e-14 | 7.27204e-13 | 3.033 |
| internal | 1.50788e-16 | 1.26083e-14 | 8.79227e-15 | 1.13059e-13 | 2.689 |

参考界仍<=1e-12；原point/J 1e-11、时间阶>=1.8、FD-J 5e-7及C1/event门槛全保留。永久double字面量的时间阶因点值舍入稍异，暂存完整原契约实测3.03385/2.69974，仍通过。

## 可复算控制及完整测试

- 用旧source语义控制运行同一独立参考，其point/J转换到double与原两组可信reference逐位相同；原连续流uncertainty预算覆盖差异，旧完整characteristic契约也退出0。
- 30/34阶的reference差最大约1.16e-57/8.65e-54，处于解析余项界内；这项agreement只作控制，不替代上述误差界。
- 污染一个原系数，query一致性明确拒绝；将实际map的x污染1e-4，reference本身保持逐字不变、原精度gate明确失败。参考没有使用待验证的map当真值。
- 暂存 `characteristic_reference_candidate_contract.cpp` 只替换原两个reference tuple；**原完整契约退出0，MOTION_RESULT,0,PASS**。原24时间族、FD、C1、FirstJet skeleton、透明Picard等分组全部保留，没有少跑或放宽阈值。

暂存数据为：

```cpp
const double refs[2][8]={
 {0.43436894113158286,0.41137,1.0008305859651514,-0.024193325947587913,0,1,1.231042115119329e-16,1.4250292077736276e-14},
 {0.4369693706869503,0.36137,0.9999910118882219,-0.03035256960438956,0,1,1.5078805650408715e-16,1.2608319572785833e-14}
};
```

## 范围、可比性与采用前验收

获范围补充后，仅维护上述一个已有测试的两组reference和来源说明，标准Release构建及原完整契约复测。保存旧reference及其旧source语义结果；不把旧数值用于新场、不改判旧结果，也不要求新map逐位恢复旧场。

本同步不改任何数值内核、未知量、accepted state、F/R/uR/KEP/M/C、source准备、AMR或physical gate，不是独立r11精修任务。参考更新不证明静滴、毛细波、低黏或H13整体通过；那些原验收仍须完成。
