# I12a15r3 / 权威 Q4 的相关 pullback 证书

见 [TODO](I12a15r3_Authoritative_Q4_Affine_Pullback_TODO.md) 与 [r2设计](I12_Affine_Face_Characteristic_Design.md)。新接口保留原snapshot权威，只适配同一多项式到相同q参数下的affine范围。

## 权威系数

`nativeQuarticLocalDofs` 是 rounded display数据，不能重新定义phase。实际source保存 `local.numerator[25] / local.denominator` 的精确齐次数据，包含shared/hanging/sector约束。已有 `QuarticCell::bernsteinControlsCleared` 返回与原多项式严格相同的张量Bernstein分子，每轴清除12，因此公共正分母为144*原denominator。

在每个轴上用精确Expansion计算幂系数 `a_k=C(4,k)*sum((-1)^(k-j)*C(k,j)*b_j,j=0..k)`，再作用第二轴。有限差分中的精确零保持零；没有拟合、截断、小量删除或浮点cache代用。两轴变换是同一Q4的代数恒等式，不增加新的geometry DOF。共同幂二平衡不会改变商；每个分子/正分母分别outward包围，再以正分母两端除法给出物理phi单位的系数区间。

## 仿射复合

导入r2的两个physical forms，不改变其共享q或剩余带。原叶chart为 `(B_x-xlo)/(xhi-xlo)` 和相同的y式；cell宽度由真实double端点的实数差向外包围，不能把约束所用的真实宽度换成一个四舍五入常量。正区间除法逐affine系数/剩余项包围。

固定25项二维Horner复用既有strict affine乘加；它包含所有非线性项和足够的逐操作舍入带。输出仍是一个中心、两个一次参数系数和剩余区间。对原plane/Q2，高次系数若精确为零，就按原代数保持零；不会因Bernstein重复依赖人为生成高次多项式。

## 叶与来源边界

snapshot入口指定一个active native Q4叶，直接取得它的原齐次多项式。API允许传入范围比该叶更大的B forms，是为了完整closed donor覆盖时包围同一叶的多项式延拓；**只有真实B位于该叶的参数点才能把这个符号用于phase**。不从该接口选择donor，不把叶外延拓当作外部液体或physical BC。

后续消费者必须沿原branch box枚举所有相交active叶并证明closed覆盖；可以在各叶对应参数子集使用各自证书，或采用对所有候选叶均成立的共同符号子区域。不同叶的系数都恢复到同一phi单位，但仍不能平均重构或生成第二几何。inventory wl/wg仍属于同一旧source，M没有参与这个多项式求值。

本项不实现液体多边形面积或gross体积，也不替代a14现有strict sign入口。实际face flux、完整driver、时间阶、动态AMR及低黏/无黏仍按原gate推进。