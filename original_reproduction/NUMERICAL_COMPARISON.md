> **Historical audit, 2026-09-28.** The text below is preserved as written. “Current” and “this run” refer to that audit, and its equation numbers, source lines, page numbers and revision recommendations are historical. For present installation and results, use [README.md](README.md) and the package [MANUSCRIPT_MAPPING.md](../MANUSCRIPT_MAPPING.md). The original Eqs. (38)/(39) now correspond to Main Eqs. (49)/(50), with stable labels `eq:d4-example` and `eq:d6-example`. This file is not a receipt for a new run.
>
> **历史身份说明（2026-09-28）：**以下正文按当时记录原样保留。“当前”“本轮”和“本次”均指当时的审计，文中的公式编号、源码行号、页码及修改建议不代表现稿。现稿对应主文式（49）和（50）；安装与运行请以本目录 README 和软件包根目录的 MANUSCRIPT_MAPPING 为准。

# 原稿数值与本轮复现值

本轮实际执行保留原样的历史脚本，并另写独立实现复核。当前主文 Methods（TeX 1926–1936、PDF pp14–15）和 SI 的 numerical cross-checks（TeX 1373–1381）报告的对象与本包 Eq. (38)/(39) 输入一致。

| 例子 | 原主文/SI 残差 | 历史 README 示例值 | 本次原样脚本真实输出 | 本次独立实现接受块残差 |
|---|---:|---:|---:|---:|
| d=4, n=4, Eq. (38) | 1.8e-16 | 1.732e-16 | 1.732e-16 | 1.224647e-16 |
| d=6, n=3, Eq. (39) | 6.6e-16 | 6.523e-16 | 6.505e-16 | 8.329022e-16 |

这里的残差均为各颜色块中最大的矩阵 operator two-norm 差；独立实现比较的是“直接 D 共轭理想投影”和“从完整相位函数差分构造的分支效果”，历史脚本比较的是“循环生成元投影乘积”和“单项式条件相位分支效果”。二者数学上检验同一接受算符恒等式，但浮点运算路径不同。

**结论是实现恒等式和谱隙成功复现；不是原稿每个末位小数被逐字重现。** 所有量均在约 1e-16 至 1e-15 的数值误差范围。缺少旧运行环境及完整原始日志，不能断言原稿 1.8e-16 和 6.6e-16 的具体差异只由哪一种舍入或库版本造成，也不应强行调参数凑旧值。

本次 Python 3.12.14 / NumPy 2.3.5 / complex128 环境中，历史脚本的 d=3 回归残差为 5.680e-16，而历史 README 为 5.473e-16；这也是末位数值环境差异，谱隙与全部检查结论不变。原样脚本自身断言容差为 1e-10，新独立入口使用 1e-12，均通过。

独立完整联合测量检查还得到：

- d=4 全部联合扇区效果的最大残差 1.731912e-16；d=6 为 1.097546e-15。
- 所有本地/全局基正交性、完备性、接受块幂等性、Hermiticity、目标接受与谱检查中，最大误差约 1.999e-15。
- d=4 验证谱 {0, 1/3, 2/3, 1} 的重数分别为 {135, 99, 21, 1}，谱隙为 1/3。
- d=6 验证谱 {0, 1/2, 1} 的重数分别为 {175, 40, 1}，谱隙为 1/2。

主文与 SI 修订时建议统一报告本次真实输出和所用方法，或报告清楚容差下的残差上界，并把完整环境/原始输出放入补充软件。不要同时保留互不说明的旧日志值、新日志值和另一实现的值。

本包没有修改主文或 SI，未自动替换原数值。它仅闭合有限原稿例子及 d=2,3 回归的本地可复现性；不代表所有理论、所有 SI 基准或新六维路径家族已复现。
