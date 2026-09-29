# 第三方声明

PhyloSlicer 是一份独立实现：仓库内不含任何派生自其他软件包的源代码，也不打包、
再分发或链接任何第三方代码。仓库内确实携带的、由第三方软件*生成*的文件，在此如
实登记。

本文件是 [`NOTICE.md`](NOTICE.md) 的中文译文，正文以英文文件为准。

## treesliceR

`validation/reference/*.csv` 与 `validation/reference/*.tre` 是参考实现
**treesliceR 1.1.0**（R 包）——即本包所重新实现的时间切片范式——的输出：

> treesliceR — MIT License，© 2023 treesliceR authors（见该包自带 `LICENSE`：
> `YEAR: 2023`、`COPYRIGHT HOLDER: treesliceR authors`）。
> 作者：Matheus Lima Araujo（aut、cre、cph；ORCID 0000-0002-9111-725X）、
> Luiz Gabriel Souza e Souza Ferreira（ORCID 0009-0002-5881-9791）、
> Gabriel Nakamura（ORCID 0000-0002-5144-5312）、
> Marco Tulio Pacheco Coelho（ORCID 0000-0002-7831-3053）、
> Thiago Fernando Rangel（ORCID 0000-0002-2001-7382）。
> 软件包：<https://github.com/AraujoMat/treesliceR> ·
> 论文：<https://doi.org/10.1111/ecog.07364>

这些文件是数值交叉验证的**期望值**一侧，而不是依赖项：PhyloSlicer 无需安装 R 即
可运行与安装，`validation/compare.py` 只读取这些 CSV。`benchmarks/*_r.R` 会调用已
安装的 treesliceR 来复现比对表格，它们是驱动脚本，其中不含任何 treesliceR 代码。

由于参考值由另一个浮点路径不同的外部程序产生，
[`validation/NOTES_CN.md`](validation/NOTES_CN.md) 中的容差本身就是这批资产的一
部分：一条速率行在 1e-4 上相符，说的是两个不同的优化器求解同一个目标函数，而不是
说二者逐位相同。

## 运行时依赖

必需依赖为 numpy、scipy、pandas、matplotlib；dendropy、biopython、geopandas、
numba 为可选扩展。PhyloSlicer 在调用时才导入它们，不做任何拷贝，并且即使缺少全部
可选依赖，功能依然可用。

## 文档许可

`README*.md`、`docs/USER_GUIDE*.md`、`benchmarks/README*.md`、
`validation/NOTES*.md`、`testing_suite/*.md` 的文字部分以 CC-BY 4.0 发布；本仓库
的代码以 BSD-3-Clause 发布（见 [`LICENSE`](LICENSE)）。
