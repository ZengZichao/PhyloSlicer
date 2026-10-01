# 安全策略

*本文件为 SECURITY.md 的中文对应版本；**以英文版（SECURITY.md）为准**。*

## 受支持的版本

PhyloSlicer 是科研软件包，只有最新的 minor 版本线会获得安全修复。

| 版本   | 是否受支持 |
| ------ | ---------- |
| 0.1.x  | ✅         |

## 报告安全漏洞

**首选渠道：** GitHub 私密漏洞报告 —— 仓库页 **Security → Report a
vulnerability**，仅维护者可见。

**备选渠道：** 邮件 <zengzichao@sjtu.edu.cn>（标题请注明
"PhyloSlicer security"）。

任何疑似安全问题**请不要**开公开 issue。

如可能请附上：

- PhyloSlicer 版本（`python -c "import phyloslicer; print(phyloslicer.__version__)"`）
- Python 版本与操作系统
- 最小复现（命令、输入文件、完整报错栈）
- 你的影响评估（数据损坏？结果错误？代码执行？）

## 处理流程

- **7 天内**确认收到
- 与你协作修复或缓解，并致谢（可选择匿名）
- 修复后发版，并在 CHANGELOG 中以公告形式记录

## 范围说明

PhyloSlicer 处理本地的系统发育与空间数据文件。重点关注：`[io]`
适配器中的路径处理、树/位点文件的不安全反序列化，以及任何可能因构造
输入文件而执行代码的路径。
