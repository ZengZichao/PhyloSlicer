# Security Policy / 安全策略

## Supported versions / 受支持的版本

PhyloSlicer is a research package; only the latest minor line receives
security fixes.

PhyloSlicer 是科研软件包，只有最新的 minor 版本线会获得安全修复。

| Version | Supported / 受支持 |
| ------- | ------------------ |
| 0.1.x   | ✅                 |

## Reporting a vulnerability / 报告安全漏洞

**Preferred channel / 首选渠道:** GitHub 私密漏洞报告 —— 仓库页
**Security → Report a vulnerability**，仅维护者可见。

**Alternative / 备选渠道:** email <zengzichao@sjtu.edu.cn>
（邮件标题请注明 "PhyloSlicer security"）。

Please **do not open a public GitHub issue** for anything you suspect to be a
security problem.

任何疑似安全问题请**不要**开公开 issue。

Please include, where possible / 如可能请附上：

- the PhyloSlicer version (`python -c "import phyloslicer; print(phyloslicer.__version__)"`) / 版本号
- Python version and operating system / Python 版本与操作系统
- a minimal reproduction (command, input file, stack trace) / 最小复现（命令、输入文件、报错栈）
- your assessment of impact (data corruption? wrong results? code execution?) / 影响评估

## What to expect / 处理流程

- Acknowledgement within **7 days** / **7 天内**确认收到
- A fix or mitigation coordinated with you, with credit unless you prefer
  anonymity / 与你协作修复并致谢（可选择匿名）
- A release + CVE-style advisory entry in the CHANGELOG once fixed / 修复后发版并在 CHANGELOG 记录

## Scope notes / 范围说明

PhyloSlicer processes local phylogenetic and spatial data files. Areas of
particular interest: path handling in the `[io]` adapters, unsafe deserialisation
of tree/site files, and any code path where a crafted input file executes code.

PhyloSlicer 处理本地的系统发育与空间数据文件。重点关注：`[io]` 适配器中的
路径处理、树/位点文件的不安全反序列化，以及任何可能因构造输入文件而执行代码的路径。
