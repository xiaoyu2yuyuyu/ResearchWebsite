# 分子动力学科研学习笔记

这是一个以 Markdown 为主要内容源的个人科研学习网站，使用 MkDocs 和 Material for MkDocs 构建。

## 平时主要看哪里

- `docs/`：网站正文。日常新增和修改笔记，主要在这里进行。
- `lammps-cases/`：以后存放小型、必要、可复现的 LAMMPS 输入脚本和分析脚本。
- `_网站配置/`：网站构建配置，通常不需要手工修改。
- `.github/workflows/`：GitHub Pages 自动部署所需的固定目录，通常不需要手工修改。

## 最简单的更新方式

1. 在 `docs/` 中新建或修改 Markdown 文件。
2. 在本地检查内容。
3. 使用 Git 记录修改并推送到 GitHub。
4. GitHub Pages 自动重新构建网站。

## 当前地址

- 在线网站：https://xiaoyu2yuyuyu.github.io/ResearchWebsite/
- GitHub 仓库：https://github.com/xiaoyu2yuyuyu/ResearchWebsite

## 文件收录原则

适合纳入仓库：

- Markdown 学习笔记
- 笔记实际使用的少量图片
- LAMMPS 输入脚本和小型配置文件
- 自己编写的 Python 分析脚本
- 小型、关键、难以替代的复现结果

通常不纳入仓库：

- dump、轨迹、restart 和大型日志
- 可重新生成的中间结果
- 虚拟环境、缓存和临时文件
- 密码、访问令牌和未确认可公开的资料
