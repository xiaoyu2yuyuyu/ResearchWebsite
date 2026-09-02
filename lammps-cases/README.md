# LAMMPS 复现文件

这个目录以后用于保存体积小、确实有复现价值的 LAMMPS 文件，而不是保存完整模拟输出。

建议每个案例使用一个独立子目录，例如：

```text
lammps-cases/
└─ melting-example/
   ├─ README.md
   ├─ in.melting
   └─ analyze.py
```

适合保存输入脚本、必要的小型数据、分析脚本和少量关键结果。轨迹、restart、大型日志和可重新生成的中间文件会被忽略。

