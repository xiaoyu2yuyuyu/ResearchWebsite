# LAMMPS Tutorial 4 — Nanosheared Electrolyte

> 文档职责：Tutorial 4 的唯一案例主文档，集中保存案例事实、复现证据、代码教学、理解检查和关闭边界。  


## 第一部分：案例事实、运行结果与验收

### 1. 概况

本教程模拟两面铝墙之间的水和 NaCl 电解质。体系首先被创建并平衡，然后上下墙沿 `x` 方向以相反速度运动，使中间的电解质产生剪切流动。最终关注水和离子沿 `z` 方向的密度分布、沿 `x` 方向的速度剖面、墙面切向力以及由此估算的黏度。

---

### 2. 相对官方文件的修改

教学版保留了官方教程的盒子、粒子数、力场、平衡温度、时间步长、剪切速度、运行步数和空间分箱参数。主要修改是把服务器不方便使用的 LAMMPS 图像输出换成 OVITO 可以读取的轨迹，并增加便于学习的诊断信息。

| 文件 | 教学版修改 | 修改理由 | 是否改变物理模型 |
|---|---|---|---|
| `create.lmp` | 用 `dump custom` 输出 `trajectory-create.lammpstrj`，代替 PPM 图像；增加 `thermo_style` | 服务器没有图形界面，需要在本机用 OVITO 查看；同时方便检查初始能量和压力 | 否，只改变输出 |
| `equilibrate.lmp` | 分别输出最小化轨迹和平衡轨迹；最小化后使用 `undump` 关闭第一条轨迹 | 区分能量最小化和 NVT 平衡两个阶段，便于学习和检查 | 否，只改变输出 |
| `shearing.lmp` | 输出 OVITO 轨迹；在日志中明确打印流体温度、墙温度和两面墙的 `x` 方向受力；增加 `shearing-final.data` | 便于观察剪切、诊断恒温和后续继续计算 | 否，只改变输出和诊断 |
| 三个轨迹 | 同时输出数字 `type` 和文本 `element`，并使用 `sort id` | OVITO 自动显示 O、H、Na、Cl、Al；保持每帧原子顺序稳定 | 否 |

轨迹中使用：

```lammps
id mol type element q x y z ix iy iz vx vy vz
```

- `type` 保留 LAMMPS 的数字原子类型，方便程序分析；
- `element` 由 `dump_modify ... element O H Na Cl Al` 映射为元素名称，方便 OVITO 识别；
- `sort id` 只让原子按 ID 排序，不是模拟成功的必要条件；
- `x y z` 是盒子内的坐标，`ix iy iz` 是周期映像计数，可用于重建跨越周期边界后的连续位置。

剪切运行完成后，服务器上的 `run-shearing.sh` 被精简为只生成一个 `run-info-shearing.txt`，不再生成 `metadata`、`time`、`status` 和重复的 `screen` TXT。已完成运行所用旧脚本的 SHA256 保存在 run-info 中，本机原版脚本保持不变；本次剪切日志保存为 `log-shearing.lammps`。

---

### 3. 当前案例中的输入文件作用

| 文件 | 在本案例中的作用 | 关键产物 |
|---|---|---|
| create.lmp | 定义盒子、上下铝墙、水和 NaCl，设定力场与分组；只用 run 0 检查初始构型。 | create.data、初始轨迹与日志 |
| equilibrate.lmp | 读取初始结构，先最小化消除严重不合理接触，再在 NVT 下平衡 30000 步。 | equilibrate.data、最小化/平衡轨迹与日志 |
| shearing.lmp | 读取平衡结构；给上下墙相反的 x 速度；对流体与墙分别控温；按 z 分箱并对后期窗口平均。 | 剪切轨迹、water/wall/ions 剖面、最终 data 与日志 |

### 4. 本案例的日志诊断与数值判断

#### 4.1 create 阶段实际警告

##### `Ignoring 'compress yes' for molecular system`

删除原子后，LAMMPS 默认可能尝试压缩原子 ID；分子体系涉及键角拓扑，因此忽略这种压缩。这里的 `compress yes` 不是脚本中 `fraction 0.15 yes` 的那个 `yes`：后者表示尽量精确选择 15% 的候选原子，前者是删除后重新编号 ID 的行为。

本次命令是：

```lammps
delete_atoms random fraction 0.15 yes H2O NULL 482793 mol yes
```

`mol yes` 才是本次安全性的关键：若随机选中某个水分子的任一原子，LAMMPS 会删除这个水分子的全部原子，避免留下缺氧或缺氢的残缺水分子。原子 ID 有空缺并不改变坐标、元素、键或角；本次应检查的是水分子拓扑是否完整，以及后续是否出现 `missing atoms`、`Lost atoms` 或键角错误，而不是 ID 是否连续。

如果后续脚本错误地假设原子 ID 必须连续，则需要额外留意；LAMMPS 本身并不要求原子 ID 连续。

##### `No fixes with time integration, atoms won't move`

create 阶段只有 `run 0`，本来就不希望原子移动，所以这是符合预期的提醒。

如果同样信息出现在本应运行 NVT 的平衡或剪切阶段，则说明没有时间积分 fix，属于必须修正的问题。

##### `Communication cutoff ... shorter than a bond length based estimate`

LAMMPS 先提醒当前通信范围可能不足以覆盖分子拓扑，随后又打印：

```text
Increasing communication cutoff to 15.1118 for TIP4P pair style
```

说明 TIP4P 样式已经自动把通信范围提高到需要的值。本次没有 Lost atoms，后续也正常，因此可以接受。判断重点不是只看第一条 warning，而是看软件是否进行了合理修正以及后续是否出现真正异常。

#### 4.2 equilibrate 阶段实际警告

##### `Using fix shake with minimization`

SHAKE 在正常 MD 中会在每个时间步把 O–H 键长和 H–O–H 角拉回目标几何；但最小化没有真实时间步，不能执行这种严格的“移动后修正”。因此本阶段的

```lammps
fix myshk H2O shake 1.0e-5 200 0 b O-H a H-O-H kbond 2000
```

会把严格 SHAKE 近似为 `kbond 2000` 所定义的强谐和弹簧约束。最小化后脚本重新定义 SHAKE，再进入 NVT；此时才使用正常 MD 的几何约束。故这条警告表示 LAMMPS 正在采用预期的近似处理，不表示 SHAKE 失败。

需要进一步处理的情形包括：最小化后出现 NaN 或 Lost atoms、SHAKE 无法收敛、出现 `Shake determinant < 0`，或强约束导致最小化明显不稳定。本次未发生这些现象；但最小化因最大力评估次数停止，所以只能说明“足以继续教学计算”，不能声称严格收敛。

##### `Stopping criterion = max force evaluations`

这不是 WARNING，而是最小化停止原因。它说明力计算次数达到 1000 的上限。判断时应同时查看：

- 初始和最终能量；
- 初始和最终力；
- 后续动力学是否稳定；
- 当前目标是教学继续还是科研严格收敛。

本次只能得出“明显改善并可继续教学运行”，不能写成“最小化已经严格收敛”。

#### 4.3 本次平衡阶段的综合判断

```text
程序退出：正常
ERROR：无
NaN/Inf：无
Lost atoms：无
Dangerous builds：两段均为 0
最终步数：30000
温度：总体围绕 300 K 波动
能量：没有爆炸或持续发散
墙间距：仍有明显调整和振荡，未证明稳定平台
```

因此结论分两层：

1. **技术结论**：平衡阶段运行成功，可以继续 Tutorial 4 的教学剪切；
2. **科研结论**：30 ps 数据不足以证明充分平衡或统计收敛。

#### 4.4 shearing 阶段的温度与压力相关 warning

剪切日志出现两次：

```text
Temperature for fix modify is not for group all
```

以及一次：

```text
Temperature for thermo pressure is not for group all
```

它们警告的是**统计组不一致**，不是在警告 `x` 方向没有恒温。剪切脚本中：

```lammps
compute Tfluid fluid temp/partial 0 1 1
fix mynvt1 fluid nvt temp 300 300 100
fix_modify mynvt1 temp Tfluid

compute Twall wall temp/partial 0 1 1
fix mynvt2 wall nvt temp 300 300 100
fix_modify mynvt2 temp Twall
thermo_modify temp Tfluid
```

`temp/partial 0 1 1` 故意不把 x 方向的平均剪切流速当成热运动，以免恒温器错误地削弱剪切流。`Tfluid` 只统计流体，`Twall` 只统计墙；它们分别用于对应的 NVT 恒温器是合理的。

警告的对象是全体系压力：压力的位力部分默认来自 `all`，而若动能部分来自只含 fluid 或 wall 的局部温度，二者的统计范围不一致。因此：

- `Tfluid` 只包含 fluid，`Twall` 只包含 wall；
- LAMMPS 的某些压力计算默认按 all 处理；
- 如果把局部温度的动能贡献和全体系的位力贡献组合成压力，统计范围可能不一致。

本次两个恒温器的 fix 组与温度 compute 组分别匹配，且使用 NVT、不用压力调节盒子，thermo 也没有输出 `press`，所以这些 warning 不影响本次温度控制、墙力和速度/密度统计。但是以后若要分析严格的全体系压力、应力张量，或使用 NPT 控压，必须另定义覆盖 `all` 且与压力统计一致的温度 compute，不能直接使用当前的 `Tfluid`。

另两条通信范围 warning 与前面阶段相同：LAMMPS 发现初始通信 cutoff 偏短，随后自动提高到 15.1118 Å；本次没有 Lost atoms 且 `Dangerous builds = 0`，可以接受。

#### 4.5 本次剪切阶段的综合判断

```text
程序退出：正常，exit_code=0
最终步数：200000
热力学时刻：801
轨迹：201 帧，Step 0 到 200000，每帧 1680 个原子
ERROR：无
NaN/Inf：无
Lost atoms：无
Dangerous builds：0
流体温度：300.29 ± 10.23 K（全部801个时刻）
墙温度：299.22 ± 14.70 K（全部801个时刻）
后约150 ps上墙平均力：-2.735 ± 7.509 kcal/(mol·Å)
后约150 ps下墙平均力：+2.678 ± 7.990 kcal/(mol·Å)
两墙平均合力：-0.0567 kcal/(mol·Å)
```

上墙向 `+x` 运动却受到平均 `-x` 力，下墙向 `-x` 运动却受到平均 `+x` 力；水和离子的下半区域平均向负方向运动，上半区域平均向正方向运动。墙速、墙面反作用力和流体速度梯度三者方向一致，共同说明剪切动量已经传递给流体。

墙力的瞬时波动远大于平均值，因此当前只能确认定性剪切趋势，仍需画出墙力时间序列、分块平均并估计不确定性，才能讨论科研统计是否稳定。

上面的“±”目前表示热力学采样值的标准差，只描述瞬时波动大小，不是独立重复模拟得到的误差棒，也不能直接当作平均值的不确定性。

---

### 5. 数据、绘图和论文价值

#### 5.1 当前已经完成的数据链

本次已完成四张图：

1. `equilibrate-pressure-wall-distance-publication`：短平衡的诊断图；
2. `shearing-density-velocity-publication`：教程核心的结构与剪切流结果图；
3. `shearing-wall-force-publication`：剪切应力传递与时间分块诊断图；
4. `shearing-temperature-diagnostics`：恒温与运行健康检查图。


| 图文件名 | 直接数据来源 | 绘图脚本与处理 | 能回答的问题 | 当前判断 |
|---|---|---|---|---|
| `equilibrate-pressure-wall-distance-publication` | `extracted_thermo/log-equilibrate_thermo_02.dat`；该文件从 `log-equilibrate.lammps` 的 NVT thermo 数据块提取，列为 `Step`、`Press`、`v_deltaz` | `plot-pressure-wall-distance.py`；Step 转 ps，压力转 katm，距离从 Å 转 nm，10 点移动平均 | 压力和墙间距是否仍在明显调整？ | 程序正常运行并发生松弛；不能单独证明科研收敛。 |
| `shearing-density-velocity-publication` | `shearing-water.dat`、`shearing-wall.dat`、`shearing-ions.dat` | `plot-shearing-density-velocity.py`；沿 z 的质量密度和 x 方向速度剖面；速度图排除 `Ncount < 0.2` 的低占有率分箱，并记录筛选规则 | 水、离子和墙怎样沿孔道分布？剪切是否传到流体？ | 水在墙附近出现分层；水速度随 z 近似线性变化，符合 Couette 流的定性预期。离子只有 30 个，不能据此作强吸附或精确输运结论。 |
| `shearing-wall-force-publication` | `log-shearing.lammps` 的 `Step`、`f_mysf1[1]`、`f_mysf2[1]` | `plot-shearing-wall-force.py`；原始墙力、10 ps 移动平均、后期 150 ps 累计平均和六个 25 ps 非重叠时间块 | 上下墙是否传递方向相反、大小接近的剪切应力？ | 两面墙平均受力方向相反，大小相差约 1.82%，支持剪切应力传递；但六个时间块仍明显波动，不能证明统计收敛。 |
| `shearing-temperature-diagnostics` | `log-shearing.lammps` 的 `Step`、`c_Tfluid`、`c_Twall` | `plot-shearing-temperature.py`；全程温度、约 10 ps 移动平均、后期六个 25 ps 时间块平均 | 剪切时恒温是否明显失效？ | 后期流体为 `300.33 K`、墙为 `299.23 K`，没有持续升温；这是运行健康检查，不是黏度收敛证明。 |

三份 `.dat` 都只有一个第 200000 步的时间平均数据块，对应第 50010 到 200000 步的 15000 次采样。它们足以支持本教程的单个最终剖面，不能自身比较多个独立时间窗口。

#### 5.2 剪切应力与表观黏度：数据来源和计算规则

| 量 | 数值 | 来源或计算 |
|---|---:|---|
| 盒子横向长度 | `Lx = Ly = 24.24 Å` | 来源：`create.lmp` 的 `lattice fcc 4.04` 与 `region box block -3 3 -3 3 -5 5`。x、y 的跨度均为 6 个晶格单位：`Lx = Ly = 6 × 4.04 = 24.24 Å`。 |
| 每面墙面积 | `A = 5.8758 × 10⁻¹⁸ m²` | 计算：`A = Lx × Ly = 24.24² = 587.5776 Å² = 5.8758 × 10⁻¹⁸ m²`。 |
| 每面墙平均力大小 | `F̄wall = 2.7027 kcal mol⁻¹ Å⁻¹` | 原始数据：`log-shearing.lammps` 的 `f_mysf1[1]`、`f_mysf2[1]`。计算：`plot-shearing-wall-force.py` 取 50–200 ps，先把下墙力反号，再合并两墙所有 thermo 样本求算术平均；公式为 `F̄signed = mean(Ftop, −Fbottom)`，本次 `|F̄signed| = 2.7027`。 |
| 水速度梯度 | `γ̇ = 2.014 × 10¹⁰ s⁻¹` | 原始数据：`shearing-water.dat`。计算：`plot-shearing-density-velocity.py` 保留 `Ncount ≥ 0.2` 的水分箱，按 Ncount 加权拟合 `vx = az + b`；得到 `a = 20.14 m s⁻¹ nm⁻¹`，并换算 `γ̇ = a × 10⁹ s⁻¹`。 |
| 单墙剪切应力 | `τ = 31.96 MPa` | 计算：`τ = |F̄wall| / A`。力单位换算为 `1 kcal mol⁻¹ Å⁻¹ = 6.9477 × 10⁻¹¹ N`，因此 `τ = (2.7027 × 6.9477 × 10⁻¹¹) / (5.8758 × 10⁻¹⁸) Pa`。 |
| 单墙定义下的表观黏度 | `ηapp = 1.59 mPa·s` | 计算：`ηapp = τ / γ̇ = 31.96 × 10⁶ / (2.014 × 10¹⁰) Pa·s`。 |

上下墙是同一片流体两侧对同一剪切应力的两次测量，因此应使用

```text
τ = (|Ftop| + |Fbottom|) / (2A)
```

而不是把两面墙的力相加后仍除以单面墙面积。后一种两倍计数会给出约 `3.17 mPa·s`。官方网页写出“每面墙约 2.7”且给出约 `3.1 mPa·s`，数值上与这个两倍计数接近；本案例如实记录此定义差异，不把官方数值当作必须强行匹配的目标。

#### 5.3 与官方教程和官方仓库附带结果的比较

| 比较项目 | 官方参考 | 本次结果 | 判断 |
|---|---|---|---|
| 初始构型 | 1680 原子、总能量 `501411.03 kcal/mol`、压力 `6165751.9 atm` | 三项逐值相同 | 初始体系生成完全一致。 |
| 平衡后期（20–30 ps） | 温度 `300.70 K`、墙间距均值 `27.92 Å` | 温度 `299.53 K`、墙间距均值 `27.75 Å` | 数值接近；短平衡本身仍不足以证明科研稳态。 |
| 剪切最终步数和粒子数 | 200000 步；水 1218 原子、墙 432 原子、离子 30 个 | 相同 | 程序流程和组分一致。 |
| 每面墙平均力大小 | `2.7143 (kcal/mol)/Å` | `2.7027 (kcal/mol)/Å` | 相差约 0.43%，定量上非常接近。 |
| 水速度梯度 | 约 `2.157×10¹⁰ s⁻¹`（仓库参考 `.dat`） | `2.014×10¹⁰ s⁻¹` | 相差约 6.6%；两者都支持网页所述的近线性 Couette 流。 |
| 密度剖面逐点比较 | 不适用 | 不适用 | 仓库附带 `.dat` 为 42 个分箱，但同目录当前 `shearing.lmp` 写的是 0.25 Å 分箱，按本盒子应有 162 个分箱；参考数据与输入文件内部不一致，不能把峰高逐点当作严格金标准。 |
| 剪切温度逐点比较 | 不适用 | 不适用 | 仓库参考剪切日志后期平均温度约 264 K，与输入的 300 K 目标不一致；原因无法仅凭现有文件确定，不应以它否定本次后期约 300 K 的正常控温结果。 |

结论分四层：

1. **程序成功运行：是。** 三阶段均达到预期最终步数，无 ERROR、NaN、Lost atoms，且 `Dangerous builds = 0`。
2. **定性趋势一致：是。** 界面分层、反向墙运动、近线性水速度剖面和方向相反的墙面阻力均符合教程的物理图像。
3. **定量一致：部分成立。** 初始构型逐值一致，平衡后期与墙面平均力接近；速度梯度量级一致。密度峰高和官方黏度不能作严格逐值复现，因为官方参考包存在分箱/输入不一致，且黏度定义存在两倍约定差异。
4. **科研统计收敛：尚未证明。** 平衡仅 30 ps、剪切仅 200 ps、离子数少、仅一个剪切速率、没有独立重复；墙力时间块仍有明显波动。

---

### 6. 当前局限和科研升级方向

本案例已经完成教学级运行、检查、绘图和官方比较；不需要再为“跑通教程”追加模拟。若研究目标从学习转为发表，需要重新设计模拟，而不是直接延长一个输出文件：

1. 平衡时间应显著增加。教程网页估计离子跨越约 1.2 nm 孔道的扩散平衡时间量级约为 1 ns；30 ps 只能作为教学演示。
2. 剪切生产阶段应输出多个独立或近似独立的时间窗口，用于比较密度、速度梯度、墙力、应力和黏度是否稳定。
3. 应使用多个初始速度种子，并报告均值和不确定度；单条轨迹的瞬时标准差不是最终黏度误差棒。
4. 应检查剪切速率、孔宽、体系横向尺寸、离子浓度、壁面相互作用和分箱宽度的敏感性。
5. 若研究离子选择性或电流，应把 Na 与 Cl 分开输出，并计算带符号电荷通量，而不是只统计合并的 ions 组。
6. 发表前应明确应力和黏度的面积/壁面约定，并与独立文献或实验范围比较。

---

---

## 第二部分：代码教学、易错点与理解检查

### 1. `create.lmp`：建立初始体系

这一阶段的任务是建立模拟盒子、上下铝墙、水分子和离子，定义相互作用参数和原子分组，然后把初始状态写入 `create.data`。它不进行真正的动力学演化。

#### 1.1 基本体系设置

```lammps
boundary p p f
units real
atom_style full
bond_style harmonic
angle_style harmonic
pair_style lj/cut/tip4p/long O H O-H H-O-H 0.1546 12.0
kspace_style pppm/tip4p 1.0e-5
kspace_modify slab 3.0
```

- `units real`：采用适合分子体系的 real 单位制。距离通常用 Å，时间用 fs，能量用 kcal/mol，温度用 K，压力用 atm，电荷用基本电荷 (e)。
- `atom_style full`：每个原子除位置和类型外，还可以保存电荷、分子 ID、键和角等信息，适合带电的分子体系。
- `bond_style harmonic`、`angle_style harmonic`：声明键和角使用谐形式。水在后续通过 SHAKE 保持刚性，但 LAMMPS 仍需要知道水分子的键类型、角类型及目标几何。
- `pair_style lj/cut/tip4p/long ...`：使用带长程静电的 TIP4P 类刚性水模型。`O H O-H H-O-H` 指定氧类型、氢类型、O–H 键类型和 H–O–H 角类型；`0.1546` 是 TIP4P 虚拟电荷点相对氧的位置参数，`12.0` Å 是实空间截断距离。
- `kspace_style pppm/tip4p 1.0e-5`：用 PPPM 计算截断距离以外的长程静电作用，目标相对精度为 (10^{-5})。只写短程 `pair_style` 不能完整处理离子和带电水分子的长程库仑作用。
- `kspace_modify slab 3.0`：对 `x`、`y` 周期而 `z` 非周期的薄膜/狭缝几何进行长程静电修正，减小不同 `z` 周期镜像之间的虚假作用。

#### 1.2 晶格、区域和空盒子

```lammps
lattice fcc 4.04
region box block -3 3 -3 3 -5 5
create_box 5 box &
  bond/types 1 &
  angle/types 1 &
  extra/bond/per/atom 2 &
  extra/angle/per/atom 1 &
  extra/special/per/atom 2
```


- `lattice fcc 4.04`：定义晶格常数为 4.04 Å 的面心立方晶格。这里既用它生成铝墙，也让后续默认采用 `units lattice` 的 `region` 尺寸以 4.04 Å 为比例。
- `region box block -3 3 -3 3 -5 5`：只定义一个长方体区域。因为没有写 `units box`，坐标按当前晶格单位解释，所以盒子范围是 `x、y = -12.12 到 12.12 Å`，`z = -20.20 到 20.20 Å`。
- `create_box 5 box`：按照 `box` 区域创建空的模拟盒子，并声明最多有 5 种原子类型。
- `bond/types 1`、`angle/types 1`：为后面的水分子预留 1 种键类型和 1 种角类型。
- `extra/...`：为每个原子预留足够的键、角和特殊邻居存储空间。它们是容量声明，不会自动生成任何键或角。


- `region` 只是“划出一块区域”；
- `create_box` 只是“建立空盒子并声明容量”；
- `create_atoms` 才真正向盒子中放入原子。

如果希望坐标直接使用 Å，可以写：

```lammps
region box block -12.12 12.12 -12.12 12.12 -20.20 20.20 units box
```

#### 1.3 类型标签

```lammps
labelmap atom 1 O 2 H 3 Na+ 4 Cl- 5 WALL
labelmap bond 1 O-H
labelmap angle 1 H-O-H
```

这些命令给数字类型添加可读名称。之后可以写 `type O`、`create_atoms Na+`，而不必到处记忆数字 1、2、3。标签方便阅读，但 LAMMPS 内部仍然使用数字类型。

#### 1.4 创建上下铝墙

```lammps
region rbotwall block -3 3 -3 3 -4 -3
region rtopwall block -3 3 -3 3 3 4
region rwall union 2 rbotwall rtopwall
create_atoms WALL region rwall
```

逐行解释：

- `rbotwall` 和 `rtopwall` 分别是下墙和上墙的区域；坐标仍然采用晶格单位。
- `region rwall union 2 ...` 把两个互不连接的区域合并为一个区域名称。
- `create_atoms WALL region rwall` 按照 FCC 晶格点在这两个区域中创建 WALL 类型的铝原子。

#### 1.5 水分子模板与水的插入

```lammps
region rliquid block INF INF INF INF -2 2
molecule h2omol water.mol
create_atoms 0 region rliquid mol h2omol 482793
```

逐行解释：

- `rliquid` 定义初始液体填充区域；`INF` 表示沿 `x`、`y` 使用整个盒子范围，`z = -2 到 2` 仍按晶格单位解释，即约 `-8.08 到 8.08 Å`。
- `molecule h2omol water.mol` 读取水分子模板，并把模板命名为 `h2omol`。
- `create_atoms 0 ... mol h2omol 482793` 在液体区域的晶格点插入水分子。这里的 `0` 表示原子类型由分子模板本身提供；`482793` 是随机取向使用的随机种子，使该次构型可以重复生成。

`water.mol` 包含：

- 3 个原子、2 根 O–H 键和 1 个 H–O–H 角；
- 原子坐标、O/H 类型和部分电荷；
- SHAKE 所需的约束信息；
- 分子内部特殊邻居关系。

`water.mol` 不是可以不加检查地用于所有代码的“通用水模板”。只有当原子类型、部分电荷、分子几何、键角标签、SHAKE 设置和所用水模型相匹配时才能复用。

#### 1.6 创建离子并设置电荷

```lammps
create_atoms Na+ random 15 5802 rliquid overlap 0.3 maxtry 500
create_atoms Cl- random 15 9012 rliquid overlap 0.3 maxtry 500
set type Na+ charge 1
set type Cl- charge -1
```

- 在 `rliquid` 中随机创建 15 个 Na 和 15 个 Cl；两个不同随机种子使位置可重复。
- `overlap 0.3` 要求新原子与已有原子的距离不能小于 0.3 Å，避免完全重叠；这个距离很小，初始构型仍可能具有较大排斥能，因此后面必须最小化。
- `maxtry 500` 限制每个随机原子的最多尝试次数。
- 两条 `set` 命令分别赋予 (+1e) 和 (-1e) 电荷，使体系加入等量正负离子并保持总电荷中性。


#### 1.7 力场参数和原子分组

```lammps
include parameters.inc
include groups.inc
```

`include` 不是单纯读取文本数据，而是立即执行被包含文件中的 LAMMPS 命令。

`parameters.inc` 主要包含：

- O、H、Na、Cl、Al 的质量；
- Lennard-Jones 相互作用参数；
- O–H 的目标键长和 H–O–H 的目标角度。

其中 `pair_coeff` 的前两个数在 `real` 单位下分别对应 Lennard-Jones 能量尺度和距离尺度。没有单独写出的交叉项通常按 LAMMPS 的混合规则生成；`O WALL` 被显式指定，用来控制水氧与铝墙的作用。

```lammps
bond_coeff O-H 0 0.9572
angle_coeff H-O-H 0 104.52
```

在 harmonic 样式中，第二个数分别是目标 O–H 键长 0.9572 Å 和目标 H–O–H 角 104.52°。前面的力常数为 0，是因为本模型在动力学中由 SHAKE 维持刚性几何；最小化阶段则由 `fix shake ... kbond 2000` 提供约束近似。不能把这两个 0 理解成“水没有键和角”。

`groups.inc` 建立以下分组：

```text
H2O   = O + H
ions  = Na + Cl
fluid = H2O + ions
wall  = WALL
walltop = 位于 z > 0 的墙原子
wallbot = 位于 z < 0 的墙原子
```

`top`、`bot` 区域本身也可能包含流体原子，但最终使用 `intersect wall top` 和 `intersect wall bot`，所以 `walltop`、`wallbot` 只包含墙原子。


#### 1.8 按完整分子删除部分水

```lammps
delete_atoms random fraction 0.15 yes H2O NULL 482793 mol yes
```

这条命令从 `H2O` 组的原子中随机选择约 15%，但 `mol yes` 要求只要选中某个水分子的一个原子，就删除这个水分子的全部原子。因此最终删除的原子比例会明显大于 15%，但不会留下只有 O 或只有 H 的残缺水分子。

本次结果为：

```text
初始水：648 个分子 = 1944 个原子
最终水：406 个分子 = 1218 个原子
墙：432 个原子
离子：30 个原子
总计：1680 个原子
```


#### 1.9 初始输出和 `run 0`

```lammps
dump traj all custom 1 trajectory-create.lammpstrj &
  id mol type element q x y z ix iy iz vx vy vz
dump_modify traj element O H Na Cl Al sort id

thermo 1
thermo_style custom step atoms temp pe etotal press

run 0
write_data create.data nocoeff
```

逐行解释：

- `dump` 定义一条轨迹输出规则；`dump_modify` 为数字类型指定元素名并按原子 ID 排序。
- `thermo 1` 表示每 1 步打印一次热力学信息；由于这里是 `run 0`，只打印第 0 步。
- `run 0` 不进行时间积分，只让 LAMMPS 完成邻居表、能量、力和输出的初始化计算。因此轨迹只有一帧，原子不会运动。
- `write_data create.data nocoeff` 把当前结构写入 data 文件；`nocoeff` 不把完整力场系数写入，所以后续仍要 `include parameters.inc`。

初始势能和压力很高不代表整个教程失败，因为粒子刚被几何方式放入，还没有经过最小化和动力学平衡。它只说明下一阶段确实必要。

---

### 2. `equilibrate.lmp`：消除严重重叠并进行短时间平衡

这一阶段先通过能量最小化降低初始构型中的强排斥和大力，再在 300 K 下运行 30000 步 NVT 分子动力学。它的输出 `equilibrate.data` 是剪切阶段的起点。

输入文件开头再次出现以下全局设置：

```lammps
boundary p p f
units real
atom_style full
bond_style harmonic
angle_style harmonic
pair_style lj/cut/tip4p/long O H O-H H-O-H 0.1546 12.0
kspace_style pppm/tip4p 1.0e-5
kspace_modify slab 3.0
```

这是 **【重复】** 设置，含义与 create 阶段相同。每个新的 LAMMPS 输入进程都要先声明如何解释 data 文件和计算相互作用；`create.data` 并不能代替这些全局模型命令。

#### 2.1 读取前一阶段并恢复模型

```lammps
read_data create.data
include parameters.inc
include groups.inc
```

`read_data` 读取原子、盒子、速度、电荷、键和角等体系状态，但 `create.data` 使用了 `nocoeff`，而且 data 文件不会保存 `fix`、`compute`、`dump`、分组逻辑等运行命令，因此必须重新包含参数文件和分组文件。


#### 2.2 SHAKE 与能量最小化

```lammps
fix myshk H2O shake 1.0e-5 200 0 b O-H a H-O-H kbond 2000

thermo 1
thermo_style custom step temp etotal press

minimize 1.0e-6 1.0e-6 1000 1000
```


- `fix myshk H2O shake ...`：保持水分子的 O–H 键长和 H–O–H 角度接近刚性几何。`1.0e-5` 是约束容差，`200` 是每步最大约束迭代次数，`0` 表示不定期打印 SHAKE 统计。
- `b O-H a H-O-H` 指定约束哪种键和角。
- 能量最小化不是普通时间演化，LAMMPS 在最小化时不能按动力学方式执行 SHAKE，因此 `kbond 2000` 提供用于最小化的强谐约束。这就是日志出现 `Using fix shake with minimization` 的原因。
- `thermo 1` 让最小化的每次迭代都打印信息，便于观察能量和力下降。
- `minimize etol ftol maxiter maxeval` 的四个数字依次是能量容差、力容差、最大迭代次数和最大力评估次数。

本次实际结果：

```text
Stopping criterion = max force evaluations
Iterations, force evaluations = 508, 1001
Energy: 501416.44 → -39112.10 kcal/mol
Force two-norm: 2825414.6 → 117.25
Force max component: 1248444.9 → 43.06
```

它因为达到最大力评估次数而停止，没有达到设置的 (10^{-6}) 力容差，因此不能声称“严格最小化收敛”。但是能量和力已经降低多个数量级，随后 30000 步动力学没有出现 NaN、Lost atoms 或能量爆炸，所以对教学复现来说可以继续。用于科研时则应检查延长最小化、改变算法或放宽/重新定义合理容差后的敏感性。


#### 2.3 最小化轨迹与时间步重置

```lammps
dump trajmin all custom 100 trajectory-equilibrate-minimize.lammpstrj ...
dump_modify trajmin element O H Na Cl Al sort id

minimize ...
undump trajmin
reset_timestep 0
```

`dump trajmin ... 100` 在最小化过程中约每 100 次迭代输出一帧。本次最小化进行了 508 次迭代，实际轨迹有 7 帧：第 0、100、200、300、400、500 次迭代和最小化结束时的第 508 次迭代。OVITO 显示的最后索引为 6，但索引从 0 开始，因此总数是 7。

`undump trajmin` 关闭最小化轨迹，避免后面的 NVT 结果继续写进同一文件。`reset_timestep 0` 把后面的动力学阶段重新从第 0 步计数，使时间轴容易解释；它不会把坐标恢复成最小化前的状态。


#### 2.4 NVT、SHAKE、重心校正和时间步长

```lammps
fix mynvt all nvt temp 300 300 100
fix myshk H2O shake 1.0e-5 200 0 b O-H a H-O-H
fix myrct all recenter NULL NULL 0
timestep 1.0
```

- `fix mynvt all nvt temp 300 300 100`：对全部原子进行 NVT 时间积分，把目标温度从 300 K 保持到 300 K；`100` fs 是恒温阻尼时间，不是运行100步。
- 第二次定义同名 `fix myshk`：从最小化使用的 `kbond` 约束切换为正常分子动力学中的 SHAKE 约束。
- `fix myrct all recenter NULL NULL 0`：不调整 `x`、`y`，把整个体系的整体中心保持在 `z = 0` 附近，防止体系整体沿非周期方向漂移。它不等于固定上下墙之间的距离。
- `timestep 1.0`：每个动力学时间步为 1 fs。


#### 2.5 墙间距变量

```lammps
variable walltopz equal xcm(walltop,z)
variable wallbotz equal xcm(wallbot,z)
variable deltaz equal v_walltopz-v_wallbotz
```

- `xcm(group,z)` 计算指定组在 `z` 方向的质心坐标；
- `walltopz` 和 `wallbotz` 分别是上下墙的 `z` 质心；
- `deltaz` 是两者之差，用来监测墙间距。

#### 2.6 平衡输出和运行时间

```lammps
dump traj all custom 250 trajectory-equilibrate.lammpstrj ...
thermo 250
thermo_style custom step temp etotal press v_deltaz
run 30000
write_data equilibrate.data nocoeff
```

- 运行 30000 步，每步 1 fs，总时间为 30 ps。
- `thermo 250`：从第 0 步到第 30000 步，每 250 步打印一次，共 (30000/250+1=121) 个热力学时刻。
- `dump ... 250`：本次实际轨迹有 121 帧，对应第 0、250、500、…、30000 步。OVITO 显示的最后索引是 120，但从索引 0 到 120 一共有 121 帧。
- `write_data equilibrate.data nocoeff` 保存第 30000 步状态，供剪切阶段读取。

技术上运行成功不等于科研上充分平衡。本次温度总体围绕 300 K 波动，没有 NaN、Lost atoms 或能量爆炸，因此可以继续教学剪切；但墙间距在 30 ps 结束时仍有明显振荡和趋势，不能据此宣称已经达到科研意义上的稳态或收敛。


### 3. `shearing.lmp`：通过反向移动墙产生剪切

这一阶段读取 `equilibrate.data`，让上下墙沿 `x` 方向反向运动，同时统计水、墙和离子沿 `z` 方向的质量密度与平均 `x` 速度。本阶段已经由我亲自运行到第 200000 步，以下代码预期已与实际输出核对。

文件开头的 `boundary`、`units`、`atom_style`、键角样式、TIP4P `pair_style`、PPPM 和 slab 修正均为 **【重复】** 设置，继续沿用前两个阶段的边界、单位和物理模型。剪切阶段改变的是外部驱动、恒温方式和统计方法，而不是重新换一套力场。

#### 3.1 读取平衡状态

```lammps
read_data equilibrate.data
include parameters.inc
include groups.inc
```

剪切不能直接从随机创建的高能构型开始，因此读取短平衡后的结构和速度。`equilibrate.data` 传递体系状态，参数、分组、恒温、剪切和输出命令仍需要在当前输入文件中重新定义。


#### 3.2 剪切流体和墙的恒温

```lammps
compute Tfluid fluid temp/partial 0 1 1
fix mynvt1 fluid nvt temp 300 300 100
fix_modify mynvt1 temp Tfluid

compute Twall wall temp/partial 0 1 1
fix mynvt2 wall nvt temp 300 300 100
fix_modify mynvt2 temp Twall
```

逐行解释：

- `fluid` 包含水和离子，`wall` 包含铝墙；两者分开恒温和监测。
- `temp/partial 0 1 1` 只使用 `y`、`z` 速度分量计算温度，完全排除 `x` 方向。
- 剪切产生的宏观流动位于 `x` 方向。如果把流动速度直接当作热速度参加恒温，恒温器可能削弱我们主动施加的剪切流。
- `fix_modify ... temp ...` 告诉相应 NVT 使用这个部分温度，而不是默认的三维温度。

它不是“把 `x` 方向温度设为 0”，也不是“仍然观察 `x` 方向温度但不控制”。这里的温度计算本身就不包含 `v_x`。更复杂的科研模拟可能需要先扣除局部速度剖面再计算热运动温度，本教程采用排除流动方向的简化办法。

#### 3.3 保持水几何、整体居中和时间步

```lammps
fix myshk H2O shake 1.0e-5 200 0 b O-H a H-O-H
fix myrct all recenter NULL NULL 0
timestep 1.0
```

与平衡阶段相同：SHAKE 保持水的刚性几何，`recenter` 抑制整个体系沿 `z` 漂移，时间步为 1 fs。

#### 3.4 产生剪切的关键命令

```lammps
fix mysf1 walltop setforce 0 NULL NULL
fix mysf2 wallbot setforce 0 NULL NULL
velocity wallbot set -2e-4 NULL NULL
velocity walltop set 2e-4 NULL NULL
```


- `setforce 0 NULL NULL` 把对应墙原子的 `x` 方向力改为 0，而 `NULL` 表示 `y`、`z` 力保持原值。
- `velocity ... set` 把下墙 `v_x` 设为 `-2 × 10⁻⁴ Å/fs`，上墙设为 `+2 × 10⁻⁴ Å/fs`。
- 在 real 单位下，两面墙速率大小约为 20 m/s，相对速度约为 40 m/s。
- 上下墙反向运动，通过墙—流体相互作用拖动中间的水和离子，形成类似平板 Couette 流的剪切速度梯度。
- 分子模拟常使用远高于宏观实验的速度或剪切速率，以便在有限计算时间内得到可测信号；这也是将教程升级为论文研究时必须检查的参数敏感性之一。



#### 3.5 轨迹和热力学输出

```lammps
dump traj all custom 1000 trajectory-shearing.lammpstrj ...
thermo 250
thermo_style custom step temp c_Tfluid c_Twall etotal f_mysf1[1] f_mysf2[1]
thermo_modify temp Tfluid
```

- 轨迹每 1000 步输出一次。从 Step 0 运行到 200000，实际得到 `0, 1000, ..., 200000` 共 201 帧，对应 1 ps 的帧间隔。`200000/1000=200` 是时间区间数，加上第 0 步初始帧才是总帧数。
- thermo 每 250 步打印一次，连同第 0 步实际有 `200000/250+1=801` 个热力学时刻。
- `temp` 在 `thermo_modify temp Tfluid` 后使用流体部分温度，因此它应与 `c_Tfluid` 对应；显式保留两列是为了教学检查。
- `c_Twall` 是墙的部分温度。
- `f_mysf1[1]`、`f_mysf2[1]` 是两个 `setforce` fix 在把 `x` 力清零前记录的墙面总 `x` 力，可用于检查上下墙切向受力，并进一步估算剪切应力。


#### 3.6 沿 `z` 方向划分空间薄片

```lammps
compute cc1 H2O chunk/atom bin/1d z 0.0 0.25
compute cc2 wall chunk/atom bin/1d z 0.0 0.25
compute cc3 ions chunk/atom bin/1d z 0.0 0.25
```

`compute chunk/atom bin/1d` 不直接计算密度。它先按照原子的 `z` 坐标，把原子分配到厚度为 0.25 Å 的一维空间薄片中。

- `cc1`：给水原子分箱；
- `cc2`：给墙原子分箱；
- `cc3`：给 Na 和 Cl 的合并离子组分箱。

`cc1`、`cc2`、`cc3` 只是三个 compute ID，可以换成其他合法名称；真正重要的是它们选择的组、方向、原点和分箱厚度。


#### 3.7 对薄片进行时间平均

```lammps
fix myac1 H2O ave/chunk 10 15000 200000 &
  cc1 density/mass vx file shearing-water.dat
fix myac2 wall ave/chunk 10 15000 200000 &
  cc2 density/mass vx file shearing-wall.dat
fix myac3 ions ave/chunk 10 15000 200000 &
  cc3 density/mass vx file shearing-ions.dat
```

`density/mass` 是一个完整关键词，意思是质量密度，不是“密度除以质量”。`vx` 是每个空间薄片中相应原子的平均 `x` 速度。

三个数字的含义是：

```text
Nevery  = 10      每 10 步采样一次
Nrepeat = 15000   每次平均使用 15000 个样本
Nfreq   = 200000  每到 200000 步生成一组平均结果
```

LAMMPS 先固定输出时刻，再向前选择样本。即输出时刻为第20万步，采用10×15000=15万步，所以第 20万 步之前的 15万 步都是输出的样本。由于之后的run 200000，即一共运行 200000 步，所以最开始的 50000 步不用于最终平均的输出。后 150000 步，每 10 步采样一次。

总运行步数由后面的 `run` 决定。

- 如果 `run 150000`：还没到第 200000 步，不会得到该完整平均结果；
- 如果 `run 200000`：在第 200000 步得到一组结果；
- 如果 `run 400000`：理论上在第 200000 和 400000 步各得到一组结果。

低 `Ncount` 分箱的速度可能没有统计意义

`fix ave/chunk` 输出的 `Ncount` 是多次采样后的平均占据数，所以可以是小数。如果某个薄片在 15000 次采样中只有 3 次出现 1 个原子，它的平均 `Ncount` 就是 `3/15000=0.0002`。此时 `vx` 只由极少数瞬时速度决定，一个偶然的高速原子就能产生非常极端的平均值。

密度图中的低密度区域本身可能有界面意义，不能随意删除；速度图则应根据 `Ncount` 或密度识别统计不足的薄片，并明确报告筛选条件，不能为了得到漂亮直线而隐藏数据。


当前 `ions` 把 Na 和 Cl 合并，因此 `shearing-ions.dat` 不能分别给出阳离子和阴离子分布。如果研究目标是比较 Na 和 Cl 的界面富集，需要为两个离子组分别建立 `compute chunk/atom` 和 `fix ave/chunk`，或从轨迹单独分析。

#### 3.8 运行并保存最终状态

```lammps
run 200000
write_data shearing-final.data nocoeff
```

`run 200000` 才真正规定总运行 200000 步。时间步为 1 fs，因此总剪切时间为 200 ps。`write_data` 保存剪切结束时的结构和速度，便于后续继续计算或复查最终状态。


---
