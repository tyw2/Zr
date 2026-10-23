# Zr Potential Validation

面向锆（Zr）机器学习原子间势的模块化验证工作流。保留原有 **33 项测试**，覆盖 GSFE（11）、NEB（4）、QSD（3）、熔点估算（3）、EOS（2）、结合能（3）、弹性常数（3）和缺陷形成能（4）。

## 项目结构

```text
zr_potential_validation/
├── bin/zrvalidate                统一命令入口（Shell）
├── config/paths.json             所有机器相关路径与 Slurm 账户
├── src/zr_validation/
│   ├── cli.py                    命令行入口
│   ├── settings.py               数值参数与配置加载
│   ├── config_check.py           路径占位符检查
│   ├── core/                     任务生成、资源调度、晶体工具与汇总
│   └── workflows/
│       ├── gsfe.py
│       ├── neb.py
│       ├── qsd.py
│       ├── eos.py
│       ├── cohesive.py
│       ├── elastic/              BCC / FCC / HCP + post
│       ├── formation/            结构 / 计算 / 后处理
│       └── melting/              结构制备 / 温度扫描 / 拟合 / 作业管理
├── tests/                        无 LAMMPS 依赖的生成器冒烟测试
└── pyproject.toml
```

## 配置

打开 `config/paths.json`，将所需路径和 Slurm 账号中的 `***` 替换为真实值。**代码和 README 中未保留以前任何个人绝对路径。** 路径集中在一处管理，数值参数在 `src/zr_validation/settings.py` 中设置。

必填：`potential_root`、`output_root`、`lammps_executable`、`nlh_table`、`slurm_account`。只有使用 GSFE 或 NEB 时才需要填写各自源配置路径。`gsfe_configurations` 对应已有 `input.pos*` 文件夹；`neb_replicas` 对应各 NEB 复制模板文件夹。

`***` 是占位符，不能直接用于真实计算。执行生成前会检查对应必要参数。

## 使用

在项目根目录执行：

```bash
./bin/zrvalidate list
./bin/zrvalidate check --only elastic --phase HCP
./bin/zrvalidate generate --only elastic --phase HCP --seed 0
./bin/zrvalidate generate --only formation --case BO --seed 0
./bin/zrvalidate generate --only eos --seed 0
./bin/zrvalidate generate --only melting --phase BCC --seed 0
./bin/zrvalidate generate --only all --seed 0
```

`generate` **仅生成** LAMMPS / Slurm 输入文件，并不自动执行长时间计算。提交短作业时进入所设 `output_root`：

```bash
./submit_elastic.sh
./submit_formation.sh
./submit_eos.sh
./submit_all.sh
```

熔点扫描采用独立节点作业和依赖管理，示例：

```bash
python3 submit_melting.py --seed 0 --phase BCC
python3 submit_melting.py --seed 0 --phase BCC --execute
```

分析文件会随已选择的任务一起写入输出目录。静态物性生成 `collect_static.py`，用于合并多势能的结果。

## 依赖与验证

输入文件生成需要 Python 3.9+。温度扫描的 RDF 后处理还需要 `numpy`、`scipy`、`matplotlib`；若要在集群计算，应使用含 MLIP3 与 NLH 表势的 LAMMPS 和 Slurm，并按对应模块设置 MPI。

```bash
python3 -m compileall -q src
python3 -m unittest discover -s tests -v
```

这些测试只检查文件生成和基本路径逻辑；不代替真实 LAMMPS 势能计算的收敛性或物理准确性验证。计算前建议先用一个 seed 验证生成结构及 LAMMPS 命令。

## 兼容性说明

保留原有输出文件夹及输出文件名（如 `seed_0/Elastic-HCP`、`elastic_result.csv`），但不再提供平铺的 `generate_*.py` 入口。已生成的 Slurm 与 LAMMPS 文件仍可按原方式执行。迁移已有输出时务必先备份；切换 `output_root` 可以避免任务跳过逻辑误判。
