# 代码快照清单

- 上游仓库: https://github.com/JuliaMolSim/Molly.jl
- 快照日期: 未知
- 归档策略: `repo_copy=code` — **代码摘录,不是完整仓库**。保留源码、构建/配置文件、README、文档;数据、二进制、图片、权重不保留。

## 保留 (134 个文件, 10.7 MB)

| 类别 | 文件数 | 大小 |
|---|---:|---:|
| 源码 | 123 | 10.6 MB |
| 文档 | 4 | 50.6 KB |
| 构建/配置 | 6 | 9.3 KB |
| README | 1 | 8.6 KB |

按后缀: `.ipynb` 2 (7.6 MB)、`.jl` 114 (3.0 MB)、`.md` 4 (58.9 KB)、`.py` 7 (46.2 KB)、`.toml` 3 (6.3 KB)、`.yml` 3 (3.0 KB)、`.bib` 1 (310 B)

> README / 文档中的相对图片与链接已改写为 GitHub 绝对地址。

## 省略 (177 个文件, 46.3 MB)

| 原因 | 文件数 | 大小 |
|---|---:|---:|
| 后缀不在白名单 | 177 | 46.3 MB |

体积最大的省略文件:

| 路径 | 大小 |
|---|---:|
| [`data/force_fields/charmm36.xml`](https://github.com/JuliaMolSim/Molly.jl/raw/HEAD/data/force_fields/charmm36.xml) | 5.6 MB |
| [`data/tyk2_ejm50.pdb`](https://github.com/JuliaMolSim/Molly.jl/raw/HEAD/data/tyk2_ejm50.pdb) | 3.8 MB |
| [`data/tyk2_ejm31.pdb`](https://github.com/JuliaMolSim/Molly.jl/raw/HEAD/data/tyk2_ejm31.pdb) | 3.8 MB |
| [`data/force_fields/amber14/GLYCAM_06j-1.xml`](https://github.com/JuliaMolSim/Molly.jl/raw/HEAD/data/force_fields/amber14/GLYCAM_06j-1.xml) | 3.2 MB |
| [`data/6mrr_equil.pdb`](https://github.com/JuliaMolSim/Molly.jl/raw/HEAD/data/6mrr_equil.pdb) | 1.2 MB |
| [`data/1ssu.pdb`](https://github.com/JuliaMolSim/Molly.jl/raw/HEAD/data/1ssu.pdb) | 1.2 MB |
| [`data/openmm_6mrr/charmm/velocities_100steps.txt`](https://github.com/JuliaMolSim/Molly.jl/raw/HEAD/data/openmm_6mrr/charmm/velocities_100steps.txt) | 922.7 KB |
| [`data/openmm_6mrr/amber/velocities_100steps.txt`](https://github.com/JuliaMolSim/Molly.jl/raw/HEAD/data/openmm_6mrr/amber/velocities_100steps.txt) | 918.2 KB |
| [`data/openmm_6mrr/velocities_300K.txt`](https://github.com/JuliaMolSim/Molly.jl/raw/HEAD/data/openmm_6mrr/velocities_300K.txt) | 917.9 KB |
| [`data/openmm_6mrr/amber/forces_bond_only.txt`](https://github.com/JuliaMolSim/Molly.jl/raw/HEAD/data/openmm_6mrr/amber/forces_bond_only.txt) | 885.7 KB |
| [`data/openmm_6mrr/amber/forces_angle_only.txt`](https://github.com/JuliaMolSim/Molly.jl/raw/HEAD/data/openmm_6mrr/amber/forces_angle_only.txt) | 884.5 KB |
| [`data/openmm_6mrr/amber/forces_all_cut.txt`](https://github.com/JuliaMolSim/Molly.jl/raw/HEAD/data/openmm_6mrr/amber/forces_all_cut.txt) | 884.4 KB |
| [`data/openmm_6mrr/amber/forces_all_pme_exact.txt`](https://github.com/JuliaMolSim/Molly.jl/raw/HEAD/data/openmm_6mrr/amber/forces_all_pme_exact.txt) | 884.4 KB |
| [`data/openmm_6mrr/amber/forces_all_pme.txt`](https://github.com/JuliaMolSim/Molly.jl/raw/HEAD/data/openmm_6mrr/amber/forces_all_pme.txt) | 884.4 KB |
| [`data/openmm_6mrr/amber/forces_coul_only.txt`](https://github.com/JuliaMolSim/Molly.jl/raw/HEAD/data/openmm_6mrr/amber/forces_coul_only.txt) | 882.2 KB |
| [`data/openmm_6mrr/charmm/forces.txt`](https://github.com/JuliaMolSim/Molly.jl/raw/HEAD/data/openmm_6mrr/charmm/forces.txt) | 882.1 KB |
| [`data/openmm_6mrr/charmm/coordinates_100steps.txt`](https://github.com/JuliaMolSim/Molly.jl/raw/HEAD/data/openmm_6mrr/charmm/coordinates_100steps.txt) | 866.6 KB |
| [`data/openmm_6mrr/amber/coordinates_100steps.txt`](https://github.com/JuliaMolSim/Molly.jl/raw/HEAD/data/openmm_6mrr/amber/coordinates_100steps.txt) | 866.5 KB |
| [`data/force_fields/gaff.xml`](https://github.com/JuliaMolSim/Molly.jl/raw/HEAD/data/force_fields/gaff.xml) | 634.7 KB |
| [`data/ethanol_garnet.pdb`](https://github.com/JuliaMolSim/Molly.jl/raw/HEAD/data/ethanol_garnet.pdb) | 552.2 KB |

> 省略规则:普通源码 >2 MB、文档/构建配置 >1 MB(README 与 .ipynb 不设限),或后缀不在白名单(图片/权重/数据表/二进制/通用 .txt/.log/.csv/.json/.xml 等)、以及子模块。需要全文时按上表链接直接取回,或用 git 克隆上游仓库。

> 本清单由 `scripts/archive.py` 每次归档时自动生成并覆盖。
