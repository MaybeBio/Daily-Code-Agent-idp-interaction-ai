# 代码快照清单

- 上游仓库: https://github.com/KULL-Centre/CALVADOS
- 快照日期: 2026-10-28
- 快照 commit: [`866321930af4`](https://github.com/KULL-Centre/CALVADOS/commit/866321930af416dedaeee2e3ec20ebcb5bd3b473)
- 归档策略: `repo_copy=code` — **代码摘录,不是完整仓库**。保留源码、构建/配置文件、README、文档;数据、二进制、图片、权重不保留。

## 保留 (81 个文件, 336.2 KB)

| 类别 | 文件数 | 大小 |
|---|---:|---:|
| 源码 | 41 | 283.0 KB |
| 文档 | 11 | 35.9 KB |
| README | 15 | 11.9 KB |
| 构建/配置 | 14 | 5.4 KB |

按后缀: `.py` 39 (277.1 KB)、`无后缀` 1 (34.3 KB)、`.md` 15 (11.9 KB)、`.ipynb` 2 (5.9 KB)、`.yml` 3 (1.7 KB)、`.rst` 10 (1.5 KB)、`pyproject.toml` 1 (1.3 KB)、`.yaml` 7 (1.3 KB)、`Makefile` 1 (634 B)、`.ini` 1 (265 B)、`.txt` 1 (130 B)

> README / 文档中的相对图片与链接已改写为 GitHub 绝对地址。

## 省略 (50 个文件, 2.5 MB)

| 原因 | 文件数 | 大小 |
|---|---:|---:|
| 后缀不在白名单 | 50 | 2.5 MB |

体积最大的省略文件:

| 路径 | 大小 |
|---|---:|
| [`examples/custom_restraints/input/hnRNPA1S.pdb`](https://github.com/KULL-Centre/CALVADOS/raw/HEAD/examples/custom_restraints/input/hnRNPA1S.pdb) | 361.7 KB |
| [`examples/single_MDP/input/hnRNPA1S.pdb`](https://github.com/KULL-Centre/CALVADOS/raw/HEAD/examples/single_MDP/input/hnRNPA1S.pdb) | 352.7 KB |
| [`examples/slab_MDP/input/hnRNPA1S.pdb`](https://github.com/KULL-Centre/CALVADOS/raw/HEAD/examples/slab_MDP/input/hnRNPA1S.pdb) | 352.7 KB |
| [`examples/single_MDP/input/TIA1.pdb`](https://github.com/KULL-Centre/CALVADOS/raw/HEAD/examples/single_MDP/input/TIA1.pdb) | 325.7 KB |
| [`examples/slab_IDR_MDP/input/TIA1.pdb`](https://github.com/KULL-Centre/CALVADOS/raw/HEAD/examples/slab_IDR_MDP/input/TIA1.pdb) | 325.7 KB |
| [`examples/slab_MDP/input/TIA1.pdb`](https://github.com/KULL-Centre/CALVADOS/raw/HEAD/examples/slab_MDP/input/TIA1.pdb) | 325.7 KB |
| [`examples/two_IDR_MDP/input/TIA1.pdb`](https://github.com/KULL-Centre/CALVADOS/raw/HEAD/examples/two_IDR_MDP/input/TIA1.pdb) | 325.7 KB |
| [`examples/single_AF_CALVADOS/input/AF-P00698.pdb`](https://github.com/KULL-Centre/CALVADOS/raw/HEAD/examples/single_AF_CALVADOS/input/AF-P00698.pdb) | 94.5 KB |
| [`examples/single_dsRNA/input/dspolyR12.pdb`](https://github.com/KULL-Centre/CALVADOS/raw/HEAD/examples/single_dsRNA/input/dspolyR12.pdb) | 50.4 KB |
| [`examples/single_AF_CALVADOS/input/AF-P00698.json`](https://github.com/KULL-Centre/CALVADOS/raw/HEAD/examples/single_AF_CALVADOS/input/AF-P00698.json) | 46.5 KB |
| [`.gitignore`](https://github.com/KULL-Centre/CALVADOS/raw/HEAD/.gitignore) | 4.5 KB |
| [`tests/data/residues_C2RNA.csv`](https://github.com/KULL-Centre/CALVADOS/raw/HEAD/tests/data/residues_C2RNA.csv) | 1.2 KB |
| [`examples/single_dsRNA/input/residues_C2RNA.csv`](https://github.com/KULL-Centre/CALVADOS/raw/HEAD/examples/single_dsRNA/input/residues_C2RNA.csv) | 1.1 KB |
| [`examples/single_RNA/residues_C2RNA.csv`](https://github.com/KULL-Centre/CALVADOS/raw/HEAD/examples/single_RNA/residues_C2RNA.csv) | 1.1 KB |
| [`examples/slab_mixed/input/residues_C2RNA.csv`](https://github.com/KULL-Centre/CALVADOS/raw/HEAD/examples/slab_mixed/input/residues_C2RNA.csv) | 1.1 KB |
| [`calvados/data/residues.csv`](https://github.com/KULL-Centre/CALVADOS/raw/HEAD/calvados/data/residues.csv) | 1.0 KB |
| [`examples/single_pIDR/input/residues_pCALVADOS2.csv`](https://github.com/KULL-Centre/CALVADOS/raw/HEAD/examples/single_pIDR/input/residues_pCALVADOS2.csv) | 996 B |
| [`examples/custom_restraints/input/residues_CALVADOS3.csv`](https://github.com/KULL-Centre/CALVADOS/raw/HEAD/examples/custom_restraints/input/residues_CALVADOS3.csv) | 931 B |
| [`examples/single_AF_CALVADOS/input/residues.csv`](https://github.com/KULL-Centre/CALVADOS/raw/HEAD/examples/single_AF_CALVADOS/input/residues.csv) | 931 B |
| [`examples/single_MDP/input/residues_CALVADOS3.csv`](https://github.com/KULL-Centre/CALVADOS/raw/HEAD/examples/single_MDP/input/residues_CALVADOS3.csv) | 931 B |

> 省略规则:普通源码 >2 MB、文档/构建配置 >1 MB(README 与 .ipynb 不设限),或后缀不在白名单(图片/权重/数据表/二进制/通用 .txt/.log/.csv/.json/.xml 等)、以及子模块。需要全文时按上表链接直接取回,或用 git 克隆上游仓库。

> 本清单由 `scripts/archive.py` 每次归档时自动生成并覆盖。
