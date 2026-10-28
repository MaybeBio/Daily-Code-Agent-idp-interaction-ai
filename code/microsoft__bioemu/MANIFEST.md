# 代码快照清单

- 上游仓库: https://github.com/microsoft/bioemu
- 快照日期: 2026-10-28
- 快照 commit: [`28d6bac7bb10`](https://github.com/microsoft/bioemu/commit/28d6bac7bb102ee2948a816da2ca280aa2c6ce6b)
- 归档策略: `repo_copy=code` — **代码摘录,不是完整仓库**。保留源码、构建/配置文件、README、文档;数据、二进制、图片、权重不保留。

## 保留 (128 个文件, 1.7 MB)

| 类别 | 文件数 | 大小 |
|---|---:|---:|
| 源码 | 98 | 858.3 KB |
| 构建/配置 | 15 | 761.9 KB |
| 文档 | 7 | 52.4 KB |
| README | 8 | 20.1 KB |

按后缀: `.py` 96 (841.9 KB)、`.lock` 1 (743.5 KB)、`.md` 13 (60.3 KB)、`.ipynb` 1 (15.5 KB)、`无后缀` 2 (12.2 KB)、`.yml` 4 (8.9 KB)、`.yaml` 9 (5.5 KB)、`pyproject.toml` 1 (4.0 KB)、`.sh` 1 (918 B)

> README / 文档中的相对图片与链接已改写为 GitHub 绝对地址。

## 省略 (22 个文件, 35.1 MB)

| 原因 | 文件数 | 大小 |
|---|---:|---:|
| 后缀不在白名单 | 22 | 35.1 MB |

体积最大的省略文件:

| 路径 | 大小 |
|---|---:|
| [`assets/bioemu_paper.pdf`](https://github.com/microsoft/bioemu/raw/HEAD/assets/bioemu_paper.pdf) | 29.8 MB |
| [`dev/afdb_cluster_definitions/afdb_cluster_definitions.json.gz`](https://github.com/microsoft/bioemu/raw/HEAD/dev/afdb_cluster_definitions/afdb_cluster_definitions.json.gz) | 4.5 MB |
| [`assets/emu.png`](https://github.com/microsoft/bioemu/raw/HEAD/assets/emu.png) | 446.4 KB |
| [`tests/expected.npz`](https://github.com/microsoft/bioemu/raw/HEAD/tests/expected.npz) | 136.2 KB |
| [`tests/test_colabfold_inline/golden/trpcage_ref_pair.npy`](https://github.com/microsoft/bioemu/raw/HEAD/tests/test_colabfold_inline/golden/trpcage_ref_pair.npy) | 100.1 KB |
| [`notebooks/PPFT/HHH_rd1_0335.pdb`](https://github.com/microsoft/bioemu/raw/HEAD/notebooks/PPFT/HHH_rd1_0335.pdb) | 61.6 KB |
| [`tests/state_dict.ptkeep`](https://github.com/microsoft/bioemu/raw/HEAD/tests/state_dict.ptkeep) | 27.3 KB |
| [`tests/test_colabfold_inline/golden/chignolin_ref_pair.npy`](https://github.com/microsoft/bioemu/raw/HEAD/tests/test_colabfold_inline/golden/chignolin_ref_pair.npy) | 25.1 KB |
| [`tests/test_colabfold_inline/golden/trpcage_ref_single.npy`](https://github.com/microsoft/bioemu/raw/HEAD/tests/test_colabfold_inline/golden/trpcage_ref_single.npy) | 15.1 KB |
| [`src/_vendor/alphafold/common/stereo_chemical_props.txt`](https://github.com/microsoft/bioemu/raw/HEAD/src/_vendor/alphafold/common/stereo_chemical_props.txt) | 8.9 KB |
| [`tests/test_colabfold_inline/golden/chignolin_ref_single.npy`](https://github.com/microsoft/bioemu/raw/HEAD/tests/test_colabfold_inline/golden/chignolin_ref_single.npy) | 7.6 KB |
| [`tests/test_data/cln_bad_sample_hpacked.pdb`](https://github.com/microsoft/bioemu/raw/HEAD/tests/test_data/cln_bad_sample_hpacked.pdb) | 7.3 KB |
| [`tests/expected.pdb`](https://github.com/microsoft/bioemu/raw/HEAD/tests/expected.pdb) | 3.9 KB |
| [`tests/test_data/cln_bad_sample.pdb`](https://github.com/microsoft/bioemu/raw/HEAD/tests/test_data/cln_bad_sample.pdb) | 3.8 KB |
| [`tests/training/chignolin.pdb`](https://github.com/microsoft/bioemu/raw/HEAD/tests/training/chignolin.pdb) | 3.8 KB |
| [`.gitignore`](https://github.com/microsoft/bioemu/raw/HEAD/.gitignore) | 1.6 KB |
| [`src/bioemu/colabfold_inline/LICENSES/COLABFOLD-MIT`](https://github.com/microsoft/bioemu/raw/HEAD/src/bioemu/colabfold_inline/LICENSES/COLABFOLD-MIT) | 1.0 KB |
| [`cgmanifest.json`](https://github.com/microsoft/bioemu/raw/HEAD/cgmanifest.json) | 976 B |
| [`MANIFEST.in`](https://github.com/microsoft/bioemu/raw/HEAD/MANIFEST.in) | 207 B |
| [`tests/test_colabfold_inline/golden/trpcage.a3m`](https://github.com/microsoft/bioemu/raw/HEAD/tests/test_colabfold_inline/golden/trpcage.a3m) | 192 B |

> 省略规则:普通源码 >2 MB、文档/构建配置 >1 MB(README 与 .ipynb 不设限),或后缀不在白名单(图片/权重/数据表/二进制/通用 .txt/.log/.csv/.json/.xml 等)、以及子模块。需要全文时按上表链接直接取回,或用 git 克隆上游仓库。

> 本清单由 `scripts/archive.py` 每次归档时自动生成并覆盖。
