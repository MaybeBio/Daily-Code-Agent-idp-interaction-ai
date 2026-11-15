# 代码快照清单

- 上游仓库: https://github.com/c2fd/LLPS-GPU-Simulator
- 快照日期: 2026-11-15
- 快照 commit: [`dd7c4d4d9fd1`](https://github.com/c2fd/LLPS-GPU-Simulator/commit/dd7c4d4d9fd16a359327e88765be9137df662eb9)
- 归档策略: `repo_copy=code` — **代码摘录,不是完整仓库**。保留源码、构建/配置文件、README、文档;数据、二进制、图片、权重不保留。

## 保留 (186 个文件, 1.1 MB)

| 类别 | 文件数 | 大小 |
|---|---:|---:|
| 源码 | 183 | 1.1 MB |
| 构建/配置 | 1 | 1.9 KB |
| README | 1 | 565 B |
| 文档 | 1 | 11 B |

按后缀: `.h` 165 (1.0 MB)、`.cu` 10 (31.7 KB)、`.cpp` 6 (10.6 KB)、`Makefile` 1 (1.9 KB)、`.sh` 2 (1.1 KB)、`.md` 2 (576 B)

> README / 文档中的相对图片与链接已改写为 GitHub 绝对地址。

## 省略 (75 个文件, 416.5 KB)

| 原因 | 文件数 | 大小 |
|---|---:|---:|
| 后缀不在白名单 | 75 | 416.5 KB |

体积最大的省略文件:

| 路径 | 大小 |
|---|---:|
| [`cusp/detail/profiler.inl`](https://github.com/c2fd/LLPS-GPU-Simulator/raw/HEAD/cusp/detail/profiler.inl) | 33.9 KB |
| [`cusp/opengl/spy/detail/matrix_canvas.inl`](https://github.com/c2fd/LLPS-GPU-Simulator/raw/HEAD/cusp/opengl/spy/detail/matrix_canvas.inl) | 27.0 KB |
| [`cusp/krylov/detail/bicgstab_m.inl`](https://github.com/c2fd/LLPS-GPU-Simulator/raw/HEAD/cusp/krylov/detail/bicgstab_m.inl) | 23.7 KB |
| [`cusp/detail/blas.inl`](https://github.com/c2fd/LLPS-GPU-Simulator/raw/HEAD/cusp/detail/blas.inl) | 21.7 KB |
| [`cusp/precond/detail/ainv.inl`](https://github.com/c2fd/LLPS-GPU-Simulator/raw/HEAD/cusp/precond/detail/ainv.inl) | 21.2 KB |
| [`cusp/krylov/detail/cg_m.inl`](https://github.com/c2fd/LLPS-GPU-Simulator/raw/HEAD/cusp/krylov/detail/cg_m.inl) | 15.6 KB |
| [`cusp/io/detail/matrix_market.inl`](https://github.com/c2fd/LLPS-GPU-Simulator/raw/HEAD/cusp/io/detail/matrix_market.inl) | 15.1 KB |
| [`cusp/detail/verify.inl`](https://github.com/c2fd/LLPS-GPU-Simulator/raw/HEAD/cusp/detail/verify.inl) | 13.4 KB |
| [`cusp/opengl/spy/detail/gl_util.inl`](https://github.com/c2fd/LLPS-GPU-Simulator/raw/HEAD/cusp/opengl/spy/detail/gl_util.inl) | 10.2 KB |
| [`cusp/detail/format_utils.inl`](https://github.com/c2fd/LLPS-GPU-Simulator/raw/HEAD/cusp/detail/format_utils.inl) | 9.9 KB |
| [`cusp/precond/aggregation/detail/aggregate.inl`](https://github.com/c2fd/LLPS-GPU-Simulator/raw/HEAD/cusp/precond/aggregation/detail/aggregate.inl) | 9.7 KB |
| [`cusp/.DS_Store`](https://github.com/c2fd/LLPS-GPU-Simulator/raw/HEAD/cusp/.DS_Store) | 8.0 KB |
| [`cusp/opengl/spy/detail/matrix_canvas_svg_output.inl`](https://github.com/c2fd/LLPS-GPU-Simulator/raw/HEAD/cusp/opengl/spy/detail/matrix_canvas_svg_output.inl) | 7.4 KB |
| [`cusp/precond/aggregation/detail/smooth.inl`](https://github.com/c2fd/LLPS-GPU-Simulator/raw/HEAD/cusp/precond/aggregation/detail/smooth.inl) | 7.1 KB |
| [`cusp/precond/aggregation/detail/strength.inl`](https://github.com/c2fd/LLPS-GPU-Simulator/raw/HEAD/cusp/precond/aggregation/detail/strength.inl) | 7.0 KB |
| [`cusp/gallery/stencil.inl`](https://github.com/c2fd/LLPS-GPU-Simulator/raw/HEAD/cusp/gallery/stencil.inl) | 6.5 KB |
| [`cusp/relaxation/detail/polynomial.inl`](https://github.com/c2fd/LLPS-GPU-Simulator/raw/HEAD/cusp/relaxation/detail/polynomial.inl) | 6.3 KB |
| [`cusp/precond/aggregation/detail/smoothed_aggregation.inl`](https://github.com/c2fd/LLPS-GPU-Simulator/raw/HEAD/cusp/precond/aggregation/detail/smoothed_aggregation.inl) | 6.3 KB |
| [`.DS_Store`](https://github.com/c2fd/LLPS-GPU-Simulator/raw/HEAD/.DS_Store) | 6.0 KB |
| [`cusp/detail/.DS_Store`](https://github.com/c2fd/LLPS-GPU-Simulator/raw/HEAD/cusp/detail/.DS_Store) | 6.0 KB |

> 省略规则:普通源码 >2 MB、文档/构建配置 >1 MB(README 与 .ipynb 不设限),或后缀不在白名单(图片/权重/数据表/二进制/通用 .txt/.log/.csv/.json/.xml 等)、以及子模块。需要全文时按上表链接直接取回,或用 git 克隆上游仓库。

> 本清单由 `scripts/archive.py` 每次归档时自动生成并覆盖。
