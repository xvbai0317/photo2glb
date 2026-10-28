---
name: photo-to-glb
description: 从单张产品照片重建程序化网格并导出米制 GLB，自带纯 numpy 的 glb 写入器、导出前自检与真实渲染验证流程。当用户要求把图片或照片做成 GLB、glTF、3D 模型文件时使用。不用于绑骨、动画、或下载现成模型。
---

# Photo to GLB

把一张零件照片变成可用的 GLB：**先测量、再分层建模、最后用真实渲染器验证**。
不要一上来就写网格；照片里量不到的尺寸靠推断，必须在交付时说明。

配套资源（按需读取，不要一次全读）：

- `scripts/glb_writer.py` —— 可直接 import 的写入器 + 通用几何生成器，**不要重写**
- `assets/verify_viewer.html` —— three.js 验证页面（带影棚 IBL），截图用
- `references/geometry-patterns.md` —— 车削/六角/滚花/螺纹/沉孔/带孔平面 的建模套路

## 1. 测量（动手建模之前）

用 Pillow 裁切放大局部，逐块量像素：

1. 找一个已知尺寸做标尺（铭牌、标准螺纹、把手直径），得到 `px_per_mm`。
2. 量主轮廓关键点：总宽、总高、厚、倒角起点、孔中心与直径、各层高度。
3. 写出"毫米尺寸表"，后面所有常量都来自它，不要现编。
4. 单独记录**特征清单**（决定像不像）：颜色环位置、倒角数量、接口位置、贴纸位置、把手形状。

照片只有一个视角：背面、内腔、流道一律按对称/常规结构推断，并在交付说明里写清。

## 2. 分层

- 单位毫米、Y 向上、+Z 朝观察者；导出时根节点乘 `diag(0.001)` 得到米制资产。
- 节点树按装配关系命名（本体 / 铭牌 / 各功能件），每个节点挂若干 `(材质序号, Geo)`。
- 材质先定表（不锈钢/机加工面/内腔/色环/贴纸…），后面每个零件只选材质，不再调参。

## 3. 建模（复用生成器，别手搓三角形）

在 `references/geometry-patterns.md` 里挑现成套路：

| 需求 | 生成器 |
|---|---|
| 回转体、阀体堆叠、阀杆 | `g_revolve`（轮廓 + 平滑/硬边法线） |
| 六角螺母、六角底座 | `g_ngon`（对边尺寸） |
| 滚花色环 | `g_fluted` |
| 螺纹颈 | `g_threaded` |
| 内螺纹接口（沉孔） | `g_bore`（内法线 + 封底） |
| 平面挖一个圆孔 | `emit_face_hole`（放射框架 + 射线求交） |
| 平板挖多个圆孔（手轮） | `strip_plate`（沿长轴条带 + 圆弧取样，孔洞边缘顶点精确落在圆上） |
| 倒角块体 | 内缩多边形 `inset_poly` + 逐边倒角带 + 长侧面 |

## 4. 五条硬约束（每一条都真实踩过坑）

1. **绕序必须与法线同向**。`Geo.flat_quad` / `flat_tri` / `smooth_quad` 都要用"几何法线 vs 顶点法线"判方向再翻绕序；只按 `outward` 参数判会在锥面/内腔上判反 —— 结果是正面被剔除，孔看起来像"透空"。
2. **内腔（孔、沉孔）**要做成封闭体：法线朝轴心、绕序跟随法线、端面封底。缺一个就只能看见背景色。
3. **材质键名必须是 `pbrMetallicRoughness`**（在 material 顶层）。写成 `pbr` 之类会静默退化成默认白色金属，全模型发白。
4. **POSITION accessor 必须给 `min`/`max`**；bufferView 的 `target` 用 34962/34963；JSON 块与 BIN 块都要 4 字节对齐（JSON 补空格，BIN 补 0）。
5. **baseColorFactor 是线性值**：把 sRGB 取 `pow 2.2` 再写。贴图 PNG 保持 sRGB，UV 的 v=0 在图像顶部。

## 5. 写文件

```bash
python -c "import sys; sys.path.insert(0,r'<skill>/scripts'); import glb_writer"
```

```python
import glb_writer as G
root = G.Node("Part"); body = G.Node("Body"); root.children.append(body)
G.build_something(body)                      # 在毫米空间造几何
G.write_glb("out.glb", root, mats, texture_png=label_png,
            root_matrix=G.MM_TO_M)           # 根节点 0.001 缩放
```

写完立刻打印尺寸/三角数/节点数，任何一项为 0 都是错的。
`scripts/glb_writer.py --self-test <dir>` 会生成一个最小 GLB 并自校验头部。

## 6. 验证（必须用真实渲染器）

自写的软光栅只能当参考，**不要用它下结论**（会误报"缺面""透空"）。

```bash
# 1) 本地起服务
python -m http.server 8777
# 2) 无头 Chrome 截图（输出 1000x1000）
chrome --headless=new --use-angle=swiftshader --hide-scrollbars \
       --window-size=1000,1000 --virtual-time-budget=20000 \
       --screenshot=shot.png "http://127.0.0.1:8777/assets/verify_viewer.html?az=38&el=18&dist=0.30"
```

看这几项，逐项和照片比对，任何一项明显不符就改代码重出图：

- 金属是否发亮（只给 IBL，金属在暗环境里会发黑 → 用页面自带的影棚环境，必要时调整曝光）
- 颜色环、贴纸文字方向与位置
- 接口是不是"暗孔"而不是"白盘/透空"（白盘通常是内腔法线反了或材质高光过曝）
- 把手轮廓是否与照片同形（长圆板 vs 圆盘、孔数、齿边）
- 三个角度：`az=0` 正视、`az=±38` 斜视、`el` 20 左右俯视

## 7. 交付

给出 GLB 链接、包围盒（米）、三角面数、节点结构，并说明：
隐藏面为推断、模型为程序化重建（非摄影测量）、靠 IBL 照明所以不同查看器观感会有差异。

## 完成检查表

- [ ] 尺寸表来自照片测量，常量集中在一处
- [ ] 每个零件都是封闭体，没有只剩单面的地方
- [ ] 孔洞从每个角度都是暗腔，不是背景色
- [ ] 材质键名与线性色正确，贴纸方向正确
- [ ] 根节点 0.001 缩放，包围盒是真实毫米尺寸/1000
- [ ] 至少 3 个角度截图与照片逐项比对通过
- [ ] 交付说明含推断项与限制
