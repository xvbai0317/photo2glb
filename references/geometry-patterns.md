# 几何套路与坑（photo-to-glb）

生成器全部在 `scripts/glb_writer.py` 里，直接用，不要重写。下面只讲**怎么用**和**哪里会翻车**。
所有坐标都是毫米、Y 向上、+Z 朝观察者。

---

## 0. 通用：Geo 与绕序

```python
g = Geo()
g.flat_quad(pts4, outward, uv=None)      # 硬边四边形，自动按 outward 翻绕序
g.flat_tri(pts3, outward, uv=None)
g.smooth_quad(pts4, nrms4, outward=None) # 逐顶点法线，绕序按法线和判定
g.extend(other, matrix=None)             # 把另一个 Geo 变换后并入（法线用逆转置）
g.check()                                # 返回绕序/法线不一致的三角形数，必须为 0
```

**坑 1：`outward` 只是提示。** `flat_quad` 用几何叉积和 `outward` 点乘决定是否翻面；在锥面、内腔、倾斜侧壁上 `outward` 很容易给反。
**坑 2：`smooth_quad` 不要传错法线方向。** 它用「几何法线 vs 四个顶点法线之和」判绕序，只要顶点法线是朝外的，绕序就自动正确。这是修掉「孔看起来是白盘/透空」的关键。

内腔一律这样构造：法线朝轴心 → 绕序自动跟随 → 端面补封底三角扇。缺封底就会看到背景色。

---

## 1. `g_revolve(profile, seg=64, smooth=True, cap_top=False, cap_bottom=False, flip=False)`

按 `(r, y)` 轮廓车削。`profile` 是从下往上的折线；每个相邻点对生成一圈四边形。

- 平滑：`smooth=True`，法线由相邻段法线求平均，适合圆柱、圆角、球头。
- 硬边：`smooth=False`，每段一个法线，适合螺纹、阶梯、滚花。
- 封端：`cap_top / cap_bottom`，用三角扇封住首末圈（`r=0` 的轮廓端点不要再封）。
- 翻转：`flip=True` 反向所有法线并反转绕序 —— 做**内腔**时用。
- 轮廓要求：`r` 单调不要写成 0 宽段；相邻点不要重合，否则法线求平均会得到 0 向量。

**坑 3：平滑分支的四边形点序必须是环形序** `(i,j) → (i,j+1) → (i+1,j+1) → (i+1,j)`。写成 `(i,j) → (i,j+1) → (i+1,j) → (i+1,j+1)` 会让四边形自相交，两个三角形里有一个朝内 —— 表现为零件表面出现一片片黑洞。用 `Geo.check()` 一次就能查出来。

**坑 4：`g_revolve` 的 `y` 是局部轴**，需要朝别的方向就用 `frame(origin, ydir)` 变换：

```python
g.extend(g_cyl(2.15, 0, 5.4, seg=24), frame((0, 28.8, 0), (0, 1, 0)))
```

---

## 2. `g_cyl(r, y0, y1, seg=64, chamfer_top=0, chamfer_bot=0, cap=True)`

最常用。`g_cyl` = `g_revolve` 的两点轮廓 + 可选倒角环。倒角值不要超过 `(y1-y0)/2` 和 `r`，否则倒角环自交。

---

## 3. `g_ngon(n, af, y0, y1, phase=0, chamfer_top=0, chamfer_bot=0, cap=True)`

正 n 棱柱，`af` 是**对边距离**（螺母规格就是这么量的，不是外接圆直径）。
内部按 `r = (af/2)/cos(pi/n)` 换成外接圆。六角螺母就是 `g_ngon(6, 20*s, 0, 5*s, chamfer_top=1.1*s)`。

`phase` 控制转角：默认让一对边朝 ±X，正面看起来是"平头"；需要"尖头"就加 `pi/n`。

---

## 4. `g_fluted(r, y0, y1, flutes=18, amp=0.4, seg_per=7, taper=0.55)`

滚花/齿环（色环、旋钮）。半径按 `r + amp*cos(flutes*a)` 起伏，`seg_per` 是每个齿的平台数（7 足够，太少齿顶会变尖）。
上下端各 0.5mm 收边用 `taper`，避免锐边在高光下发白。
**用 `flat_quad` 保证每个面是硬的**，不要改成 smooth —— 滚花看起来就是要有清晰棱线。

---

## 5. `g_threaded(r, y0, y1, rings=5, amp=0.35, taper=0, seg=56, ...)`

外螺纹颈。用 `(0.5+0.5*cos)^2` 的凸起形状，条纹不会倒扣。`rings` 是螺纹圈数，按照片数；`seg=56` 左右够圆。

---

## 6. `g_bore(r_mouth, r_wall, depth, rings=6, amp=0.45, taper=0.5, seg=64)`

**内螺纹沉孔**，局部 +Y 是钻进方向，孔底已封。做法：`g_revolve(..., flip=True)` + 底部三角扇。
用 `frame(mouth, -outn)` 把它贴到孔口，`outn` 是该面的外法线。

```python
bore = g_bore(5.8, 4.4, 12.0, rings=6, amp=0.35, taper=0.35)
g.extend(bore, T(-1.5, -8.0, Z_F) @ Rx(-90 * D2R))
```

`r_wall` 给成 `r_mouth*0.78` 左右像螺纹底孔；`amp/taper` 都按半径的比例给（不要写死 0.4，小孔会螺纹过粗）。

---

## 7. `emit_face_hole(g, origin, u, v, quad2d, hole_c, hole_r, outward, seg=72)`

在平面上的多边形里挖一个圆孔，输出从孔边到外轮廓的环形三角带。
`quad2d` 是该平面内的二维轮廓（逆时针），`hole_c/hole_r` 也是二维。外轮廓顶点靠 `ray_poly_hit` 从孔心射线求交得到 —— 所以 `quad2d` 必须是**凸**的，凹多边形要拆。

```python
emit_face_hole(g, (0, 0, zf), (1, 0, 0), (0, 1, 0), ins, (-1.5, -8.0), 5.8, (0, 0, 1))
```

**坑 5：孔一定要配一个 `g_bore`。** 只挖洞不装内腔，等于在壳体上开了个能看见内部的洞；正确的"暗孔"= 孔口环 + 内腔 + 封底。

---

## 8. `strip_plate(L, W, thick, holes, lobes, amp, nx, bev)`

长圆板(stadium) + 多个圆孔，用于手轮、压板。沿长轴切 `nx` 条带，每条带按 `holes` 扣除上下缺口，孔壁用同一批弦端点接缝，所以不会裂缝。

- `holes=[(x, r), ...]`，`x` 沿长轴，`r` 孔半径；三个孔的顺序无所谓，内部按位置取最近孔。
- `lobes/amp` 给外缘加余弦齿边（`lobes` 是**半周期数**，视觉齿数≈`lobes`）。
- `bev` 是上下倒角宽度，会生成三带（上倒角 / 直壁 / 下倒角）。

**坑 6：`nx` 太小时孔壁会变成多边形。** 孔径 `r=3` 时 `nx≈76`（每 0.43mm 一段）看起来是圆的。孔半径越大越要加 `nx`。
**坑 7：孔之间的距离要留够。** 两孔边缘间距 < `2*bev` 时，倒角带会互相穿插。

---

## 9. `chamfered_prism(poly, zf, zb, bevel, front_hole=None, side_ports=(), seg=72, port_depth=14.0)`

凸多边形棱柱 + 前后倒角带 + 长侧面，可选正面一个圆孔与若干侧面圆孔。

返回 `[(role, Geo)]`：`role ∈ {'face','side','cavity'}`，直接映射到材质序号。

```python
POLY = [(-33, -20), (33, -20), (33, 6.5), (13, 20), (-13, 20), (-33, 6.5)]  # 逆时针
parts = chamfered_prism(POLY, 20, -20, 1.6,
                        front_hole=((-1.5, -8.0), 5.8),
                        side_ports=[(0, 13.0, 18.4, 7.0), ...])
for role, geo in parts:
    body.add(MAT[role], geo)
```

- `poly` 必须是凸的、逆时针（`inset_poly` 是 miter 内缩，凹角会炸）。
- `side_ports` 的 `edge_index` 对应 `poly` 的边 `i → i+1`；`t` 是从首顶点沿边的距离；`zd` 是孔心在该边深度方向上的位置（从**后**倒角边起算）。
- 侧面孔会自动附带一个 `cavity` 内腔，不需要自己再装。

**坑 8：边索引弄错一个，孔就跑到别的面上。** 建模前先把 `poly` 的边列成表（每条边的两个端点、长度、朝向），对着照片点。
**坑 9：`zd` 是"从后倒角边起算"，不是从中心。** 深度 = `(zf-bevel)-(zb+bevel)`，孔心别给到边界外。

---

## 10. 材质与贴图

```python
STEEL    = mat("StainlessSteel",     (0.90, 0.905, 0.915), 0.95, 0.30)
MACHINED = mat("StainlessMachined", (0.89, 0.895, 0.905), 0.95, 0.36)
CAVITY   = mat("ThreadCavity",      (0.30, 0.31, 0.33),   0.70, 0.50)
PAPER    = mat("NamePlate",         (1.0, 1.0, 1.0),      0.0,  0.55, tex=True)
```

- `mat()` 内部已经把 sRGB 取 `pow 2.2` 转成线性，**传入的是 sRGB 观感值**。
- 金属件 `metallic≈0.95`，`roughness 0.30~0.40`；抛太亮（<0.2）在无环境时像塑料。
- `tex=True` 表示该材质需要 `TEXCOORD_0`，几何必须带 UV（`flat_quad(..., uv=...)`）。
- 贴图方向：glTF 的 `v=0` 在图像**顶部**，写 UV 时别按 OpenGL 习惯翻转。

铭牌做法：用 Pillow 的 `Image.transform(size, Image.QUAD, data=quad)` 从照片里透视校正出一张矩形贴图，再做亮度归一化去掉现场光照；挂在 `NamePlate` 子节点上，用一个略高于本体的薄板 + 边缘倒角。

---

## 11. 装配与导出

```python
root = Node("Part")
body = Node("Body");            root.children.append(body)
valve = Node("ValveCentre", T(0, 19, 0));  body.children.append(valve)

write_glb("out.glb", root, [STEEL, MACHINED, CAVITY, ...],
          texture_png=label_png, root_matrix=MM_TO_M)
```

- 节点矩阵用 4x4 `numpy`：`T(x,y,z) @ Rx(deg*D2R) @ Rz(...) @ S(s)`。
- 斜置件统一用 `Rx/Rz` 转到目标朝向，再在局部空间里按 Y 轴堆：底座 → 鼓 → 螺母 → 螺纹 → 色环 → 压帽 → 阀杆 → 手轮。
- 根节点用 `MM_TO_M`（0.001）转米制；**不要**在子节点里重复缩放。
- 导出后立刻检查：`verify_glb(path)` 返回 `bytes/nodes/meshes/materials/tris/problems`，`problems` 必须为空，`tris` 不能为 0。
