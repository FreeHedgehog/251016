import math
from pathlib import Path

# 目标原子数（若想按床厚度完全填满，可设为 None）
num_atoms = 5000  # 或将其设为 None 以按 bed_height 自动填充
atom_type = 1

DIAMETER = 1.7e-4
DENSITY = 2650.0

# 模拟区域边界（需与 OpenFOAM/LAMMPS 案例一致）
xlo, xhi = 0.0, 0.05
ylo, yhi = 0.0, 0.02
zlo, zhi = 0.0, 0.01

# 沙床在 x 方向铺到哪里（通常就是到出口）
bed_x_hi = xhi

# 沙床厚度（从底部开始向上，单位米）。
# 也可用层数控制：将 num_layers 设为整数会覆盖 bed_height。
bed_height = 0.004  # 例如铺 4mm 厚的床
num_layers = None   # 例如设为 5 覆盖 bed_height；保持 None 则按 bed_height 计算

R = DIAMETER / 2.0
margin = R
x_min = xlo + margin
x_max = bed_x_hi - margin
y_min = ylo + margin
y_max = yhi - margin

# 层内采用 2D 六角紧密近似：x 间距 = 直径，y 间距 = 直径 * sqrt(3)/2，行间交错半个直径
DX = DIAMETER
DY = DIAMETER * math.sqrt(3) / 2

Nx = max(1, int(math.floor((x_max - x_min) / DX)) + 1)
Ny = max(1, int(math.floor((y_max - y_min) / DY)) + 1)

# 层间 z 间距：用简单立方堆叠避免重叠（每层相隔 1 个直径）
DZ = DIAMETER

# 计算可铺到的 z 上限（由床厚决定，或由几何上限决定）
z_top_limit = zlo + (bed_height if bed_height is not None else (zhi - zlo))
if num_layers is None:
    max_layers = max(1, int(math.floor((z_top_limit - (zlo + R)) / DZ)) + 1)
else:
    max_layers = max(1, num_layers)

coords = []
for k in range(max_layers):
    z = zlo + R + k * DZ
    if z + R > z_top_limit + 1e-12:
        break
    for j in range(Ny):
        y = y_min + j * DY
        if y > y_max + 1e-12:
            break
        # 行间交错半个直径，层内形成近似六角紧密排列
        x_offset = 0.5 * DX if (j % 2 == 1) else 0.0
        for i in range(Nx):
            x = x_min + i * DX + x_offset
            if x > x_max + 1e-12:
                break
            coords.append((x, y, z))
            if num_atoms is not None and len(coords) >= num_atoms:
                break
        if num_atoms is not None and len(coords) >= num_atoms:
            break
    if num_atoms is not None and len(coords) >= num_atoms:
        break

# 写到案例根目录，便于 in.lammps 直接读取
out_path = Path("./initial_bed.lammps")
with out_path.open("w") as f:
    f.write(" Granular Flow simulation\n\n")
    f.write(f"{len(coords):12d}  atoms\n")
    f.write(f"{1:12d}  atom types\n\n")
    f.write(f"{xlo:11.6f}  {xhi:11.6f}  xlo xhi\n")
    f.write(f"{ylo:11.6f}  {yhi:11.6f}  ylo yhi\n")
    f.write(f"{zlo:11.6f}  {zhi:11.6f}  zlo zhi\n\n")
    f.write(" Atoms\n\n")
    for idx, (x, y, z) in enumerate(coords, start=1):
        f.write(f"{idx:6d}   {atom_type:d} {DIAMETER:.5E}  {DENSITY:.5E}  {x:.6f}  {y:.6f}  {z:.6f}\n")

print(f"Wrote {len(coords)} atoms to {out_path}")
print(f"Layers used: {min(max_layers, (len(coords) // max(1, Nx*Ny)) + (1 if len(coords) % max(1, Nx*Ny) else 0))}, Nx={Nx}, Ny={Ny}, DX={DX:.6e}, DY={DY:.6e}, DZ={DZ:.6e}")
