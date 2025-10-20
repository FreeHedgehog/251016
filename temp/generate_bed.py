import math
from pathlib import Path

# 若设为 None，则按 bed_height/num_layers 填满；否则按目标数量截断
num_atoms = 5000  # None 或整数
atom_type = 1

DIAMETER = 1.7e-4
DENSITY = 2650.0

xlo, xhi = 0.0, 0.05
ylo, yhi = 0.0, 0.02
zlo, zhi = 0.0, 0.01

bed_x_hi = xhi

bed_height = 0.004
num_layers = None

R = DIAMETER / 2.0
x_min = xlo + R
x_max = bed_x_hi - R
y_min = ylo + R
y_max = yhi - R

# 层内近似六角紧密：x 间距 = 直径；y 间距 = 直径 * sqrt(3)/2
DX = DIAMETER
DY = DIAMETER * math.sqrt(3) / 2

Nx = max(1, int(math.floor((x_max - x_min) / DX)) + 1)
Ny = max(1, int(math.floor((y_max - y_min) / DY)) + 1)

# 层间距：每层相隔一个直径，避免重叠
DZ = DIAMETER

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

out_path = Path("./initial_bed.lammps")
with out_path.open("w") as f:
    # 标题行（注释），避免触发 Unknown identifier
    f.write("LAMMPS data file via generate_bed.py\n\n")
    # 计数与边界（标准头部）
    f.write(f"{len(coords)} atoms\n")
    f.write(f"{1} atom types\n\n")
    f.write(f"{xlo} {xhi} xlo xhi\n")
    f.write(f"{ylo} {yhi} ylo yhi\n")
    f.write(f"{zlo} {zhi} zlo zhi\n\n")
    # Atoms 段，明确风格为 sphere；字段顺序必须为：id type x y z diameter density
    f.write("Atoms # sphere\n\n")
    for idx, (x, y, z) in enumerate(coords, start=1):
        f.write(f"{idx} {atom_type} {x:.6f} {y:.6f} {z:.6f} {DIAMETER:.6e} {DENSITY:.6e}\n")

print(f"Wrote {len(coords)} atoms to {out_path}")
print(f"Layers used: {min(max_layers, (len(coords) // max(1, Nx*Ny)) + (1 if len(coords) % max(1, Nx*Ny) else 0))}, Nx={Nx}, Ny={Ny}, DX={DX:.6e}, DY={DY:.6e}, DZ={DZ:.6e}")
