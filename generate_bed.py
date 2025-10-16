import math
from pathlib import Path


num_atoms = 100
atom_type = 1

DIAMETER = 1.7e-4
DENSITY = 2650.0

xlo, xhi = 0.0, 0.05
ylo, yhi = 0.0, 0.02
zlo, zhi = 0.0, 0.01

bed_x_hi = 0.05


R = DIAMETER / 2.0
margin = R
x_min = xlo + margin
x_max = bed_x_hi - margin
y_min = ylo + margin
y_max = yhi - margin
z_val = zlo + R

spacing = DIAMETER * 1.05

Nx = max(1, int(math.floor((x_max - x_min) / spacing)) + 1)
Ny = max(1, int(math.floor((y_max - y_min) / spacing)) + 1)

xs = [x_min + i * spacing for i in range(Nx) if x_min + i * spacing <= x_max]
ys = [y_min + j * spacing for j in range(Ny) if y_min + j * spacing <= y_max]

coords = []
for j, y in enumerate(ys):
   
    x_offset = 0.0
    for i, x in enumerate(xs):
        x_pos = x + x_offset
        if x_pos > x_max:
            continue
        coords.append((x_pos, y, z_val))
        if len(coords) >= num_atoms:
            break
    if len(coords) >= num_atoms:
        break

if len(coords) < num_atoms:
    extra_needed = num_atoms - len(coords)
    spacing2 = DIAMETER * 1.02
    xs2 = [x_min + i * spacing2 for i in range(int(math.floor((x_max - x_min) / spacing2)) + 1)
           if x_min + i * spacing2 <= x_max]
    ys2 = [y_min + j * spacing2 for j in range(int(math.floor((y_max - y_min) / spacing2)) + 1)
           if y_min + j * spacing2 <= y_max]
    for y in ys2:
        for x in xs2:
            coords.append((x, y, z_val))
            if len(coords) >= num_atoms:
                break
        if len(coords) >= num_atoms:
            break

assert len(coords) >= num_atoms, f"Only {len(coords)} positions available"
coords = coords[:num_atoms]

out_path = Path("./IC_uniform_100.in")
with out_path.open("w") as f:
    f.write(" Granular Flow simulation\n\n")
    f.write(f"{num_atoms:12d}  atoms\n")
    f.write(f"{1:12d}  atom types\n\n")
    f.write(f"{xlo:11.6f}  {xhi:11.6f}  xlo xhi\n")
    f.write(f"{ylo:11.6f}  {yhi:11.6f}  ylo yhi\n")
    f.write(f"{zlo:11.6f}  {zhi:11.6f}  zlo zhi\n\n")
    f.write(" Atoms\n\n")
    for idx, (x, y, z) in enumerate(coords, start=1):
        f.write(f"{idx:6d}   {atom_type:d} {DIAMETER:.5E}  {DENSITY:.5E}  {x:.5f}  {y:.5f}  {z:.5f}\n")

print(f"Wrote {num_atoms} atoms to {out_path}")
