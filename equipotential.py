import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import sys

# ── Usage ─────────────────────────────────────────────────────────────────────
# python3 equipotential.py                  → uses default E_col = 2747 cm-1
# python3 equipotential.py 1800             → uses E_col = 1800 cm-1
# ─────────────────────────────────────────────────────────────────────────────

E_col = float(sys.argv[1]) if len(sys.argv) > 1 else 2747.0

# ── 1. Load PES file ──────────────────────────────────────────────────────────
data = []
with open('kpdata_ne.pot', 'r') as f:
    for line in f.readlines()[2:]:        # skip 2 header lines
        parts = line.split()
        if len(parts) == 4:
            data.append(list(map(float, parts)))
data = np.array(data)

r_vals     = np.unique(data[:, 0])
R_vals     = np.unique(data[:, 1])
theta_vals = np.unique(data[:, 2])

# ── 2. Build 3D array V[r_idx, R_idx, theta_idx] ─────────────────────────────
V = np.zeros((len(r_vals), len(R_vals), len(theta_vals)))
for row in data:
    ri = np.where(r_vals     == row[0])[0][0]
    Ri = np.where(R_vals     == row[1])[0][0]
    ti = np.where(theta_vals == row[2])[0][0]
    V[ri, Ri, ti] = row[3]

r_eq_idx = np.where(r_vals == 5.873)[0][0]   # r_eq = 5.873 bohr

# ── 3. Find turning point R_turn(theta) for given E_col ──────────────────────
def find_turning_point(E, r_idx):
    """
    For each theta, find R where V(R, r_vals[r_idx], theta) = E
    by linear interpolation between bracketing grid points.
    Returns dict {theta: R_turn} for thetas where a crossing exists.
    """
    turning = {}
    for ti, theta in enumerate(theta_vals):
        v_vs_R = V[r_idx, :, ti]
        for Ri in range(len(R_vals) - 1):
            if v_vs_R[Ri] >= E >= v_vs_R[Ri + 1]:
                frac   = (E - v_vs_R[Ri]) / (v_vs_R[Ri + 1] - v_vs_R[Ri])
                R_turn = R_vals[Ri] + frac * (R_vals[Ri + 1] - R_vals[Ri])
                turning[theta] = R_turn
                break
    return turning

turning = find_turning_point(E_col, r_eq_idx)

if 0.0 not in turning or 90.0 not in turning:
    print(f"ERROR: no turning point found at theta=0 or theta=90 for E={E_col}")
    print(f"Available turning points: {list(turning.keys())}")
    sys.exit(1)

a = turning[0.0]    # end-on semi-axis
c = turning[90.0]   # side-on semi-axis

print(f"E_col = {E_col:.1f} cm-1")
print(f"a (end-on,  theta=0°)  = {a:.4f} bohr")
print(f"c (side-on, theta=90°) = {c:.4f} bohr")
print(f"a/c = {a/c:.4f}")
print()
print("All turning points:")
for theta, R_t in turning.items():
    print(f"  theta={theta:5.1f}  R_turn={R_t:.4f} bohr")

# ── 4. Convert turning points to Cartesian for equipotential curve ────────────
# x = R_turn * cos(theta)  (along bond axis)
# y = R_turn * sin(theta)  (perpendicular)
# reflect to get full closed curve

thetas_rad = np.radians(list(turning.keys()))
R_turns    = np.array(list(turning.values()))

# sort by theta
sort_idx   = np.argsort(thetas_rad)
thetas_rad = thetas_rad[sort_idx]
R_turns    = R_turns[sort_idx]

# build full 360-degree curve by reflecting
theta_full = np.concatenate([
    thetas_rad,                        # 0 to 90
    np.pi - thetas_rad[::-1][1:],      # 90 to 180 (mirror x)
    np.pi + thetas_rad[1:],            # 180 to 270
    2*np.pi - thetas_rad[::-1][1:-1],  # 270 to 360
])
R_full = np.concatenate([
    R_turns,
    R_turns[::-1][1:],
    R_turns[1:],
    R_turns[::-1][1:-1],
])

x_equip = R_full * np.cos(theta_full)
y_equip = R_full * np.sin(theta_full)

# ── 5. Plot ───────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(12, 5))
gs  = gridspec.GridSpec(1, 2, width_ratios=[1, 1.1], wspace=0.35)
ax1 = fig.add_subplot(gs[0])   # equipotential
ax2 = fig.add_subplot(gs[1])   # V vs R

# --- Left panel: equipotential ---
ax1.fill(x_equip, y_equip, color='#cce4f7', alpha=0.5, zorder=1)
ax1.plot(x_equip, y_equip, color='steelblue', linewidth=1.8, zorder=2,
         label=f'Equipotential  V = {E_col:.0f} cm⁻¹')

# Li₂ bond and atoms
r_eq = r_vals[r_eq_idx]
ax1.plot([-r_eq/2, r_eq/2], [0, 0], 'k-', linewidth=2.5, zorder=3)
for xpos, label in [(-r_eq/2, 'Li'), (r_eq/2, 'Li')]:
    ax1.add_patch(plt.Circle((xpos, 0), 0.25, color='#6baed6', zorder=4))
    ax1.text(xpos, 0, label, ha='center', va='center', fontsize=8, zorder=5)

# Ne atoms at turning points
ax1.add_patch(plt.Circle((a, 0), 0.2, color='#fc8d59', zorder=4))
ax1.text(a, 0, 'Ne', ha='center', va='center', fontsize=7, zorder=5)
ax1.add_patch(plt.Circle((0, c), 0.2, color='#fc8d59', zorder=4))
ax1.text(0, c, 'Ne', ha='center', va='center', fontsize=7, zorder=5)

# a and c arrows
ax1.annotate('', xy=(a, 0), xytext=(0, 0),
             arrowprops=dict(arrowstyle='->', color='steelblue', lw=1.5))
ax1.annotate('', xy=(0, c), xytext=(0, 0),
             arrowprops=dict(arrowstyle='->', color='firebrick', lw=1.5))

ax1.text(a/2, 0.35, f'a = {a:.2f} bohr', ha='center', color='steelblue', fontsize=9)
ax1.text(0.4, c/2,  f'c = {c:.2f} bohr', ha='left',   color='firebrick',  fontsize=9)
ax1.text(a + 0.3, 0.1, 'end-on',  color='steelblue', fontsize=8)
ax1.text(0.2, c + 0.3, 'side-on', color='firebrick',  fontsize=8)

# reference axes
lim = a * 1.25
ax1.axhline(0, color='gray', linewidth=0.5, linestyle='--', alpha=0.5)
ax1.axvline(0, color='gray', linewidth=0.5, linestyle='--', alpha=0.5)
ax1.set_xlim(-lim, lim)
ax1.set_ylim(-lim * 0.7, lim * 0.7)
ax1.set_aspect('equal')
ax1.set_xlabel('x  (bohr) — along bond axis')
ax1.set_ylabel('y  (bohr) — perpendicular')
ax1.set_title(f'Equipotential surface  (E = {E_col:.0f} cm⁻¹)')
ax1.legend(fontsize=8, loc='upper right')

# --- Right panel: V vs R ---
ti0  = np.where(theta_vals == 0.)[0][0]
ti90 = np.where(theta_vals == 90.)[0][0]

V0   = V[r_eq_idx, :, ti0]    # V at theta=0
V90  = V[r_eq_idx, :, ti90]   # V at theta=90

# clip for readability
V_max = max(E_col * 3, 1000)
mask0  = V0  < V_max
mask90 = V90 < V_max

ax2.plot(R_vals[mask0],  V0[mask0],  'o-', color='steelblue', markersize=4,
         linewidth=1.8, label='θ = 0°  (end-on)')
ax2.plot(R_vals[mask90], V90[mask90], 's-', color='firebrick', markersize=4,
         linewidth=1.8, label='θ = 90°  (side-on)')

# E_col horizontal line
ax2.axhline(E_col, color='#d4a017', linewidth=1.5, linestyle='--',
            label=f'E_col = {E_col:.0f} cm⁻¹')

# turning point vertical lines
if 0.0 in turning:
    ax2.axvline(a, color='steelblue', linewidth=1, linestyle=':')
    ax2.plot(a, E_col, 'o', color='steelblue', markersize=8, zorder=5)
    ax2.text(a + 0.05, E_col * 1.05, f'a = {a:.2f}', color='steelblue', fontsize=8)

if 90.0 in turning:
    ax2.axvline(c, color='firebrick', linewidth=1, linestyle=':')
    ax2.plot(c, E_col, 's', color='firebrick', markersize=8, zorder=5)
    ax2.text(c + 0.05, E_col * 1.05, f'c = {c:.2f}', color='firebrick', fontsize=8)

ax2.set_xlabel('R  (bohr)')
ax2.set_ylabel('V  (cm⁻¹)')
ax2.set_title(f'V(R, r_eq, θ)  →  turning points at E = {E_col:.0f} cm⁻¹')
ax2.set_xlim(R_vals[0], R_vals[-3])
ax2.set_ylim(-200, V_max)
ax2.legend(fontsize=8)
ax2.grid(True, alpha=0.2)

fig.suptitle(f'Li₂–Ne equipotential at E_col = {E_col:.0f} cm⁻¹'
             f'     a = {a:.3f} bohr,  c = {c:.3f} bohr,  a/c = {a/c:.3f}',
             fontsize=10)

plt.savefig(f'equipotential_{int(E_col)}cm.png', dpi=150, bbox_inches='tight')
plt.close()
print(f"\nSaved equipotential_{int(E_col)}cm.png")
