"""Merge grid features (Model 2) with ERA5 climate (Model 1) and draw maps. usage: python3 -I merge_and_maps.py OUT CLIMATE_CSV BUILDINGS_GEOJSON"""
import sys, os
import numpy as np, pandas as pd, geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

out, climate_fp, bld_fp = sys.argv[1:4]

# ---------------- merge: grid features + climate summary ----------------
clim = pd.read_csv(climate_fp, parse_dates=["Date"])
summ = {
    "clim_temp_mean_C": clim["Temperature_C"].mean(),
    "clim_wind_mean_ms": clim["WindSpeed_ms"].mean(),
    "clim_dewpoint_mean_C": clim["DewPoint_C"].mean(),
    "clim_CSI_partial_mean": clim["CSI_partial"].mean(),
    "clim_CSI_partial_max": clim["CSI_partial"].max(),
    "clim_CSI_partial_premonsoon_mean": clim.loc[clim["Month"].isin([4, 5, 6]), "CSI_partial"].mean(),
}
for cell in (500, 250):
    g = pd.read_csv(os.path.join(out, f"grid_fractal_{cell}m.csv"))
    for k, v in summ.items():
        g[k] = round(v, 4)
    g.to_csv(os.path.join(out, f"Data_master_{cell}m.csv"), index=False)

# ---------------- maps (500 m) ----------------
INK, MUTED = "#1f2937", "#6b7280"
grid = gpd.read_file(os.path.join(out, "grid_fractal_500m.geojson"))
roads = gpd.read_file(os.path.join(out, "roads_classified.geojson"))
bld = gpd.read_file(bld_fp)
# plot in UTM so the map is not stretched
for gdf_name in ("grid", "roads", "bld"):
    globals()[gdf_name] = globals()[gdf_name].to_crs(32645)

THICK_C, THIN_C = "#8c1c13", "#e08a2b"

# show only what the analysis actually covers: clip everything to the grid extent
from shapely.geometry import box as _box
gx0, gy0, gx1, gy1 = grid.total_bounds
extent = _box(gx0, gy0, gx1, gy1)
bld = bld.clip(extent)
roads = roads.clip(extent)

# Fig 1: thick vs thin road extraction (mentor's key figure)
fig, ax = plt.subplots(figsize=(9, 8.2), dpi=200)
bld.plot(ax=ax, color="#d7dce2", linewidth=0)
grid.boundary.plot(ax=ax, color="#9aa3af", linewidth=0.6, linestyle=(0, (4, 3)))
roads[roads.road_class == "thin"].plot(ax=ax, color=THIN_C, linewidth=0.8)
roads[roads.road_class == "thick"].plot(ax=ax, color=THICK_C, linewidth=2.6)
for _, r in grid.iterrows():
    c = r.geometry.centroid
    x0, y0, x1, y1 = r.geometry.bounds
    ax.text(x0 + 12, y1 - 12, r.grid_id[1:], fontsize=6, color="#374151", ha="left", va="top",
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.8))
ax.set_axis_off()
ax.set_title("Road extraction: thick vs thin roads — Central Kolkata (OSM)", fontsize=12.5, color=INK, fontweight="bold")
ax.legend(handles=[
    Line2D([0], [0], color=THICK_C, lw=2.6, label="Thick road (trunk / primary / secondary)"),
    Line2D([0], [0], color=THIN_C, lw=0.9, label="Thin road (tertiary / residential / service / other)"),
    plt.Rectangle((0, 0), 1, 1, fc="#d7dce2", label="Building footprints"),
], loc="upper center", bbox_to_anchor=(0.5, -0.01), ncol=3, frameon=False, fontsize=8)
ax.text(0.0, -0.075, "Data: © OpenStreetMap contributors. Dashed squares = 500 m analysis grid (label = row/col). North is up.", transform=ax.transAxes, fontsize=7, color=MUTED)
fig.savefig(os.path.join(out, "map_thin_thick_roads.png"), bbox_inches="tight")
plt.close(fig)


# Fig 2: RFD and BFD choropleths (single-hue sequential, value annotated; flagged cells hatched)
def choropleth(ax, col, flagcol, cmap, title, label):
    g = grid.copy()
    vmin, vmax = np.nanmin(g[col]), np.nanmax(g[col])
    g.plot(ax=ax, column=col, cmap=cmap, vmin=vmin, vmax=vmax, edgecolor="white", linewidth=1.2,
           missing_kwds={"color": "#eef0f3"})
    flagged = g[g[flagcol].astype(bool) | g[col].isna()]
    flagged.plot(ax=ax, facecolor="none", edgecolor="#6b7280", hatch="///", linewidth=0)
    norm = plt.Normalize(vmin, vmax)
    for _, r in g.iterrows():
        c = r.geometry.centroid
        if np.isnan(r[col]):
            txt, colr = "n/a", MUTED
        else:
            txt = f"{r[col]:.2f}"
            colr = "white" if norm(r[col]) > 0.55 else INK
        ax.text(c.x, c.y, txt, ha="center", va="center", fontsize=7.5, color=colr, fontweight="bold")
    ax.set_axis_off()
    ax.set_title(title, fontsize=11.5, color=INK, loc="left", fontweight="bold")
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = plt.colorbar(sm, ax=ax, fraction=0.04, pad=0.02, shrink=0.8)
    cb.set_label(label, fontsize=8, color=MUTED)
    cb.outline.set_visible(False)


fig, axes = plt.subplots(1, 2, figsize=(15, 6.4), dpi=200)
choropleth(axes[0], "RFD", "low_data_road", "YlOrBr", "Road Fractal Dimension (RFD)", "RFD (box-counting)")
choropleth(axes[1], "BFD", "low_data_built", "Blues", "Built-up Fractal Dimension (BFD)", "BFD (box-counting)")
fig.text(0.01, 0.01, "500 m grid. Hatched cells = too little road / building data for a reliable value (RFD < 1000 m of road; BFD < 20 buildings). "
         "Data: © OpenStreetMap contributors.", fontsize=7, color=MUTED)
fig.tight_layout(rect=[0, 0.03, 1, 1])
fig.savefig(os.path.join(out, "map_RFD_BFD.png"), bbox_inches="tight")
plt.close(fig)
print("done")
