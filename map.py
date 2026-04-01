
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np

# 读取数据
df = pd.read_csv("C:/Users/zhangqi/Desktop/map-gutMAG.csv")
world = gpd.read_file("C:/Users/zhangqi/Desktop/map-country/ne_110m_admin_0_countries.shp")

# 国家名匹配
name_map = {
    "USA": "United States of America",
    "UK": "United Kingdom",
    "Czech": "Czechia",
    "Korea": "South Korea"
}
df["country"] = df["country"].replace(name_map)

# 合并
world = world.merge(df, how="left", left_on="ADMIN", right_on="country")

# log10 转换
world["log_risk"] = world["size"]

# 色标范围
vmin = world["log_risk"].quantile(0.02)
vmax = world["log_risk"].quantile(0.98)

# 作图
fig, ax = plt.subplots(1, 1, figsize=(18, 8))

# 无数据国家底图
world.plot(
    ax=ax,
    color="#f2f2f2",
    edgecolor="none"
)

# 有数据国家
world.dropna(subset=["log_risk"]).plot(
    column="log_risk",
    cmap="inferno",
    linewidth=0,
    ax=ax,
    legend=True,
    vmin=vmin,
    vmax=vmax,
    legend_kwds={"label": "log10(ARG risk)", "shrink": 0.75}
)

# -------------------------
# 加坐标轴与边框
# -------------------------
ax.set_xlim(-180, 180)
ax.set_ylim(-90, 90)

# 经度刻度
ax.set_xticks([-180, -120, -60, 0, 60, 120, 180])
ax.set_xticklabels(
    ["180°W", "120°W", "60°W", "0°", "60°E", "120°E", "180°E"],
    fontsize=11
)

# 纬度刻度
ax.set_yticks([-60, -30, 0, 30, 60])
ax.set_yticklabels(
    ["60°S", "30°S", "0°", "30°N", "60°N"],
    fontsize=11
)

# 坐标轴标题
ax.set_xlabel("Longitude", fontsize=13)
ax.set_ylabel("Latitude", fontsize=13)

# 外边框
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(1.0)
    spine.set_color("black")

# 刻度线
ax.tick_params(axis='both', which='both', length=4, width=0.8, color='black')

# 可选：浅灰经纬网
ax.grid(True, linestyle="--", linewidth=0.4, color="lightgray", alpha=0.7)

plt.tight_layout()
plt.show()