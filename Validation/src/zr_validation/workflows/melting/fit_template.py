#!/usr/bin/env python3


MELT_ANALYSIS_FUNCTIONS = r'''def parse_rdf(path):
    """
    读取 LAMMPS fix ave/time mode vector 的 RDF 文件。
    返回最后一个 block 的 r 和 g(r)。
    """
    with open(path) as f:
        lines = [ln for ln in f if not ln.lstrip().startswith("#")]

    blocks = []
    i = 0

    while i < len(lines):
        parts = lines[i].split()

        if len(parts) == 2:
            try:
                nrows = int(parts[1])
            except ValueError:
                i += 1
                continue

            data = []
            for j in range(i + 1, min(i + 1 + nrows, len(lines))):
                p = lines[j].split()
                if len(p) >= 3:
                    data.append((float(p[1]), float(p[2])))

            if data:
                blocks.append(np.array(data))

            i += 1 + nrows
        else:
            i += 1

    if not blocks:
        raise ValueError(f"无法解析 RDF 文件: {path}")

    arr = blocks[-1]
    return arr[:, 0], arr[:, 1]


def _cross(x0, y0, x1, y1, level):
    if y1 == y0:
        return x0
    return x0 + (level - y0) * (x1 - x0) / (y1 - y0)


def quarter_max_width(r, g, peak_search_max=PEAK_SEARCH_MAX, frac=QUARTER_FRAC):
    """
    计算 RDF 第一峰在 frac * peak height 处的全宽。
    默认 frac = 0.25，即四分之一最大值全宽。
    """
    r = np.asarray(r, float)
    g = np.asarray(g, float)

    region = np.where(r <= peak_search_max)[0]
    if len(region) == 0:
        raise ValueError("PEAK_SEARCH_MAX 太小，RDF 中没有可搜索的点。")

    ip = region[int(np.argmax(g[region]))]
    level = frac * g[ip]

    rl = r[0]
    for i in range(ip, 0, -1):
        if (g[i] - level) * (g[i - 1] - level) <= 0:
            rl = _cross(r[i - 1], g[i - 1], r[i], g[i], level)
            break

    rr = r[-1]
    for i in range(ip, len(r) - 1):
        if (g[i] - level) * (g[i + 1] - level) <= 0:
            rr = _cross(r[i], g[i], r[i + 1], g[i + 1], level)
            break

    return rr - rl


def moving_average(y, window):
    """
    居中滑动平均。
    边缘窗口自动收缩，避免平台值被 padding 影响。
    """
    y = np.asarray(y, float)
    half = window // 2
    n = len(y)
    out = np.empty(n)

    for i in range(n):
        lo = max(0, i - half)
        hi = min(n, i + half + 1)
        out[i] = y[lo:hi].mean()

    return out


def linear_fit(T, y, window):
    lo, hi = window
    m = (T >= lo) & (T <= hi)

    if m.sum() < 2:
        raise ValueError(
            f"拟合窗口 {window} K 内只有 {int(m.sum())} 个数据点；"
            f"请调整窗口或扩大扫描范围。当前数据范围: {T.min():.0f}-{T.max():.0f} K。"
        )

    a, b = np.polyfit(T[m], y[m], 1)
    return a, b


def vertical_distance(T, y, a, b):
    """
    点到直线 y = a*T + b 的竖直距离。
    """
    return np.abs(y - (a * T + b))


def gaussian(x, A, mu, sigma, c):
    return A * np.exp(-(x - mu) ** 2 / (2.0 * sigma ** 2)) + c


def analyze():
    from scipy.optimize import curve_fit

    temps = temp_list()
    T_list = []
    fwqm = []

    for T in temps:
        path = rdf_path(T)
        if not os.path.exists(path):
            print(f"[skip] 缺少 {path}")
            continue

        try:
            r, g = parse_rdf(path)
            T_list.append(T)
            fwqm.append(quarter_max_width(r, g))
        except Exception as e:
            print(f"[warn] T={T:.0f} K 解析失败: {e}")

    if len(T_list) != len(temps):
        sys.exit(
            f"只有 {len(T_list)}/{len(temps)} 个温度的数据；拒绝对不完整的扫描给出熔点。\n"
            f"请检查 {SCAN_DIR}/ 下的 RDF 文件和 LAMMPS log。"
        )

    T = np.array(T_list)
    fwqm = np.array(fwqm)

    order = np.argsort(T)
    T = T[order]
    fwqm = fwqm[order]

    fwqm_s = moving_average(fwqm, SMOOTH_WINDOW)

    a_lo, b_lo = linear_fit(T, fwqm_s, LOW_T_WINDOW)
    a_hi, b_hi = linear_fit(T, fwqm_s, HIGH_T_WINDOW)

    dist = np.minimum(
        vertical_distance(T, fwqm_s, a_lo, b_lo),
        vertical_distance(T, fwqm_s, a_hi, b_hi),
    )

    p0 = [
        dist.max(),
        T[int(np.argmax(dist))],
        (T.max() - T.min()) / 10.0,
        0.0,
    ]

    try:
        popt, _ = curve_fit(gaussian, T, dist, p0=p0, maxfev=20000)
        Tm = popt[1]
    except Exception as e:
        print(f"[warn] 高斯拟合失败: {e}")
        print("       退回使用 min-distance 曲线最大值位置作为 Tm。")
        popt = None
        Tm = T[int(np.argmax(dist))]

    print("\n================ 结果 ================")
    print(f"低温拟合线 {LOW_T_WINDOW} K:  y = {a_lo:.3e} * T + {b_lo:.4f}")
    print(f"高温拟合线 {HIGH_T_WINDOW} K: y = {a_hi:.3e} * T + {b_hi:.4f}")
    print(f"估算熔点 Tm = {Tm:.1f} K")
    print("======================================\n")

    np.savetxt(
        OUT_CSV,
        np.column_stack([T, fwqm, fwqm_s, dist]),
        header="T[K]  FWQM_raw[A]  FWQM_smooth[A]  min_dist[A]",
        fmt="%.6g",
        comments="",
    )
    print(f"[save] {OUT_CSV}")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 8), sharex=True)
    tline = np.array([T.min(), T.max()])

    ax1.plot(T, fwqm, ".", ms=4, alpha=0.4, label="raw FWQM")
    ax1.plot(T, fwqm_s, "-", lw=2, label=f"{SMOOTH_WINDOW}-pt smooth")
    ax1.plot(tline, a_lo * tline + b_lo, "--", label="low-T fit")
    ax1.plot(tline, a_hi * tline + b_hi, "--", label="high-T fit")
    ax1.axvline(Tm, color="k", ls=":", lw=1)
    ax1.set_ylabel("first-peak quarter-max width [Å]")
    ax1.set_title("RDF first-peak broadening vs T")
    ax1.legend(fontsize=8)

    ax2.plot(T, dist, ".", ms=4, alpha=0.5, label="min distance to fits")

    if popt is not None:
        tt = np.linspace(T.min(), T.max(), 800)
        ax2.plot(tt, gaussian(tt, *popt), "-", lw=2, label="Gaussian fit")

    ax2.axvline(Tm, color="k", ls=":", lw=1, label=f"$T_m$ = {Tm:.0f} K")
    ax2.set_xlabel("Temperature [K]")
    ax2.set_ylabel("min distance [Å]")
    ax2.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(OUT_FIG, dpi=150)
    print(f"[save] {OUT_FIG}")
    result = {
        "method": "RDF first-peak quarter-maximum width / Gaussian heuristic",
        "Tm_estimate_K": float(Tm),
        "fit_status": "gaussian" if popt is not None else "peak-fallback",
        "temperature_points": len(T),
        "T_min_K": float(T.min()),
        "T_max_K": float(T.max()),
        "low_T_window": list(LOW_T_WINDOW),
        "high_T_window": list(HIGH_T_WINDOW),
        "note": "RDF-derived transition estimate, not a rigorous coexistence melting point",
    }
    with open("melt_result.json", "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)
    print("[save] melt_result.json")


'''

