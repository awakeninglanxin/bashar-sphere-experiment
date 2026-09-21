``````
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
铜球位置共振实验 · 数据分析工具包
依赖: numpy scipy matplotlib pandas (均可用 pip 安装)

用法:
    python 铜球实验分析工具.py --demo      # 跑一遍全部自检，生成示例图
    或按需 import 其中的函数使用
"""
import numpy as np
from scipy.optimize import curve_fit
import argparse, sys

# ============================================================
# 1. ring-down 拟合：从衰减波形同时提取 f0 与 Q
# ============================================================
def damped_sine(t, A, f0, Q, phi, C):
    """阻尼正弦: A*exp(-pi*f0*t/Q)*sin(2*pi*f0*t+phi)+C"""
    return A * np.exp(-np.pi * f0 * t / Q) * np.sin(2 * np.pi * f0 * t + phi) + C

def fit_ringdown(t, s, f_guess=None, Q_guess=200):
    """
    输入: t (s), s (任意单位电压)
    返回: dict(f0, Q, A, phi, C, f0_err, Q_err, r2)
    """
    t = np.asarray(t, float); s = np.asarray(s, float)
    dt = np.median(np.diff(t))
    if f_guess is None:
        # 用 FFT 粗估
        S = np.abs(np.fft.rfft(s - s.mean()))
        fr = np.fft.rfftfreq(len(s), dt)
        f_guess = fr[np.argmax(S[1:]) + 1]
    p0 = [s.std() * np.sqrt(2), f_guess, Q_guess, 0.0, s.mean()]
    # 相位粗估
    ph = np.linspace(0, 2 * np.pi, 64)
    err = [np.sum((s - damped_sine(t, p0[0], f_guess, Q_guess, p, p0[4])) ** 2) for p in ph]
    p0[3] = ph[int(np.argmin(err))]
    try:
        popt, pcov = curve_fit(damped_sine, t, s, p0=p0, maxfev=40000)
    except Exception as e:
        return {'ok': False, 'err': str(e)}
    perr = np.sqrt(np.diag(pcov))
    resid = s - damped_sine(t, *popt)
    r2 = 1 - np.sum(resid ** 2) / np.sum((s - s.mean()) ** 2)
    return {'ok': True, 'f0': popt[1], 'Q': popt[2], 'A': popt[0], 'phi': popt[3],
            'C': popt[4], 'f0_err': perr[1], 'Q_err': perr[2], 'r2': r2,
            'tau_ms': popt[2] / (np.pi * popt[1]) * 1e3}


# ============================================================
# 2. Allan 方差：找最优积分时间
# ============================================================
def allan_deviation(y, tau_list, rate=1.0):
    """
    y: 相对频差序列 (等间隔)
    tau_list: 要计算的积分时间 (单位: 采样点数 或 秒, 取决于 rate)
    rate: 采样率 (Hz)，若 tau 以秒为单位则传入
    返回: (tau_arr, sigma_arr)
    """
    y = np.asarray(y, float); y = y[~np.isnan(y)]
    out = []
    for tau in tau_list:
        n = int(round(tau * rate))
        if n < 2 or n > len(y) // 2:
            out.append(np.nan); continue
        m = len(y) // n
        if m < 2:
            out.append(np.nan); continue
        seg = y[:m * n].reshape(m, n)
        mu = seg.mean(axis=1)
        out.append(np.sqrt(0.5 * np.mean(np.diff(mu) ** 2)))
    return np.array(tau_list, float), np.array(out, float)

def optimal_tau(tau, sigma):
    """返回最优积分时间与对应 sigma"""
    mask = ~np.isnan(sigma)
    if not mask.any():
        return None, None
    i = np.nanargmin(sigma)
    return tau[mask][i], sigma[mask][i]


# ============================================================
# 3. 双球交换法：核心差分
# ============================================================
def swap_method(delta1, delta2):
    """
    delta1: 配置 I  (球1@A, 球2@B) 的差频数组  = (F1-F2) + (PA-PB)
    delta2: 配置 II (球1@B, 球2@A) 的差频数组  = (F1-F2) - (PA-PB)
    返回: (位置效应估计, 标准误, 95%CI)
    """
    d1 = np.asarray(delta1, float); d2 = np.asarray(delta2, float)
    per_trial = (d1 - d2) / 2.0            # 每轮的位置效应估计
    est = per_trial.mean()
    se = per_trial.std(ddof=1) / np.sqrt(len(per_trial))
    return est, se, (est - 1.96 * se, est + 1.96 * se)


# ============================================================
# 4. 多模温度补偿
# ============================================================
def temp_compensate(f1, f2, T):
    """
    用模态1、模态2 构造对温度一阶不敏感的组合。
    先用 T 拟合各自的温度系数，再解出 k 使 dF/dT = 0。
    返回: (k, alpha1, alpha2, F)
    """
    f1 = np.asarray(f1, float); f2 = np.asarray(f2, float); T = np.asarray(T, float)
    A = np.vstack([T, np.ones_like(T)]).T
    b1, a1 = np.linalg.lstsq(A, f1, rcond=None)[0]
    b2, a2 = np.linalg.lstsq(A, f2, rcond=None)[0]
    # F = f1 - k*f2  ->  dF/dT = b1 - k*b2 = 0
    k = b1 / b2 if abs(b2) > 1e-30 else 0.0
    F = f1 - k * f2
    return k, b1, b2, F


# ============================================================
# 5. 拍频/灵敏度换算
# ============================================================
def beat_sensitivity(T_int, f0):
    """拍频法可分辨的最小相对频差"""
    return 1.0 / (T_int * f0)

def phase_sensitivity(Q, f0, dphi=1e-6):
    """相位斜率法: df = dphi * f0 / (2Q)"""
    df = dphi * f0 / (2 * Q)
    return df, df / f0

def amp_peak_sensitivity(Q, SNR):
    """幅度峰位法 (基线): df/f ~ 1/(Q*SNR)"""
    return 1.0 / (Q * SNR)


# ============================================================
# 6. 物理常量查询
# ============================================================
MU0 = 4 * np.pi * 1e-7
C_LIGHT = 2.99792458e8
MAT = {
    '铜':    dict(sigma=5.96e7, rho=8960., E=117e9, nu=0.34),
    '铝':    dict(sigma=3.77e7, rho=2700., E=69e9,  nu=0.33),
    '不锈钢': dict(sigma=1.4e6,  rho=7900., E=193e9, nu=0.29),
}

def skin_depth(f, mat='铜'):
    """皮肤深度 (m)"""
    return np.sqrt(2.0 / (2 * np.pi * f * MU0 * MAT[mat]['sigma']))

def plate_speed(mat='铜'):
    p = MAT[mat]
    return np.sqrt(p['E'] / (p['rho'] * (1 - p['nu'] ** 2)))

def ring_freq(R, mat='铜'):
    return plate_speed(mat) / (2 * np.pi * R)

def bend_freq(R, t, C=1.0, mat='铜'):
    """薄壳弯曲模量级估计 (C 为模态常数, 1~3)"""
    return ring_freq(R, mat) * (t / R) * C

def cavity_TE011(R):
    """球形腔 TE011 频率 (Hz)"""
    return 2.7437 * C_LIGHT / (2 * np.pi * R)

def q_from_ringdown(f0, tau):
    return np.pi * f0 * tau

def tau_from_q(f0, Q):
    return Q / (np.pi * f0)


# ============================================================
# 7. 八点立方体阵列工具
# ============================================================
def cube_vertices():
    """内接立方体的 8 个球面顶点"""
    v = np.array([[sx, sy, sz] for sx in (1, -1) for sy in (1, -1) for sz in (1, -1)], float)
    return v / np.sqrt(3)

def split_tetrahedra(v):
    """按坐标符号乘积分成两个对映正四面体"""
    prod = np.prod(np.sign(v), axis=1)
    return v[prod > 0], v[prod < 0]

def symmetry_adapted(s):
    """
    s: 4 路信号 [s1,s2,s3,s4]
    返回: S(A1), D1, D2, D3 (T2)
    """
    s = np.asarray(s, float)
    return np.array([
        0.5 * ( s[0] + s[1] + s[2] + s[3]),
        0.5 * ( s[0] + s[1] - s[2] - s[3]),
        0.5 * ( s[0] - s[1] + s[2] - s[3]),
        0.5 * ( s[0] - s[1] - s[2] + s[3]),
    ])

def parity_project(A, B):
    """对映点宇称分离: 和 -> 偶 l, 差 -> 奇 l"""
    A = np.asarray(A, float); B = np.asarray(B, float)
    return A + B, A - B


# ============================================================
# 自检 / 演示
# ============================================================
def demo():
    print("=" * 62)
    print("铜球实验分析工具包 · 自检")
    print("=" * 62)

    # 1) ring-down 拟合
    print("\n[1] ring-down 拟合")
    fs = 48000.; t = np.arange(0, 1.0, 1 / fs)
    f_true, Q_true = 396.0, 300.0
    sig = damped_sine(t, 1.0, f_true, Q_true, 0.7, 0.02)
    rng = np.random.default_rng(0)
    sig += rng.normal(0, 0.02, len(sig))
    r = fit_ringdown(t, sig)
    print(f"    真值 f0={f_true:.2f} Hz, Q={Q_true}")
    print(f"    拟合 f0={r['f0']:.3f}±{r['f0_err']:.3f} Hz, Q={r['Q']:.1f}±{r['Q_err']:.1f}, R²={r['r2']:.4f}")
    print(f"    ring-down tau = {r['tau_ms']:.1f} ms")

    # 2) 双球交换
    print("\n[2] 双球交换法")
    rng = np.random.default_rng(1)
    F1mF2 = 137.5                      # 球的本征差 (Hz)，会被抵消
    PAmPB = 0.004                      # 待求的位置效应 (Hz)
    d1 = F1mF2 + PAmPB + rng.normal(0, .002, 40)
    d2 = F1mF2 - PAmPB + rng.normal(0, .002, 40)
    est, se, ci = swap_method(d1, d2)
    print(f"    真值 P_A-P_B = {PAmPB:.4f} Hz")
    print(f"    估计 = {est:.4f} ± {se:.4f} Hz, 95%CI = [{ci[0]:.4f}, {ci[1]:.4f}]")

    # 3) 温度补偿
    print("\n[3] 多模温度补偿")
    T = np.linspace(20, 25, 300) + rng.normal(0, .01, 300)
    f1 = 396.0 * (1 - 1.5e-4 * (T - 22)) + rng.normal(0, 1e-4, 300)
    f2 = 612.0 * (1 - 1.2e-4 * (T - 22)) + rng.normal(0, 1e-4, 300)
    k, b1, b2, F = temp_compensate(f1, f2, T)
    print(f"    alpha1={b1:.5f} Hz/K, alpha2={b2:.5f} Hz/K -> k={k:.4f}")
    print(f"    补偿后 F 对 T 的残余斜率 = {np.polyfit(T, F, 1)[0]:.2e} Hz/K (应≈0)")

    # 4) 灵敏度
    print("\n[4] 三种读数方式灵敏度 (f0=2.618GHz, Q=1e4, SNR=1e3)")
    print(f"    幅度峰位法 : {amp_peak_sensitivity(1e4, 1e3):.2e}")
    df, rel = phase_sensitivity(1e4, 2.618e9)
    print(f"    相位斜率法 : {rel:.2e}   (df={df:.3f} Hz)")
    for T_ in [1, 100, 1000]:
        print(f"    拍频 T={T_:<5}s : {beat_sensitivity(T_, 2.618e9):.2e}")

    # 5) 物理常量
    print("\n[5] 物理常量速查")
    print(f"    铜 @396Hz 皮肤深度 = {skin_depth(396)*1e3:.2f} mm  (壁厚 0.5mm -> 电磁透明)")
    print(f"    R=5cm 环频率 = {ring_freq(0.05):.0f} Hz")
    print(f"    R=5cm t=0.5mm 弯曲模 = {bend_freq(0.05, 5e-4):.1f} Hz")
    print(f"    R=5cm TE011 腔模 = {cavity_TE011(0.05)/1e9:.3f} GHz")
    print(f"    396Hz Q=300 -> tau = {tau_from_q(396,300)*1e3:.1f} ms")

    # 6) 阵列
    print("\n[6] 八点阵列")
    v = cube_vertices(); A, B = split_tetrahedra(v)
    print(f"    组A 4点, 组B 4点, 对映关系: B = -A -> {np.allclose(B, -A[np.argsort([np.argmin(np.linalg.norm(a+b) for b in B) for a in A])])}")
    s = np.array([1.0, 0.8, 0.6, 0.4])
    print(f"    4路信号 {s} -> 对称性适配 {np.round(symmetry_adapted(s),4)}")
    print("\n自检完成。")


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--demo', action='store_true', help='运行自检演示')
    args = ap.parse_args()
    if args.demo:
        demo()
    else:
        print(__doc__)
        print("用法: python 铜球实验分析工具.py --demo")

``````