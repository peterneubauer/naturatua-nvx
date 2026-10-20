# SWEREF 99 TM <-> WGS84 (GRS80), Krügers formler enligt Lantmäteriet. Noggrannhet ~1 mm.
import math
a = 6378137.0; f = 1 / 298.257222101; e2 = f * (2 - f); n = f / (2 - f)
k0, L0, FE = 0.9996, math.radians(15.0), 500000.0
Ah = a / (1 + n) * (1 + n**2 / 4 + n**4 / 64)
B = (n / 2 - 2 * n**2 / 3 + 5 * n**3 / 16 + 41 * n**4 / 180, 13 * n**2 / 48 - 3 * n**3 / 5 + 557 * n**4 / 1440,
     61 * n**3 / 240 - 103 * n**4 / 140, 49561 * n**4 / 161280)
D = (n / 2 - 2 * n**2 / 3 + 37 * n**3 / 96 - n**4 / 360, n**2 / 48 + n**3 / 15 - 437 * n**4 / 1440,
     17 * n**3 / 480 - 37 * n**4 / 840, 4397 * n**4 / 161280)

def to_tm(lon, lat):
    phi, lam = math.radians(lat), math.radians(lon) - L0
    A_, B_, C_, D_ = e2, (5 * e2**2 - e2**3) / 6, (104 * e2**3 - 45 * e2**4) / 120, 1237 * e2**4 / 1260
    s = math.sin(phi)
    ps = phi - s * math.cos(phi) * (A_ + B_ * s**2 + C_ * s**4 + D_ * s**6)
    xi, eta = math.atan(math.tan(ps) / math.cos(lam)), math.atanh(math.cos(ps) * math.sin(lam))
    N = k0 * Ah * (xi + sum(b * math.sin(2 * j * xi) * math.cosh(2 * j * eta) for j, b in enumerate(B, 1)))
    E = k0 * Ah * (eta + sum(b * math.cos(2 * j * xi) * math.sinh(2 * j * eta) for j, b in enumerate(B, 1))) + FE
    return E, N

def to_wgs84(E, N):
    xi, eta = N / (k0 * Ah), (E - FE) / (k0 * Ah)
    xp = xi - sum(d * math.sin(2 * j * xi) * math.cosh(2 * j * eta) for j, d in enumerate(D, 1))
    ep = eta - sum(d * math.cos(2 * j * xi) * math.sinh(2 * j * eta) for j, d in enumerate(D, 1))
    ps = math.asin(math.sin(xp) / math.cosh(ep)); lam = math.atan(math.sinh(ep) / math.cos(xp))
    As = e2 + e2**2 + e2**3 + e2**4
    Bs = -(7 * e2**2 + 17 * e2**3 + 30 * e2**4) / 6
    Cs = (224 * e2**3 + 889 * e2**4) / 120
    Ds = -(4279 * e2**4) / 1260
    s = math.sin(ps)
    phi = ps + s * math.cos(ps) * (As + Bs * s**2 + Cs * s**4 + Ds * s**6)
    return math.degrees(L0 + lam), math.degrees(phi)
