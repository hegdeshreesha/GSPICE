"""Clean-room BSIM4.8.3 Id(Vg) reference calculator (faithful transcription).

Every equation below is transcribed from the official BSIM4.8.3 CMC source
(ECL-2.0, numeric reference only).  Source locations:
  - phi/ni/Eg0/Xdep0/cdep0/vbi/theta0vb0/thetaRout : b4temp.c:211-213,1323-1375,1519-1542
  - k1ox/k2ox/vfb/vtfbphi1/vtfbphi2/coxe/coxp/toxp : b4temp.c:1488,1774-1804,2403
  - Vbseff/sqrtPhis/Xdep/Vth/n/Vgsteff/mobility    : b4ld.c:1002-1025,1034-1154,1238-1296
  - Vdsat/Vdseff/Vasat/Idl/Cclm/Id                 : b4ld.c:1313-1336,1592-2091
Constants (SI):
  MIN_EXP=1e-300, EXP_THRESHOLD=80, MAX_EXP=1e36 (b4const.h)
"""
import math
import sys

Kb = 1.380649e-23
q = 1.602176634e-19
EPS0 = 8.8541878128e-12
epsr_ox = 3.9
epsr_sub = 11.7
epssub = epsr_sub * EPS0

MIN_EXP = 1.0e-300
MAX_EXP = 1.0e36
EXP_THRESHOLD = 80.0

def fermi_ni(Tnom):
    Vtm0 = (Kb / q) * Tnom
    Eg0 = 1.16 - 7.02e-4 * Tnom * Tnom / (Tnom + 1108.0)
    ni = 1.45e10 * (Tnom / 300.15) * math.sqrt(Tnom / 300.15) * \
        math.exp(21.5565981 - Eg0 / (2.0 * Vtm0))
    return Vtm0, Eg0, ni

def prep(p):
    Tnom = p.get("T", 300.15)
    Vtm0, Eg0, ni = fermi_ni(Tnom)
    ndep = p.get("ndep", 1.7e17)
    nsub = p.get("nsub", 6.0e16)
    phin = p.get("phin", 0.0)
    phi = Vtm0 * math.log(ndep / ni) + phin + 0.4
    sqrtPhi = math.sqrt(phi)
    Xdep0 = math.sqrt(2.0 * epssub / (q * ndep * 1e6)) * sqrtPhi
    cdep0 = math.sqrt(q * epssub * ndep * 1e6 / (2.0 * phi))
    nsd = p.get("nsd", 1.0e20)
    vbi = Vtm0 * math.log(nsd * ndep / (ni * ni))
    toxe = p.get("toxe", 3.0e-9)
    toxm = p.get("toxm", toxe)
    toxp = p.get("toxp", toxe)          # b4set: toxp=toxe when toxe given
    k1 = p.get("k1", 0.0)
    k2 = p.get("k2", 0.0)
    k1ox = k1 * toxe / toxm
    k2ox = k2 * toxe / toxm
    coxe = epsr_ox * EPS0 / toxe        # b4temp.c:178
    coxp = epsr_ox * EPS0 / toxp        # b4temp.c:180
    vth0 = p.get("vth0", 0.4)
    mtype = p.get("type", 1)            # +1 nmos, -1 pmos
    # vfb / vtfbphi1 / vtfbphi2 (b4temp.c:1488,1778-1788), ngate=0 & mtrlMod=0
    vfb = mtype * vth0 - phi - k1 * sqrtPhi
    T3 = mtype * vth0 - vfb - phi
    vtfbphi1 = max(0.0, 2.0 * T3 if mtype == 1 else 2.5 * T3)
    vtfbphi2 = max(0.0, 4.0 * T3)
    # theta0vb0 (b4temp.c:1519-1529)
    leff = p.get("leff", 1.0e-6)
    dsub = p.get("dsub", 0.56)
    tmp = math.sqrt(epssub / (epsr_ox * EPS0) * toxe * Xdep0)
    T0 = dsub * leff / tmp
    if T0 < EXP_THRESHOLD:
        T1 = math.exp(T0); T2 = T1 - 1.0
        theta0vb0 = T1 / (T2 * T2 + 2.0 * T1 * MIN_EXP)
    else:
        theta0vb0 = 1.0 / (MAX_EXP - 2.0)
    # thetaRout (b4temp.c:1531-1542)
    drout = p.get("drout", 0.56)
    T0 = drout * leff / tmp
    if T0 < EXP_THRESHOLD:
        T1 = math.exp(T0); T2 = T1 - 1.0
        thetaRout = p.get("pdibl1", 0.39) * (T1 / (T2 * T2 + 2.0 * T1 * MIN_EXP)) \
            + p.get("pdibl2", 0.0086)
    else:
        thetaRout = p.get("pdibl1", 0.39) * (1.0 / (MAX_EXP - 2.0)) + p.get("pdibl2", 0.0086)
    # litl (b4temp.c:1347), mtrlMod=0
    xj = p.get("xj", 1.5e-7)
    litl = math.sqrt(3.0 * 3.9 / epsr_ox * xj * toxe)
    factor1 = math.sqrt(epssub / (epsr_ox * EPS0) * toxe)  # b4temp.c:205
    return dict(locals())

def evaluate(prep, vg, vd, vs=0.0, vb=0.0, T=None):
    p = prep["p"]
    mtype = prep["mtype"]
    vth0 = prep["vth0"]
    vtm0 = prep["Vtm0"]
    weff = p.get("weff", 1.0e-6)
    leff = p.get("leff", 1.0e-6)
    toxe = prep["toxe"]
    phi = prep["phi"]
    sqrtPhi = prep["sqrtPhi"]
    Xdep0 = prep["Xdep0"]
    cdep0 = prep["cdep0"]
    vbi = prep["vbi"]
    k1ox = prep["k1ox"]
    k2ox = prep["k2ox"]
    coxe = prep["coxe"]
    coxp = prep["coxp"]
    theta0vb0 = prep["theta0vb0"]
    thetaRout = prep["thetaRout"]
    litl = prep["litl"]
    factor1 = prep["factor1"]
    vtfbphi2 = prep["vtfbphi2"]
    toxp = prep["toxp"]
    T = p.get("T", 300.15) if T is None else T
    Vtm = (Kb / q) * T

    vgs = vg - vs
    vds = vd - vs
    vsb = max(0.0, vs - vb)

    # --- Vbseff, Phis, Xdep (b4ld.c:1002-1027) ---
    vbsc = p.get("vbsc", -30.0)
    T0 = vsb - vbsc - 0.001
    T1 = math.sqrt(T0 * T0 - 0.004 * vbsc)
    if T0 >= 0.0:
        Vbseff = vbsc + 0.5 * (T0 + T1)
    else:
        T2 = -0.002 / (T1 - T0)
        Vbseff = vbsc * (1.0 + T2)
    T9 = 0.95 * phi
    T0 = T9 - Vbseff - 0.001
    T1 = math.sqrt(T0 * T0 + 0.004 * T9)
    Vbseff = T9 - 0.5 * (T0 + T1)
    Phis = phi - Vbseff
    sqrtPhis = math.sqrt(Phis)
    Xdep = Xdep0 * sqrtPhis / sqrtPhi

    # --- Vth (b4ld.c:1034-1129) ---
    V0 = vbi - phi
    dvt0 = p.get("dvt0", 2.2); dvt1 = p.get("dvt1", 0.53); dvt2 = p.get("dvt2", -0.032)
    dvt0w = p.get("dvt0w", 0.0); dvt1w = p.get("dvt1w", 5.3e6); dvt2w = p.get("dvt2w", -0.032)
    def theta_of(dvt1_, dvt2_):
        T0 = dvt2_ * Vbseff
        if T0 >= -0.5:
            T1 = 1.0 + T0
        else:
            T4 = 1.0 / (3.0 + 8.0 * T0)
            T1 = (1.0 + 3.0 * T0) * T4
        lt = factor1 * math.sqrt(Xdep) * T1
        T0 = dvt1_ * leff / lt
        if T0 < EXP_THRESHOLD:
            T1 = math.exp(T0); T2 = T1 - 1.0
            return T1 / (T2 * T2 + 2.0 * T1 * MIN_EXP)
        return 1.0 / (MAX_EXP - 2.0)
    th0 = theta_of(dvt1, dvt2)
    th0w = theta_of(dvt1w, dvt2w)
    Delt_vth = dvt0 * th0 * V0
    T2w = dvt0w * th0w * V0

    TempRatio = T / p.get("Tnom", T) - 1.0
    lpe0 = p.get("lpe0", 1.74e-7)
    lpeb = p.get("lpeb", 0.0)
    Lpe_Vb = math.sqrt(1.0 + lpeb / leff)
    Lpe_T = k1ox * (math.sqrt(1.0 + lpe0 / leff) - 1.0) * sqrtPhi \
        + (p.get("kt1", -0.11) + p.get("kt1l", 0.0) / leff
           + p.get("kt2", 0.022) * Vbseff) * TempRatio

    w0 = p.get("w0", 2.5e-6)
    Vth_NarrowW = toxe * phi / (weff + w0)

    eta0 = p.get("eta0", 0.08); etab = p.get("etab", -0.07)
    T3 = eta0 + etab * Vbseff
    if T3 < 1.0e-4:
        T9 = 1.0 / (3.0 - 2.0e4 * T3)
        T3 = (2.0e-4 - T3) * T9
    DIBL_Sft = T3 * theta0vb0 * vds

    k3 = p.get("k3", 80.0); k3b = p.get("k3b", 0.0)
    Vth = mtype * vth0 + (k1ox * sqrtPhis - prep["k1"] * sqrtPhi) * Lpe_Vb \
        - k2ox * Vbseff - Delt_vth - T2w \
        + (k3 + k3b * Vbseff) * Vth_NarrowW + Lpe_T - DIBL_Sft

    # --- n (b4ld.c:1133-1154) ---
    nfactor = p.get("nfactor", 1.0)
    cdsc = p.get("cdsc", 2.4e-4); cdscb = p.get("cdscb", 0.0); cdscd = p.get("cdscd", 0.0)
    cit = p.get("cit", 0.0)
    tmp1 = epssub / Xdep
    tmp2 = nfactor * tmp1
    tmp3 = cdsc + cdscb * Vbseff + cdscd * vds
    tmp4 = (tmp2 + tmp3 * th0 + cit) / coxe
    if tmp4 >= -0.5:
        n = 1.0 + tmp4
    else:
        T0 = 1.0 / (3.0 + 8.0 * tmp4)
        n = (1.0 + 3.0 * tmp4) * T0

    # --- Vgsteff smoothing (b4ld.c:1236-1296) ---
    mstar = 0.5 + math.atan(p.get("minv", 0.0)) / math.pi   # b4temp.c:1425
    voffcbn = p.get("voff", -0.08) + p.get("voffl", 0.0) / leff
    Vgst = vgs - Vth
    T0 = n * Vtm
    T1 = mstar * Vgst
    T2 = T1 / T0
    if T2 > EXP_THRESHOLD:
        T10 = T1
    elif T2 < -EXP_THRESHOLD:
        T10 = n * Vtm * math.log(1.0 + MIN_EXP)
    else:
        T10 = n * Vtm * math.log(1.0 + math.exp(T2))
    T1 = voffcbn - (1.0 - mstar) * Vgst
    T2 = T1 / T0
    if T2 < -EXP_THRESHOLD:
        T3 = coxe * MIN_EXP / cdep0
        T9 = mstar + T3 * n
    elif T2 > EXP_THRESHOLD:
        T3 = coxe * MAX_EXP / cdep0
        T9 = mstar + T3 * n
    else:
        T3 = coxe / cdep0
        T9 = mstar + n * T3 * math.exp(T2)
    Vgsteff = T10 / T9

    # --- Mobility mobMod=0 (b4ld.c:1422-1439,1561-1578) ---
    u0 = p.get("u0", 0.067) * p.get("u0mult", 1.0)
    Tnom = p.get("Tnom", T)
    u0temp = u0 * math.pow(T / Tnom, p.get("ute", -1.5))
    deltaTemp = T - Tnom
    ua = p.get("ua", 1.0e-9) * (1.0 + p.get("ua1", 1.0e-9) * deltaTemp)
    ub = p.get("ub", 1.0e-19) * (1.0 + p.get("ub1", -1.0e-18) * deltaTemp)
    uc = p.get("uc", -0.0465e-9) * (1.0 + p.get("uc1", -0.056e-9) * deltaTemp)
    ud = p.get("ud", 0.0) * (1.0 + p.get("ud1", 0.0) * deltaTemp)
    T0 = Vgsteff + 2.0 * Vth
    T2 = ua + uc * Vbseff
    T3 = T0 / toxe
    T12 = math.sqrt(Vth * Vth + 0.0001)
    T9 = 1.0 / (Vgsteff + 2.0 * T12)
    T10 = T9 * toxe
    T8 = ud * T10 * T10 * Vth
    T6 = T8 * Vth
    T5 = T3 * (T2 + ub * T3) + T6
    if T5 >= -0.8:
        Denomi = 1.0 + T5
    else:
        T9 = 1.0 / (7.0 + 10.0 * T5)
        Denomi = (0.6 + T5) * T9
    ueff = u0temp / Denomi

    # --- Esat / WVCox / Rds (b4ld.c:1313-1336,1580-1589) ---
    vsattemp = p.get("vsat", 1.0e5) * (1.0 - p.get("at", 3.3e4) * deltaTemp)
    Esat = 2.0 * vsattemp / ueff
    EsatL = Esat * leff
    WVCox = weff * vsattemp * coxe
    rdsmod = p.get("rdsmod", 0)
    if rdsmod == 1:
        Rds = 0.0
    else:
        prwg = p.get("prwg", 1.0); prwb = p.get("prwb", 0.0)
        T0 = 1.0 + prwg * Vgsteff
        T1 = prwb * sqrtPhis
        T2 = 1.0 / T0 + T1
        T3 = T2 + math.sqrt(T2 * T2 + 0.01)
        rds0 = p.get("rdsw", 200.0) * math.pow(weff * 1e6, -p.get("wr", 1.0))
        Rds = p.get("rdswmin", 0.0) + T3 * (0.5 * rds0)
    WVCoxRds = WVCox * Rds

    # --- Lambda (b4ld.c:1592-1609); a1=0 -> a2 ---
    a1 = p.get("a1", 0.0); a2 = p.get("a2", 1.0)
    if a1 == 0.0:
        Lambda = a2
    else:
        raise NotImplementedError("a1 != 0 not needed")

    # --- Abulk (b4ld.c:1338-1376) ---
    xj = prep["xj"]
    T9 = 0.5 * k1ox * Lpe_Vb / sqrtPhis
    T1a = T9 + k2ox - k3b * Vth_NarrowW
    T9 = math.sqrt(xj * Xdep)
    tmp1 = leff + 2.0 * T9
    T5 = leff / tmp1
    tmp2 = p.get("a0", 1.0) * T5
    tmp3 = weff + p.get("b1", 0.0)
    tmp4 = p.get("b0", 0.0) / tmp3
    T2 = tmp2 + tmp4
    Abulk0 = 1.0 + T1a * T2
    T8 = p.get("ags", 0.0) * p.get("a0", 1.0) * T5 * T5 * T5
    Abulk = Abulk0 - T8 * Vgsteff
    if Abulk0 < 0.1:
        T9 = 1.0 / (3.0 - 20.0 * Abulk0)
        Abulk0 = (0.2 - Abulk0) * T9
    if Abulk < 0.1:
        T9 = 1.0 / (3.0 - 20.0 * Abulk)
        Abulk = (0.2 - Abulk) * T9

    # --- Vdsat (b4ld.c:1611-1679) ---
    Vgst2Vtm = Vgsteff + 2.0 * Vtm
    if (Rds == 0.0) and (Lambda == 1.0):
        T0 = 1.0 / (Abulk * EsatL + Vgst2Vtm)
        T2 = Vgst2Vtm * T0
        T3 = EsatL * Vgst2Vtm
        Vdsat = T3 * T0
    else:
        T9 = Abulk * WVCoxRds
        T8 = Abulk * T9
        T7 = Vgst2Vtm * T9
        T6 = Vgst2Vtm * WVCoxRds
        T0 = 2.0 * Abulk * (T9 - 1.0 + 1.0 / Lambda)
        T1 = Vgst2Vtm * (2.0 / Lambda - 1.0) + Abulk * EsatL + 3.0 * T7
        T2 = Vgst2Vtm * (EsatL + 2.0 * T6)
        T3 = math.sqrt(T1 * T1 - 2.0 * T0 * T2)
        Vdsat = (T1 - T3) / T0

    # --- Vdseff (b4ld.c:1682-1721) ---
    delta = p.get("delta", 0.01)
    T1 = Vdsat - vds - delta
    T2 = math.sqrt(T1 * T1 + 4.0 * delta * Vdsat)
    if T1 >= 0.0:
        Vdseff = Vdsat - 0.5 * (T1 + T2)
    else:
        T4 = 2.0 * delta / (T2 - T1)
        T5 = 1.0 - T4
        Vdseff = Vdsat * T5
    if Vdseff > vds:
        Vdseff = vds
    diffVds = vds - Vdseff

    # --- Vasat (b4ld.c:1765-1788) ---
    tmp4 = 1.0 - 0.5 * Abulk * Vdsat / Vgst2Vtm
    T9 = WVCoxRds * Vgsteff
    T8 = T9 / Vgst2Vtm
    T0 = EsatL + Vdsat + 2.0 * T9 * tmp4
    T9 = WVCoxRds * Abulk
    T1 = 2.0 / Lambda - 1.0 + T9
    Vasat = T0 / T1

    # --- Idl via Coxeff quantum correction (b4ld.c:1790-1849) ---
    T0 = (Vgsteff + vtfbphi2) / (2.0e8 * toxp)
    tmp3 = math.exp(p.get("bdos", 1.0) * 0.7 * math.log(T0))
    T1 = 1.0 + tmp3
    Tcen = p.get("ados", 1.0) * 1.9e-9 / T1
    Coxeff = epssub * coxp / (epssub + coxp * Tcen)
    beta = ueff * Coxeff * weff / leff
    T0 = 1.0 - 0.5 * Vdseff * Abulk / Vgst2Vtm
    fgche1 = Vgsteff * T0
    T9 = Vdseff / EsatL
    fgche2 = 1.0 + T9
    gche = beta * fgche1 / fgche2
    T0 = 1.0 + gche * Rds
    Idl = gche / T0

    # --- CLM / DIBL / DITS / SCBE (b4ld.c:1853-2091) ---
    fprout = p.get("fprout", 0.0)
    if fprout <= 0.0:
        FP = 1.0
    else:
        FP = 1.0 / (1.0 + fprout * math.sqrt(leff) / Vgst2Vtm)
    pvag = p.get("pvag", 0.0)
    T9 = pvag / EsatL * Vgsteff
    if T9 > -0.9:
        PvagTerm = 1.0 + T9
    else:
        T4 = 1.0 / (17.0 + 20.0 * T9)
        PvagTerm = (0.8 + T9) * T4
    pclm = p.get("pclm", 1.3)
    if (pclm > MIN_EXP) and (diffVds > 1.0e-10):
        T0 = 1.0 + Rds * Idl
        T2 = Vdsat / Esat
        T1 = leff + T2
        Cclm = FP * PvagTerm * T0 * T1 / (pclm * litl)
        VACLM = Cclm * diffVds
    else:
        VACLM = MAX_EXP
    if thetaRout > MIN_EXP:
        T8 = Abulk * Vdsat
        T0 = Vgst2Vtm * T8
        T1 = Vgst2Vtm + T8
        T2 = thetaRout
        VADIBL = (Vgst2Vtm - T0 / T1) / T2
        T7 = p.get("pdiblb", 0.0) * Vbseff
        if T7 >= -0.9:
            T3 = 1.0 / (1.0 + T7)
            VADIBL *= T3
        else:
            T4 = 1.0 / (0.8 + T7)
            VADIBL *= (17.0 + 20.0 * T7) * T4
        VADIBL *= PvagTerm
    else:
        VADIBL = MAX_EXP
    Va = Vasat + VACLM
    pdits = p.get("pdits", 0.0)
    if pdits > MIN_EXP:
        T0 = p.get("pditsd", 0.0) * vds
        if T0 > EXP_THRESHOLD:
            T1 = MAX_EXP
        else:
            T1 = math.exp(T0)
        T2 = 1.0 + p.get("pditsl", 0.0) * leff
        VADITS = (1.0 + T2 * T1) / pdits * FP
    else:
        VADITS = MAX_EXP
    pscbe1 = p.get("pscbe1", 4.24e8); pscbe2 = p.get("pscbe2", 1.0e-5)
    if (pscbe2 > 0.0) and (pscbe1 >= 0.0):
        if diffVds > pscbe1 * litl / EXP_THRESHOLD:
            T0 = pscbe1 * litl / diffVds
            VASCBE = leff * math.exp(T0) / pscbe2
        else:
            VASCBE = MAX_EXP * leff / pscbe2
    else:
        VASCBE = MAX_EXP

    Idsa = Idl * (1.0 + diffVds / VADIBL) * (1.0 + diffVds / VADITS)
    T0 = math.log(Va / Vasat)
    T1 = T0 / Cclm
    Idsa *= (1.0 + T1)
    Ids = Idsa * (1.0 + diffVds / VASCBE) * Vdseff  # cdrain = Ids*Vdseff

    return dict(Vth=Vth, Vgsteff=Vgsteff, Vdsat=Vdsat, n=n, ueff=ueff,
                Abulk=Abulk, Rds=Rds, Id=Ids, Idl=Idl, phi=phi, Xdep0=Xdep0,
                cdep0=cdep0, Vth_NarrowW=Vth_NarrowW, DIBL_Sft=DIBL_Sft,
                Vbseff=Vbseff, Vdseff=Vdseff, VACLM=VACLM, Vasat=Vasat,
                Coxeff=Coxeff, beta=beta)

if __name__ == "__main__":
    # Deck registry: params (SI) + grid points (vg, vd, vs, vb, T).
    # Mirrors the deck registry in tests/bsim4_reference_probe.cpp exactly.
    deck_a_params = dict(toxe=10e-9, u0=0.05, vth0=0.4, k1=0, k2=0, dvt0=0,
                         dvt1=0.53, dvt2=-0.032, xj=1.5e-7, ua=2.0e-9,
                         ub=5.0e-19, vsat=1.0e5, T=300.15, Tnom=300.15)
    deck_b_params = dict(toxe=5e-9, toxm=5e-9, u0=0.028, vth0=0.50,
                         ndep=2.0e17, nsd=1.0e20, k1=0.6, k2=-0.05, k3=120.0,
                         k3b=0.5, w0=0.6e-6, lpe0=1.2e-7, lpeb=0.3, dvt0=1.8,
                         dvt1=0.42, dvt2=-0.02, dvt0w=0.4, dvt1w=1.2e6,
                         dvt2w=-0.03, dsub=0.35, drout=0.45, nfactor=1.2,
                         cdsc=1.2e-4, cdscb=0.02, cdscd=0.005, cit=1.0e-5,
                         eta0=0.15, etab=-0.08, xj=1.2e-7, vsat=1.1e5,
                         a0=0.9, ags=0.2, b0=5e-8, b1=1e-8, keta=0.03,
                         a1=0.0, a2=1.0, pclm=1.2, pdibl1=0.2, pdibl2=0.003,
                         pvag=0.1, pdiblb=0.02, fprout=0.05, rdsmod=0,
                         rdsw=150.0, wr=1.0, prwg=1.0, prwb=0.1,
                         delta=0.01, voff=-0.08, minv=0.2, ua=2.0e-9,
                         ub=6.0e-19, T=300.15, Tnom=300.15)
    deck_c_params = dict(toxe=10e-9, u0=0.05, vth0=0.4, k1=0, k2=0, dvt0=0,
                         dvt1=0.53, dvt2=-0.032, xj=1.5e-7, ua=2.0e-9,
                         ub=5.0e-19, vsat=1.0e5, weff=5.0e-6, leff=0.3e-6,
                         T=300.15, Tnom=300.15)
    deck_d_params = dict(toxe=10e-9, u0=0.05, vth0=0.4, k1=0, k2=0, dvt0=0,
                         dvt1=0.53, dvt2=-0.032, xj=1.5e-7, ua=2.0e-9,
                         ub=5.0e-19, vsat=1.0e5, at=3.3e-4, Tnom=300.15)
    grids = {
        "a": [(0.1 * i, 0.8, 0.0, 0.0, 300.15) for i in range(11)],
        "b": [(vg, 0.6, 0.0, vb, 300.15)
              for vg in (0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2)
              for vb in (0.0, -0.3, 0.3)],
        "c": [(0.8, vd, 0.0, 0.0, 300.15)
              for vd in (0.02, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5)],
        "d": [(vg, 0.8, 0.0, 0.0, T)
              for vg in (0.2, 0.5, 0.8, 1.1) for T in (273.15, 373.15)],
    }
    decks = {"a": deck_a_params, "b": deck_b_params,
             "c": deck_c_params, "d": deck_d_params}

    deck = None
    if "--deck" in sys.argv:
        deck = sys.argv[sys.argv.index("--deck") + 1]
        if deck not in grids:
            print("unknown deck", deck, file=sys.stderr)
            sys.exit(2)
    if "--json" in sys.argv and deck is not None:
        p_ = prep(dict(decks[deck]))
        for vg, vd, vs, vb, T in grids[deck]:
            print("%.12e" % evaluate(p_, vg, vd, vs, vb, T)["Id"])
        sys.exit(0)
    base = decks.get(deck, decks["a"]) if deck else decks["a"]
    p_ = prep(base)
    print("phi=%.6f Xdep0=%.3e cdep0=%.3e vbi=%.4f coxe=%.3e theta0vb0=%.3e"
          % (p_["phi"], p_["Xdep0"], p_["cdep0"], p_["vbi"], p_["coxe"], p_["theta0vb0"]))
    print(" Vg     Vth      Vgsteff    Vdsat      n       Id")
    for i in range(11):
        vg = 0.1 * i
        r = evaluate(p_, vg, 0.8)
        print("%4.1f  %.4f  %9.4e  %9.4e %5.3f  %9.4e" %
              (vg, r["Vth"], r["Vgsteff"], r["Vdsat"], r["n"], r["Id"]))
