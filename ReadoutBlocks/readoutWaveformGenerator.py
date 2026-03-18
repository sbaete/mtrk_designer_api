import numpy as np
import matplotlib.pyplot as plt

## Adapted and extended from the Pulpy library by J.B. Martin (https://github.com/jonbmartin/pulpy/tree/master)


__all__ = [
    "min_trap_grad",
    "trap_grad",
    "spiral_varden",
    "spiral_arch",
    "spiral_k",
    "mtrk_epi",
    "rosette",
    "spokes_grad",
    "stack_of",
    "traj_array_to_complex",
    "traj_complex_to_array",
]


def min_trap_grad(area, gmax, dgdt, dt):
    r"""Minimal duration trapezoidal gradient designer. Design for target area
    under the flat portion (for non-ramp-sampled pulses)

    Args:
        area (float): pulse area in (g*sec)/cm
        gmax (float): maximum gradient in g/cm
        dgdt (float): max slew rate in g/cm/sec
        dt (float): sample time in sec

    Returns:
        2-element tuple containing

        - **trap** (*array*): gradient waveform in g/cm.
        - **ramppts** (*int*): number of points in ramps.

    """

    if np.abs(area) > 0:
        # we get the solution for plateau amp by setting derivative of
        # duration as a function of amplitude to zero and solving
        a = np.sqrt(dgdt * area / 2)

        # finish design with discretization
        # make a flat portion of magnitude a and enough area for the swath
        pts = np.floor(area / a / dt)
        flat = np.ones((1, int(pts)))
        flat = flat / np.sum(flat) * area / dt
        if np.max(flat) > gmax:
            flat = np.ones((1, int(np.ceil(area / gmax / dt))))
            flat = flat / np.sum(flat) * area / dt

        # make attack and decay ramps
        ramppts = int(np.ceil(np.max(flat) / dgdt / dt))
        ramp_up = np.linspace(0, ramppts, num=ramppts + 1) / ramppts * np.max(flat)
        ramp_dn = np.linspace(ramppts, 0, num=ramppts + 1) / ramppts * np.max(flat)

        trap = np.concatenate((ramp_up, np.squeeze(flat), ramp_dn))

    else:
        # negative-area trap requested?
        trap, ramppts = 0, 0

    return np.expand_dims(trap, axis=0), ramppts


def trap_grad(area, gmax, dgdt, dt, *args):
    r"""General trapezoidal gradient designer for total target area
    (for rewinders)

    Args:
        area (float): pulse area in (g*sec)/cm
        gmax (float): maximum gradient in g/cm
        dgdt (float): max slew rate in g/cm/sec
        dt (float): sample time in sec

    Returns:
        2-element tuple containing

        - **trap** (*array*): gradient waveform in g/cm.
        - **ramppts** (*int*): number of points in ramps.

    """

    if len(args) < 5:
        # in case we are making a rewinder
        rampsamp = 1

    if np.abs(area) > 0:
        if rampsamp:
            ramppts = int(np.ceil(gmax / dgdt / dt))
            triareamax = ramppts * dt * gmax
            sign = 1
            if area<0:
                area = -area
                sign = -1
            if triareamax > np.abs(area):
                # triangle pulse
                newgmax = np.sqrt(np.abs(area) * dgdt)
                ramppts = int(np.ceil(newgmax / dgdt / dt))
                ramp_up = np.linspace(0, ramppts, num=ramppts + 1) / ramppts
                ramp_dn = np.linspace(ramppts, 0, num=ramppts + 1) / ramppts
                pulse = sign * np.concatenate((ramp_up, ramp_dn))
            else:
                # trapezoid pulse
                nflat = int(np.ceil((area - triareamax) / gmax / dt / 2) * 2)
                # nflat = np.abs(nflat) # this is a test
                ramp_up = np.linspace(0, ramppts, num=ramppts + 1) / ramppts
                ramp_dn = np.linspace(ramppts, 0, num=ramppts + 1) / ramppts
                pulse = sign * np.concatenate((ramp_up, np.ones(nflat), ramp_dn))

            trap = pulse * (area / (sum(pulse) * dt))

        else:
            # make a flat portion of magnitude gmax
            # and enough area for the entire swath
            flat = np.ones(1, np.ceil(area / gmax / dt))
            flat = flat / sum(flat) * area / dt
            flat_top = np.max(flat)

            # make attack and decay ramps
            ramppts = int(np.ceil(np.max(flat) / dgdt / dt))
            ramp_up = np.linspace(0, ramppts, num=ramppts + 1) / ramppts * flat_top
            ramp_dn = np.linspace(ramppts, 0, num=ramppts + 1) / ramppts * flat_top
            trap = np.concatenate((ramp_up, flat, ramp_dn))

    else:
        trap, ramppts = 0, 0

    return np.expand_dims(trap, axis=0), ramppts


def spiral_varden(fov, res, gts, gslew, gamp, densamp, dentrans, nl, rewinder=False):
    r"""Variable density spiral designer. Produces trajectory, gradients,
    and slew rate. Gradient units returned are in g/cm, g/cm/s

    Args:
        fov (float): imaging field of view (cm).
        res (float): imaging isotropic resolution (cm).
        gts (float): gradient sample time in sec.
        gslew (float): max slew rate in g/cm/s.
        gamp (float): max gradient amplitude in g/cm.
        densamp (float):  duration of full density sampling (# of samples).
        dentrans (float): duration of transition from higher to lower
            (should be >= densamp/2).
        nl (float): degree of undersampling outer region.
        rewinder (Boolean): if True, include rewinder. If false, exclude.

    Returns:
        tuple: (g, k, t, s, dens) tuple containing

        - **g** - (array): gradient waveform [g/cm]
        - **k** - (array): exact k-space corresponding to gradient g.
        - **time** - (array):  sampled time
        - **s** - (array): slew rate [g/cm/s]
        - **dens** - (array): undersampling factor at each time point.

    References:
        Code and algorithm based on spiralgradlx6 from
        Doug Noll, U. of Michigan BME
    """
    gslew /= 100 # following actually assume that gslew is in T/m/s (or G/cm/cs)
    fsgcm = gamp  # fullscale g/cm
    risetime = gamp / gslew * 10000  # us
    ts = gts  # sampling time
    gts = gts  # gradient sampling time
    N = np.floor(fov / res)
    targetk = N / 2
    A = 32766  # output scaling of waveform (fullscale)

    max_dec_ratio = 32
    gam = 4257.0
    S = (gts / 1e-6) * A / risetime
    dr = ts / gts
    OMF = 2.0 * np.pi * fov / (1 / (gam * fsgcm * gts))
    OM = 2.0 * np.pi / nl * fov / (1 / (gam * fsgcm * gts))
    distance = 1.0 / (fov * gam * fsgcm * gts / A)

    ac = A
    loop = 1
    absk = 0
    dec_ratio = 1
    s0 = gslew * 100
    ggx, ggy = [], []
    dens = []
    kx, ky = [], []

    while loop > 0:
        loop = 0
        om = OM / dec_ratio
        omf = OMF / dec_ratio
        s = S / dec_ratio
        g0 = 0
        gx = g0
        gy = 0
        absg = np.abs(g0)
        oldkx = 0
        oldky = 0
        tkx = gx
        tky = gy
        kxt = tkx
        kyt = tky
        thetan_1 = 0
        taun = 0
        n = 0
        den1 = 0

        while absk < targetk:
            realn = n / dec_ratio
            taun_1 = taun
            taun = np.abs(tkx + 1j * tky) / A
            tauhat = taun
            if realn > densamp:
                if den1 == 0:
                    den1 = 1

                if realn > densamp + dentrans:
                    if "scthat" not in locals():
                        scthat = 0
                    scoffset = scthat
                    denoffset = taun_1
                    scthat = scoffset + om * (tauhat - denoffset)
                    fractrans = 1

                else:
                    scoffset = scthat
                    denoffset = taun_1
                    fractrans = (realn - densamp) / dentrans
                    fractrans = 1 - ((fractrans - 1) * (fractrans - 1))
                    scthat = omf + (om - omf) * fractrans
                    scthat *= tauhat - denoffset
                    scthat += scoffset

            else:
                fractrans = 0
                scthat = omf * tauhat

            theta = np.arctan2(scthat, 1.0) + scthat

            if absg < ac:
                deltheta = theta - thetan_1
                B = 1.0 / (1.0 + np.tan(deltheta) * np.tan(deltheta))
                gtilde = absg
                t1 = s * s
                t2 = gtilde * gtilde * (1 - B)

                if t2 > t1:
                    dec_ratio = dec_ratio * 2.0

                    if dec_ratio > max_dec_ratio:
                        print("k-space calculation failed.\n")
                        return

                    loop = 1
                    break

                t3 = np.sqrt(t1 - t2)
                absg = np.sqrt(B) * gtilde + t3

                if absg > ac:
                    absg = ac

            tgx = absg * np.cos(theta)
            tgy = absg * np.sin(theta)
            tkx += tgx
            tky += tgy
            thetan_1 = theta

            if np.remainder(n, dec_ratio) == 0:
                m = int(np.round(n / dec_ratio))
                gx = np.round((tkx - oldkx) / dec_ratio)
                gx = gx - np.remainder(gx, 2)
                gy = np.round((tky - oldky) / dec_ratio)
                gy = gy - np.remainder(gy, 2)
                if m > len(ggx) - 1:
                    ggx.append(gx)
                    ggy.append(gy)
                else:
                    ggx[m] = gx
                    ggy[m] = gy
                kxt = kxt + gx
                kyt = kyt + gy
                oldkx = tkx
                oldky = tky

                if np.remainder(m, dr) == 0:
                    m = int(m / dr)
                    absk = np.abs(kxt + 1j * kyt) / distance

                    if m > len(dens) - 1:
                        dens.append(omf / (omf + (om - omf) * fractrans))
                        if absk > targetk:
                            break
                        kx.append(kxt / distance)
                        ky.append(kyt / distance)
                    else:
                        dens[m] = omf / (omf + (om - omf) * fractrans)
                        if absk > targetk:
                            break
                        kx[m] = kxt / distance
                        ky[m] = kyt / distance

            n += 1

    g = []
    for i in range(len(ggx)):
        g.append(complex(ggx[i], ggy[i]) / A * fsgcm)
    dt = gts * 1000
    delk = 1 / 4.258 / fov  # (g ms)/cm

    # ramp down
    l2 = len(g) - 1
    rsteps = int(np.ceil(np.abs(g[l2]) / (s0 * 0.99) / gts))
    ind3 = l2 + np.linspace(1, rsteps, num=rsteps)
    c = g[l2] * np.linspace(rsteps, 0, num=rsteps) / rsteps
    g.extend(c)
    dens.extend([0] * len(ind3))

    # rewinder
    if rewinder:
        rewx, ramppts = trap_grad(abs(np.real(sum(g))) * gts, gamp, gslew * 50, gts)
        rewy, ramppts = trap_grad(abs(np.imag(sum(g))) * gts, gamp, gslew * 50, gts)
        
        rewx, rewy = rewx.squeeze(), rewy.squeeze()

        # append rewinder gradient
        if len(rewx) > len(rewy):
            r = -np.sign(np.real(sum(g))) * rewx
            p = np.sign(np.imag(sum(g)))
            p *= 1j * np.abs(np.imag(sum(g))) / np.real(sum(g)) * rewx
            r = r - p
        else:
            p = -np.sign(np.real(sum(g)))
            p *= np.abs(np.real(sum(g)) / np.imag(sum(g))) * rewy
            r = p - 1j * np.sign(np.imag(sum(g))) * rewy

        g = np.concatenate((g, r))

    # change from (real, imag) notation to (Nt, 2) notation
    gtemp = np.zeros((len(g), 2))
    gtemp[:, 0] = np.real(g)
    gtemp[:, 1] = np.imag(g)
    g = gtemp

    # calculate trajectory, slew rate factor from designed gradient
    k = np.cumsum(g, axis=0) * dt / delk / fov  # trajectory
    t = np.linspace(0, len(g), num=len(g) + 1)  # time vector
    s = np.diff(g, axis=0) / (gts * 1000)  # slew rate factor

    return g, k, t, s, dens

def cartesian(fov, n, dt, gamp, gslew, infos, dirx=-1, diry=1, dirz=1):
    r"""Basic cartesian single-line readout designer.

    Args:
        fov (float): imaging field of view in cm.
        n (int): # of pixels (square). resolution??
        etl (int): echo train length.
        dt (float): sample time in s.
        gamp (float): max gradient amplitude in mT/m.
        gslew (float): max slew rate in mT/m/ms.
        infos (Info): structure with sequence information
        dirx (int): x direction, -1 left to right, 1 right to left
        diry (int): y direction, -1 posterior-anterior, 1 anterior-posterior
        dirz (int): z direction, -1 feet-head, 1 head-feet

    Returns:
        tuple: (g, k, t, s) tuple containing

        - **g = [gx, gy]** - (array): gradient waveforms [mT/m] gx = [[[waveform array], [loops]], [[waveform array], [loops]], [[waveform array], [loops]], ...]
        - **k** - (array): exact k-space corresponding to gradient g.
        - **time** - (array):  sampled time
        - **s** - (array): slew rate [mT/m/ms]


    References:
        Adapted fron EPI version. 
    """
    is3D = infos.is3D
    
    s = gslew * dt * 1000
    scaley = 20 # original value 20


    # make the various gradient waveforms
    gamma = 4.2575  # kHz/Gauss
    g = (1 / (1000 * dt)) / (gamma * fov)  # Gauss/cm
    if g > gamp:
        g = gamp
        print("max g reduced to {}".format(g))

    # readout trapezoid
    adc = np.ones((1, n))
    gxro = g * np.ones((1, n))  # plateau of readout trapezoid
    areapd = np.sum(gxro) * dt

    ramp = np.expand_dims(np.linspace(0, g, int(np.ceil(g / s))+1), axis=0)
    gxro = np.concatenate((ramp, gxro, np.fliplr(ramp)),axis=1)

    # x prewinder. make sure res_kpre is even. Handle even N by changing prew.
    if n % 2 == 0:
        area = (np.sum(gxro) - dirx * g) * dt
    else:
        area = np.sum(gxro) * dt
    gxprew = dirx * trap_grad(area / 2, gamp, gslew * 1000, dt)[0]

    gxprew = np.concatenate(
        (np.zeros((1, (gxprew.size + ramp.size) % 2)), gxprew), axis=1
    )

    # phase-encode trapezoids before/after gx
    # handle even N by changing prewinder
    if n % 2 == 0:
        areayprew = areapd / 2 - g * dt
    else:
        areayprew = (areapd - g * dt) / 2 - g * dt

    gyprew = diry * trap_grad(areayprew, gamp, gslew  * 1000, dt)[0]
    gyprew = np.concatenate((np.zeros((1, gyprew.size % 2)), gyprew), axis=1)

    # add rephasers at end of gx and gy readout
    areagy = -areayprew  # units = G/cm*s
    gyrep = trap_grad(areagy, gamp, gslew * 1000, dt)[0]
    # gy = np.concatenate((gy, gyrep), axis=1)

    ##for 3D
    if is3D:
        # partition-encode trapezoids
        gslab = g * np.ones((1, infos.slices))  # test for 5 slabs
        areapdz = np.sum(gslab) * dt
        if infos.slices % 2 == 0:
            areazprew = areapdz / 2 - g * dt
        else:
            areazprew = (areapdz - g * dt) / 2 - g * dt
        gzprew = dirz * trap_grad(areazprew, gamp, gslew  * 1000, dt)[0]
        gzprew = np.concatenate((np.zeros((1, gzprew.size % 2)), gzprew), axis=1)
    ##for 3D
    
    areagx = area
    gxrep = trap_grad(-areagx, gamp, gslew * 1000, dt)[0]

    ## Prepare dynamic phase encoding
    sign = 1
    if areayprew<0:
        sign = -1
    gyprew_max_ampl = max(abs(gyprew[0]))
    step = gyprew_max_ampl / n
    gyprew_equation = str(sign) + "*(" + str(gyprew_max_ampl) + "-" + str(2*step) +"*counterPE)"

    sign = 1
    if areagy<0:
        sign = -1
    gyrep_max_ampl = max(abs(gyrep[0]))
    step = gyprew_max_ampl / n
    gyrep_equation = str(sign) + "*(" + str(gyrep_max_ampl) + "-" + str(2*step) +"*counterPE)"

    ##for 3D
    if is3D:
        # Prepare dynamic partition encoding
        signz = 1
        if areazprew<0:
            signz = -1
        gzprew_max_ampl = max(abs(gzprew[0]))
        stepz = gzprew_max_ampl / infos.slices # for part_num slabs
        gzprew_equation = str(signz) + "*(" + str(gzprew_max_ampl) + "-" + str(2*stepz) +"*counter3D)"
    ##for 3D
    
    # prepare blocks for mtrk
    gxprew_startTime = 0
    gyprew_startTime = 0
    gzprew_startTime = 0 ##for 3D
    gxro_startTime = max(gxprew.size, gyprew.size) * 10e-5
    adc1_startTime = gxro_startTime + (ramp.size + 1) * 10e-5 # 
    # gxrep_startTime = gxro_startTime + ( gxro.size * 10e-5 ) 
    gyrep_startTime = gxro_startTime + ( gxro.size * 10e-5 ) 

    block1 = [1, 
              [gxprew[0]/max(gxprew[0], key=abs), 
                gxprew.size,
                "read", 
                max(gxprew[0], key=abs), 
                gxprew_startTime],
              [gyprew[0]/max(gyprew[0], key=abs), 
               gyprew.size,
               "phase", 
               gyprew_equation, 
               gyprew_startTime],
              [gxro[0]/max(gxro[0], key=abs), 
               gxro.size,
               "read", 
               max(gxro[0], key=abs), 
               gxro_startTime],
              [adc[0],
               adc.size, 
               "adc", 
               1,
               adc1_startTime],
            #   [gxrep[0]/max(gxrep[0], key=abs), 
            #    gxrep.size,
            #    "read", 
            #    max(gxrep[0], key=abs), 
            #    gxrep_startTime],
              [gyrep[0]/max(gyrep[0], key=abs), 
               gyrep.size,
               "phase", 
               gyrep_equation, 
               gyrep_startTime]]
    if is3D:
        block1.append([gzprew[0]/max(gzprew[0], key=abs),  ##for 3D
               gzprew.size,
               "slice", 
               gzprew_equation, 
               gzprew_startTime])
    
    blocks = [block1]

    time_before_center = (adc1_startTime * 1e2 + (adc.size/2) *dt * 1e3)
    time_after_center = gyrep_startTime * 1e2 + (gyrep.size)*dt * 1e3 - time_before_center

    return blocks, time_before_center, time_after_center

def radial(fov, n_spokes, theta, dt, gamp, gslew):
    r"""Basic radial single-line readout designer.

    Args:
        fov (float): imaging field of view in cm.
        n (int): # of pixels (square). resolution??
        etl (int): echo train length.
        dt (float): sample time in s.
        gamp (float): max gradient amplitude in mT/m.
        gslew (float): max slew rate in mT/m/ms.
        offset (int): used for multi-shot EPI goes from 0 to #shots-1
        dirx (int): x direction of EPI -1 left to right, 1 right to left
        diry (int): y direction of EPI -1 bottom-top, 1 top-bottom

    Returns:
        tuple: (g, k, t, s) tuple containing

        - **g = [gx, gy]** - (array): gradient waveforms [mT/m] gx = [[[waveform array], [loops]], [[waveform array], [loops]], [[waveform array], [loops]], ...]
        - **k** - (array): exact k-space corresponding to gradient g.
        - **time** - (array):  sampled time
        - **s** - (array): slew rate [mT/m/ms]


    References:
        Adapted fron EPI version. 
    """
    s = gslew * dt * 1000
    scaley = 20 # original value 20


    # make the various gradient waveforms
    gamma = 4.2575  # kHz/Gauss
    g = (1 / (1000 * dt)) / (gamma * fov)  # Gauss/cm
    if g > gamp:
        g = gamp
        print("max g reduced to {}".format(g))

    # readout trapezoid
    adc = np.ones((1, n_spokes))
    ro = g * np.ones((1, n_spokes))  # plateau of readout trapezoid
    areapd = np.sum(ro) * dt

    ramp = np.expand_dims(np.linspace(s, g, int(g / s)), axis=0)
    ro = np.concatenate(
        (np.expand_dims(np.array([0]), axis=1), ramp, ro, np.fliplr(ramp)),
        axis=1,
    )

    # prewinder. make sure res_kpre is even. Handle even N by changing prew.
    if n_spokes % 2 == 0:
        area = (np.sum(ro) - g) * dt
    else:
        area = np.sum(ro) * dt
    prew = trap_grad(area / 2, gamp, gslew * 1000, dt)[0]

    prew = np.concatenate(
        (np.zeros((1, (prew.size + ramp.size) % 2)), prew), axis=1
    )


    # add rephasers at end of gx and gy readout
    areagx = area
    gxrep = trap_grad(-areagx, gamp, gslew * 1000, dt)[0]

    ## Prepare dynamic encoding
    sign = 1
    if areapd<0:
        sign = -1
    ro_max_ampl = max(abs(ro[0]))
    gxro_equation = str(sign) + "*" + str(ro_max_ampl) + "*cos(counter3*" + str(theta) +")"
    gyro_equation = str(sign) + "*" + str(ro_max_ampl) + "*sin(counter3*" + str(theta) +")"


    prew_max_ampl = max(abs(prew[0]))
    gxprew_equation = str(-sign) + "*" + str(prew_max_ampl) + "*cos(counter3*" + str(theta) +")"
    gyprew_equation = str(-sign) + "*" + str(prew_max_ampl) + "*sin(counter3*" + str(theta) +")"
    
    # prepare blocks for mtrk
    gxprew_startTime = 0
    gyprew_startTime = 0
    gxro_startTime = (prew.size+1) * 10e-5
    gyro_startTime = (prew.size+1) * 10e-5
    adc1_startTime = gxro_startTime + (ramp.size + 1) * 10e-5 # 
    gxrew_startTime = gyro_startTime + (ro.size+1) * 10e-5
    gyrew_startTime = gyro_startTime + (ro.size+1) * 10e-5
    # gxrep_startTime = gxro_startTime + ( gxro.size * 10e-5 ) 

    block1 = [1, 
              [prew[0]/max(prew[0], key=abs), 
                prew.size,
                "read", 
                gxprew_equation, 
                gxprew_startTime],
              [prew[0]/max(prew[0], key=abs), 
               prew.size,
               "phase", 
               gyprew_equation, 
               gyprew_startTime],
              [ro[0]/max(ro[0], key=abs), 
               ro.size,
               "read", 
               gxro_equation, 
               gxro_startTime],
              [ro[0]/max(ro[0], key=abs), 
               ro.size,
               "phase", 
               gyro_equation, 
               gyro_startTime],
              [adc[0],
               adc.size, 
               "adc", 
               1,
               adc1_startTime],
              [prew[0]/max(prew[0], key=abs), 
                prew.size,
                "read", 
                gxprew_equation, 
                gxrew_startTime],
              [prew[0]/max(prew[0], key=abs), 
               prew.size,
               "phase", 
               gyprew_equation, 
               gyrew_startTime]
              ]  
    
    blocks = [block1]

    time_before_center = (gxro_startTime + (ro.size/2) * dt) * 1e2

    return blocks, time_before_center

def spiral_arch(fov, resolution, gts, gslew, gamp):
    r"""Analytic Archimedean spiral designer. Produces trajectory, gradients,
    and slew rate. Gradient returned has units mT/m.

    Args:
        fov (float): imaging field of view in m.
        res (float): resolution, in m.
        gts (float): sample time in s.
        gslew (float): max slew rate in mT/m/ms.
        gamp (float): max gradient amplitude in mT/m.

    Returns:
        tuple: (g, k, t, s) tuple containing

        - **g** - (array): gradient waveform [mT/m]
        - **k** - (array): exact k-space corresponding to gradient g.
        - **time** - (array):  sampled time
        - **s** - (array): slew rate [mT/m/ms]

    References:
        Glover, G. H.(1999).
        Simple Analytic Spiral K-Space Algorithm.
        Magnetic resonance in medicine, 42, 412-415.

        Bernstein, M.A.; King, K.F.; amd Zhou, X.J. (2004).
        Handbook of MRI Pulse Sequences. Elsevier.
    """

    print("resolution: ", resolution, "cm")
    gslew = gslew * 1e2 ## conversion to G/cm/s
    gambar = 4257 # Hz/Gauss
    gam = gambar * (2*np.pi) 
    lam = 1 / (2 * np.pi * fov)
    beta = gambar * gslew / lam

    a_2 = (9 * beta / 4) ** (1 / 3)  # rad ** (1/3) / s ** (2/3)
    lamb = 5
    theta_max = (np.pi * resolution)
    ts = (3 * gambar * gamp / (2 * lam * a_2**2)) ** 3
    theta_s = 0.5 * beta * ts**2
    theta_s /= lamb + beta / (2 * a_2) * ts ** (4 / 3)
    t_g = np.pi * lam * (theta_max**2 - theta_s**2) / (gam * gamp)
    n_s = int(np.round(ts / gts))
    n_g = int(np.round(t_g / gts))
    
    if theta_max > theta_s:
        print(" Spiral trajectory is slewrate limited or amplitude limited")

        tacq = ts + t_g

        t_s = np.linspace(0, ts, n_s)
        t_g = np.linspace(ts + gts, tacq, n_g)

        theta_1 = beta / 2 * t_s**2
        theta_1 /= lamb + beta / (2 * a_2) * t_s ** (4 / 3)
        theta_2 = theta_s**2 + gam / (np.pi * lam) * gamp * (t_g - ts)
        theta_2 = np.sqrt(theta_2)

        k1 = lam * theta_1 * (np.cos(theta_1) + 1j * np.sin(theta_1))
        k2 = lam * theta_2 * (np.cos(theta_2) + 1j * np.sin(theta_2))
        k = np.concatenate((k1, k2), axis=0)

    else:
        tacq = (2*np.pi * fov/3) * np.sqrt(1/(2*gambar*gslew*(fov/resolution)**3))
        n_t = int(np.round(tacq / gts))
        t_s = np.linspace(0, tacq, n_t)
        theta_1 = beta / 2 * t_s**2
        theta_1 /= lamb + beta / (2 * a_2) * t_s ** (4 / 3)

        k = lam * theta_1 * (np.cos(theta_1) + 1j * np.sin(theta_1)) 

    # end of trajectory calculation; prepare outputs
    g = np.diff(k, 1, axis=0) / (gts * gambar)  # gradient
    g = np.pad(g, (0, 1), "constant")
    s = np.diff(g, 1, axis=0) / (gts * 1000)  # slew rate factor
    s = np.pad(s, (0, 1), "constant")

    # change from (real, imag) notation to (Nt, 2) notation
    k = traj_complex_to_array(k)
    g = traj_complex_to_array(g)
    s = traj_complex_to_array(s)

    t = np.linspace(0, len(g), num=len(g) + 1)  # time vector

    spiral_array = np.transpose(g)
    gx = spiral_array[0] #[::10]
    gy = spiral_array[1] #[::10]

    print("gx last: ", np.abs(gx[-2] - gx[-1])/gts, "G/cm")
    print("gy last: ", np.abs(gy[-2] - gy[-1])/gts, "G/cm")
    print(gslew)
    # if np.abs(gy[-2] - gy[-1])/gts> gslew*1e-1:
    #     print("Warning: gy last gradient is larger than gslew")
    #     print("adding a ramp to gy")
    #     print(int(np.round((gy[-2] - gy[-1])/(gslew*1e-1*gts))))
    #     # removing the last point
    #     gy = gy[:-1]
    #     # adding a ramp
    #     ramp = np.linspace(gy[-1], 0, int(np.round(gy[-1]/(gslew*1e-1*gts))+1))
    #     print("ramp ", ramp/max(gy, key=abs))
    #     gy = np.concatenate((gy, ramp))
    #     print("new gy last: ", np.abs(gy[-2] - gy[-1])/gts, "G/cm")

    gx_startTime = 0
    gy_startTime = 0
    adc_startTime = 0
    blocks = [[1,
              [gx/max(gx, key=abs), 
               gx.size,
               "read", 
               max(gx, key=abs), 
               gx_startTime],
               [gy/max(gy, key=abs), 
               gy.size,
               "phase", 
               max(gy, key=abs), 
               gy_startTime],
              [np.ones(resolution*resolution),
               gx.size, 
               "adc", 
               1,
               adc_startTime]]]
    
    # duration = gx.size * gts * 1e3 # different from epi, error?
    # print("Duration: ", duration, "ms")

    # time_before_center = (gx.size/2) * gts * 1e3
    time_before_center = 0.0

    return blocks, time_before_center, g, k, t, s


def spiral_k(fov, N, f_sampling, R, ninterleaves, alpha, gm, sm, gamma=2.678e8):
    """Generate variable density spiral trajectory ONLY. No gradient included.

    Args:
        fov (float): field of view in meters.
        N (int): effective matrix shape.
        f_sampling (float): undersampling factor in freq encoding direction.
        R (float): undersampling factor.
        ninterleaves (int): number of spiral interleaves
        alpha (float): variable density factor
        gm (float): maximum gradient amplitude (T/m)
        sm (float): maximum slew rate (T/m/s)
        gamma (float): gyromagnetic ratio in rad/T/s

    Returns:
        array: spiral coordinates.

    References:
        Dong-hyun Kim, Elfar Adalsteinsson, and Daniel M. Spielman.
        'Simple Analytic Variable Density Spiral Design.' MRM 2003.

    """
    res = fov / N

    lam = 0.5 / res  # in m**(-1)
    n = 1 / (1 - (1 - ninterleaves * R / fov / lam) ** (1 / alpha))
    w = 2 * np.pi * n
    Tea = lam * w / gamma / gm / (alpha + 1)  # in s
    Tes = np.sqrt(lam * w**2 / sm / gamma) / (alpha / 2 + 1)  # in s
    Ts2a = (
        Tes ** ((alpha + 1) / (alpha / 2 + 1)) * (alpha / 2 + 1) / Tea / (alpha + 1)
    ) ** (
        1 + 2 / alpha
    )  # in s

    if Ts2a < Tes:
        tautrans = (Ts2a / Tes) ** (1 / (alpha / 2 + 1))

        def tau(t):
            return (t / Tes) ** (1 / (alpha / 2 + 1)) * (0 <= t) * (t <= Ts2a) + (
                (t - Ts2a) / Tea + tautrans ** (alpha + 1)
            ) ** (1 / (alpha + 1)) * (t > Ts2a) * (t <= Tea) * (Tes >= Ts2a)

        Tend = Tea
    else:

        def tau(t):
            return (t / Tes) ** (1 / (alpha / 2 + 1)) * (0 <= t) * (t <= Tes)

        Tend = Tes

    def k(t):
        return lam * tau(t) ** alpha * np.exp(w * tau(t) * 1j)

    dt = Tea * 1e-4  # in s

    Dt = dt * f_sampling / fov / abs(k(Tea) - k(Tea - dt))  # in s

    t = np.linspace(0, Tend, int(Tend / Dt))
    kt = k(t)  # in rad

    # generating cloned interleaves
    k = kt
    for i in range(1, ninterleaves):
        k = np.hstack((k, kt[0:] * np.exp(2 * np.pi * 1j * i / ninterleaves)))

    k = np.stack((np.real(k), np.imag(k)), axis=1)

    return k


def mtrk_epi(fov, n, etl, dt, gamp, gslew, offset=0, dirx=-1, diry=1):
    r"""Basic EPI single-shot trajectory designer.

    Args:
        fov (float): imaging field of view in cm.
        n (int): # of pixels (square). N = etl*nl, where etl = echo-train-len
            and nl = # leaves (shots). nl default 1.
        etl (int): echo train length.
        dt (float): sample time in s.
        gamp (float): max gradient amplitude in mT/m.
        gslew (float): max slew rate in mT/m/ms.
        offset (int): used for multi-shot EPI goes from 0 to #shots-1
        dirx (int): x direction of EPI -1 left to right, 1 right to left
        diry (int): y direction of EPI -1 bottom-top, 1 top-bottom

    Returns:
        tuple: (g, k, t, s) tuple containing

        - **g = [gx, gy]** - (array): gradient waveforms [mT/m] gx = [[[waveform array], [loops]], [[waveform array], [loops]], [[waveform array], [loops]], ...]
        - **k** - (array): exact k-space corresponding to gradient g.
        - **time** - (array):  sampled time
        - **s** - (array): slew rate [mT/m/ms]


    References:
        From Antonis Matakos' contrib to Jeff Fessler's IRT.
    """
    s = gslew * dt * 1000 
    scaley = 20 # original value 20


    # make the various gradient waveforms
    gamma = 4.2575  # kHz/Gauss
    g = (1 / (1000 * dt)) / (gamma * fov)  # Gauss/cm
    if g > gamp:
        g = gamp
        print("max g reduced to {}".format(g))

    # readout trapezoid
    adc = np.ones((1, n))
    # adc = np.ones((1, n+1)) # added 1 to adc to make it even???
    gxro = g * np.ones((1, n))  # plateau of readout trapezoid
    areapd = np.sum(gxro) * dt

    ramp = np.expand_dims(np.linspace(s, g, int(g / s)), axis=0)
    gxro = np.concatenate(
        (np.expand_dims(np.array([0]), axis=1), ramp, gxro, np.fliplr(ramp)),
        axis=1,
    )

    # x prewinder. make sure res_kpre is even. Handle even N by changing prew.
    if n % 2 == 0:
        area = (np.sum(gxro) - dirx * g) * dt
    else:
        area = np.sum(gxro) * dt
    gxprew = dirx * trap_grad(area / 2, gamp, gslew * 1000, dt)[0]

    gxprew = np.concatenate(
        (np.zeros((1, (gxprew.size + ramp.size) % 2)), gxprew), axis=1
    )

    # partial dephaser (one cycle of phase across each voxel)
    gxpd = -trap_grad(areapd / 2, gamp, gslew * 1000, dt)[0]
    gxpd = np.concatenate((np.zeros((1, gxpd.size % 2)), gxpd), axis=1)

    # phase-encode trapezoids before/after gx
    # handle even N by changing prewinder
    if n % 2 == 0:
        areayprew = areapd / 2 - offset * g * dt
    else:
        areayprew = (areapd - g * dt) / 2 - offset * g * dt

    gyprew = diry * trap_grad(areayprew, gamp, gslew / scaley * 1000, dt)[0]
    gyprew = np.concatenate((np.zeros((1, gyprew.size % 2)), gyprew), axis=1)

    lx = gxpd.size
    ly = gyprew.size
    if lx > ly:
        gyprew = np.concatenate((gyprew, np.zeros((1, lx - ly))), axis=1)
    else:
        gxpd = np.concatenate((gxpd, np.zeros((1, ly - lx))), axis=1)

    # gy readout gradient elements
    # changed readout patterns to create interleaved EPIs
    areagyblip = areapd / etl
    gyblip = trap_grad(areagyblip, gamp, gslew / scaley * 1000, dt)[0]
    gyro = np.concatenate((np.zeros((1, gxro.size - gyblip.size)), gyblip), axis=1)
    gyro2 = np.expand_dims(np.array([0]), axis=1)

    # put together gx and gy

    gxro = -dirx * gxro
    gx = gxprew

    gyro = -diry * gyro
    gyro2 = -diry * gyro2
    gy = np.expand_dims(np.array([0]), axis=1)
    lx = gx.size
    ly = gy.size
    if lx > ly:
        gy = np.concatenate((gy, np.zeros((1, lx - ly))), axis=1)
    else:
        gx = np.concatenate((gx, np.zeros((1, ly - lx))), axis=1)

    gy = np.concatenate((gy, np.zeros((1, int(gyblip.size / 2)))), axis=1)

    for ee in range(1, etl):
        flip = (-1) ** (ee + 1)
        gx = np.concatenate((gx, flip * gxro), axis=1)
        gy = np.concatenate((gy, gyro), axis=1)

    if etl == 1:
        ee = 1
    else:
        ee += 1

    # concatenate with added 0 to limit max s
    gx = np.concatenate(
        (gx, (-(1 ** (ee + 1)) * gxro), np.expand_dims(np.array([0]), axis=1)),
        axis=1,
    )
    gy = np.concatenate((gy, np.zeros((1, gx.size - gy.size))), axis=1)

    # add rephasers at end of gx and gy readout
    areagx = np.sum(gx) * dt
    gxrep = trap_grad(-areagx, gamp, gslew * 1000, dt)[0]
    gx = np.concatenate((gx, gxrep), axis=1)

    areagy = np.sum(gy) * dt  # units = G/cm*s
    gyrep = trap_grad(-areagy, gamp, gslew / scaley * 1000, dt)[0]
    gy = np.concatenate((gy, gyrep), axis=1)

    # make sure length of gx and gy are same, and even
    lx = gx.size
    ly = gy.size
    if lx > ly:
        gy = np.concatenate((gy, np.zeros((1, lx - ly))), axis=1)
    else:
        gx = np.concatenate((gx, np.zeros((1, ly - lx))), axis=1)


    
    # prepare blocks for mtrk
    gxprew_startTime = 0
    gyprew_startTime = 0
    gxro1_startTime = 0
    adc1_startTime = (ramp.size + 1) * 10e-5 # 
    gyblip1_startTime = gxro.size * 10e-5 # size - 1 ???
    gxro2_startTime = gyblip1_startTime + ( gyblip.size * 10e-5 ) # size - 1 ???
    adc2_startTime = gxro2_startTime + ( (ramp.size + 1) * 10e-5 ) # size - 1 ???
    gyblip2_startTime = gxro2_startTime + ( gxro.size * 10e-5 ) # size - 1 ???
    gxrep_startTime = 0
    gyrep_startTime = 0

    block1 = [1, 
              [gxprew[0]/max(gxprew[0], key=abs), 
                gxprew.size,
                "read", 
                max(gxprew[0], key=abs), 
                gxprew_startTime],
              [gyprew[0]/max(gyprew[0], key=abs), 
               gyprew.size,
               "phase", 
               max(gyprew[0], key=abs), 
               gyprew_startTime]]
    block1_duration = gxprew_startTime + max(gxprew.size, gyprew.size) * dt 
    # print("block1_duration: ", block1_duration*1e3)

    block2 = [etl/2-1,
              [gxro[0]/max(gxro[0], key=abs), 
               gxro.size,
               "read", 
               max(gxro[0], key=abs), 
               gxro1_startTime],
              [adc[0],
               adc.size, 
               "adc", 
               1,
               adc1_startTime],
              [gyblip[0]/max(gyblip[0], key=abs), 
               gyblip.size,
               "phase", 
               max(gyblip[0], key=abs), 
               gyblip1_startTime],
              [gxro[0]/max(gxro[0], key=abs), 
               gxro.size,
               "read", 
               -max(gxro[0], key=abs), 
               gxro2_startTime],
              [adc[0],
               adc.size,  
               "adc", 
               1,
               adc2_startTime],
              [gyblip[0]/max(gyblip[0], key=abs), 
               gyblip.size,
               "phase", 
               max(gyblip[0], key=abs), 
               gyblip2_startTime]]
    block2_duration = gyblip2_startTime + gyblip.size * dt
    # print("block2_duration: ", block2_duration*1e3)

    block3 = [1,
              [gxro[0]/max(gxro[0], key=abs), 
               gxro.size,
               "read", 
               max(gxro[0], key=abs), 
               gxro1_startTime],
              [adc[0],  
               adc.size, 
               "adc",
               1,
               adc1_startTime],
              [gyblip[0]/max(gyblip[0], key=abs), 
               gyblip.size,
               "phase", 
               max(gyblip[0], key=abs), 
               gyblip1_startTime],
              [gxro[0]/max(gxro[0], key=abs), 
               gxro.size,
               "read", 
               -max(gxro[0], key=abs), 
               gxro2_startTime],
              [adc[0], 
               adc.size,
               "adc", 
               1,
               adc2_startTime]]
    block3_duration = gxro2_startTime + gxro.size * dt 
    # print("block3_duration: ", block3_duration*1e3)

    block4 = [1,
              [gxrep[0]/max(gxrep[0], key=abs), 
               gxrep.size,
               "read", 
               max(gxrep[0], key=abs), 
               gxrep_startTime],
              [gyrep[0]/max(gyrep[0], key=abs), 
               gyrep.size,
               "phase", 
               max(gyrep[0], key=abs), 
               gyrep_startTime]]
    block4_duration = gxrep_startTime + max(gxrep.size, gyrep.size) * dt
    # print("block4_duration: ", block4_duration*1e3)
    
    # returns **blocks** with the following structure: 
    # blocks = [block1, block2, ..., blockN]
    # Each block is a list of lists with the following structure:
    # block = [[waveform1, duration1, axis1, amplitude1, startTime1], ..., [waveformM, durationM, axisM, amplitudeM, startTimeM], loops]
    
    blocks = [block1, block2, block3, block4]

    duration = (block1_duration + (etl/2-1) * block2_duration + block3_duration + block4_duration) * 1e2
    # print("Duration: ", duration, "ms")

    time_before_center =  (block1_duration + (etl/4) * block2_duration) * 1e2

    # subplot, axis = plt.subplots(2, sharex=True)
    # subplot.suptitle("EPI trajectory")
    # axis[0].set_title("gy")
    # axis[0].plot(gyblip[0])
    # axis[1].set_title("gx")
    # axis[1].plot(gxro[0])
    # plt.show()    

    gx = np.concatenate((gx, np.zeros((1, gx.size % 2))), axis=1)
    gy = np.concatenate((gy, np.zeros((1, gy.size % 2))), axis=1)
    g = np.concatenate((gx, gy), axis=0)

    sx = np.diff(gx, axis=1) / (dt * 1000)
    sy = np.diff(gy, axis=1) / (dt * 1000)
    s = np.concatenate((sx, sy), axis=0)

    kx = np.cumsum(gx, axis=1) * gamma * dt * 1000
    ky = np.cumsum(gy, axis=1) * gamma * dt * 1000
    k = np.concatenate((kx, ky), axis=0)

    t = np.linspace(0, kx.size, kx.size) * dt

    return blocks, time_before_center


def rosette(kmax, w1, w2, dt, dur, gamp=None, gslew=None):
    r"""Basic rosette trajectory designer.

    Args:
        kmax (float): 1/m.
        w1 (float): rotational frequency (Hz).
        w2 (float): center sampling frequency (Hz).
        dt (float): sample time (s).
        dur (float): total duration (s).
        gamp (float): max gradient amplitude (mT/m).
        gslew (float): max slew rate (mT/m/ms).

    Returns:
        tuple: (g, k, t, s) tuple containing

        - **g** - (array): gradient waveform [mT/m]
        - **k** - (array): exact k-space corresponding to gradient g.
        - **time** - (array):  sampled time
        - **s** - (array): slew rate [mT/m/ms]

    References:
        D. C. Noll, 'Multi-shot rosette trajectories for spectrally selective
        MR imaging.' IEEE Trans. Med Imaging 16, 372-377 (1997).
    """

    # check if violates gradient or slew rate constraints
    gam = 267.522 * 1e6 / 1000  # rad/s/mT
    gambar = gam / 2 / np.pi  # Hz/mT
    if gamp is not None:
        if (1 / gambar) * kmax * w1 > gamp:
            print("gmax exceeded, decrease rosette kmax or w1")
            return
    if gslew is not None:
        if (1 / gambar) * kmax * (w1**2 + w2**2) / 1000 > gslew:
            print("smax exceeded, dcrease rosette kmax, w1, or w2")
            return
    t = np.linspace(0, dur, dur / dt)
    k = kmax * np.sin(w1 * t) * np.exp(1j * w2 * t)

    # end of trajectory calculation; prepare outputs
    g = np.diff(k, 1, axis=0) / (dt * gambar)  # gradient
    g = np.pad(g, (0, 1), "constant")
    s = np.diff(g, 1, axis=0) / (dt * 1000)  # slew rate factor
    s = np.pad(s, (0, 1), "constant")

    # change from (real, imag) notation to (Nt, 2) notation
    k = traj_complex_to_array(k)
    g = traj_complex_to_array(g)
    s = traj_complex_to_array(s)

    t = np.linspace(0, len(g), num=len(g) + 1)  # time vector

    return g, k, t, s


def spokes_grad(k, tbw, sl_thick, gmax, dgdtmax, gts):
    r"""Radial trajectory gradient designer.

    Args:
        k (array): spokes locations, [Nspokes, 2]
        tbw (int): time bandwidth product.
        sl_thick (float): slice thickness (mm).
        gmax (float): max gradient amplitude (g/cm).
        dgdtmax (float): max gradient slew (g/cm/s).
        gts (float): hardware sampling dwell time (s).

    Returns:
        g (array): gz, gy, and gz waveforms  in g/cm [3, Nt]

    References:
          Pypulseq (examples/scripts.write_radial_gre) 
    """
    fov = 260e-3
    Nx = 64  # Define FOV and resolution
    alpha = 10  # Flip angle
    slice_thickness = 3e-3  # Slice thickness
    TE = 8e-3  # Echo time
    TR = 20e-3  # Repetition time
    Nr = 60  # Number of radial spokes
    N_dummy = 20  # Number of dummy scans
    delta = np.pi / Nr  # Angular increment
    plot = False
    write_seq = False
    seq_filename = "radial_gre.seq"
    rf_spoiling_inc = 117  # RF spoiling increment

    ## System limits
    max_grad=28,
    grad_unit='mT/m',
    max_slew=120,
    slew_unit='T/m/s',
    rf_ringdown_time=20e-6,
    rf_dead_time=100e-6,
    adc_dead_time=10e-6
    gradient_raster_time=10e-6

    ## Define gradients and ADC events
    deltak = 1 / fov

    gx_area = Nx * deltak

    gx = min_trap_grad(gx_area, max_grad, max_slew, gradient_raster_time)


    # return g


def spokes_grad(k, tbw, sl_thick, gmax, dgdtmax, gts):
    r"""Spokes gradient designer. Given some chosen spoke locations k, return
    the gradients required to move between those spoke locations.

    Args:
        k (array): spokes locations, [Nspokes, 2]
        tbw (int): time bandwidth product.
        sl_thick (float): slice thickness (mm).
        gmax (float): max gradient amplitude (g/cm).
        dgdtmax (float): max gradient slew (g/cm/s).
        gts (float): hardware sampling dwell time (s).

    Returns:
        g (array): gz, gy, and gz waveforms  in g/cm [3, Nt]

    References:
           Grissom, W., Khalighi, M., Sacolick, L., Rutt, B. & Vogel, M (2012).
           Small-tip-angle spokes pulse design using interleaved greedy and
           local optimization methods. Magnetic Resonance in Medicine, 68(5),
           1553-62.

    """
    n_spokes = k.shape[0]

    area = tbw / (sl_thick / 10) / 4257  # thick * kwid = twb, kwid = gam*area
    [subgz, nramp] = min_trap_grad(area, gmax, dgdtmax, gts)

    # calc gradient, add extra 0 location at end for return to (0, 0)
    # print("k[:, 0], k[:, 1] ", k[:, 0], k[:, 1])
    gxarea = np.diff(np.concatenate((k[:, 0], np.zeros(1)))) / 4257
    gyarea = np.diff(np.concatenate((k[:, 1], np.zeros(1)))) / 4257

    gx, gy, gz = [], [], []
    gz_sign = -1
    for ii in range(n_spokes):
        # print("spoke ", ii + 1, " of ", n_spokes)
        gz_sign *= -1
        gz.extend(np.squeeze(gz_sign * subgz).tolist())  # alt sign of gz

        gx.extend([0] * np.size(subgz))  # zeros for gz duration
        if np.absolute(gxarea[ii]) > 0:
            # print("gxarea[ii] ", gxarea[ii])
            [gblip, _] = trap_grad(abs(gxarea[ii]), gmax, dgdtmax, gts)
            gxblip = int(np.sign(gxarea[ii])) * gblip
            gx = gx[: len(gx) - len(gxblip.T)]
            gx.extend(np.squeeze(gxblip).tolist())

        gy.extend([0] * np.size(subgz))
        if np.absolute(gyarea[ii]) > 0:
            # print("gyarea[ii] ", gyarea[ii])
            [gblip, _] = trap_grad(abs(gyarea[ii]), gmax, dgdtmax, gts)
            gyblip = int(np.sign(gyarea[ii])) * gblip
            gy = gy[: len(gy) - len(gyblip.T)]
            gy.extend(np.squeeze(gyblip).tolist())

    [gref, _] = trap_grad(gts * np.sum(subgz) / 2, gmax, dgdtmax, gts)
    gzref = -gref
    gz.extend(np.squeeze(gzref).tolist())
    gx.extend([0] * np.size(gzref))
    gy.extend([0] * np.size(gzref))

    # combine gradient waveforms
    gx = np.array(gx)
    g = np.vstack((np.array(gx), np.array(gy), np.array(gz)))    

    return g


def stack_of(k, num, zres):
    r"""Function for creating a 3D stack of ____ trajectory from a 2D [Nt 2]
    trajectory.

    Args:
        k (array): 2D array in [2 x Nt]. Will be bottom of stack.
        num (int): number of layers of stack.
        zres (float): spacing between stacks in cm.
    """

    z = np.linspace(-num * zres / 2, num * zres / 2, num)
    kout = np.zeros((k.shape[0] * num, 3))

    # we will be performing a complex rotation on our trajectory
    k = traj_array_to_complex(k)

    for ii in range(num):
        kr = k[0:] * np.exp(2 * np.pi * 1j * ii / num)
        z_coord = np.expand_dims(np.ones(len(kr)) * z[ii], axis=1)
        krz = np.concatenate((traj_complex_to_array(kr), z_coord), axis=1)

        kout[ii * len(krz) : (ii + 1) * len(krz), :] = krz

    return kout


def traj_complex_to_array(k):
    r"""Function to convert complex convention trajectory to [Nt 2] trajectory

    Args:
        k (complex array): Nt vector
    """
    kout = np.zeros((len(k), 2))
    kout[:, 0], kout[:, 1] = np.real(k), np.imag(k)
    return kout


def traj_array_to_complex(k):
    r"""Function to convert [Nt 2] convention traj to complex convention

    Args:
        k (complex array): Nt vector
    """
    kout = k[:, 0] + 1j * k[:, 1]
    return kout